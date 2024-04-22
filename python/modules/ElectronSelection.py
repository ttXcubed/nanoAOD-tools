import json
import math
import os
import random
import sys
from collections import OrderedDict

import numpy as np
import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, getHist, getSFXY


class ElectronSelection(Module):
    WP80 = 1
    WP90 = 2
    INV = 3
    NONE = 4

    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Electron"),
        #outputName_list = ["tightElectrons","mediumElectrons","looseElectrons"],
        triggerMatch=False,
        id_type=[],
        electronMinPt = 29.,
        electronMaxEta = 2.4,
        storeKinematics= [] ,#['pt','eta'],
        storeTruthKeys=[]
    ):

        self.inputCollection = inputCollection
        self.id_type = id_type
        #self.outputName_list = outputName_list
        self.electronMinPt = electronMinPt
        self.electronMaxEta = electronMaxEta
        self.storeKinematics = storeKinematics
        self.triggerMatch = triggerMatch
        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []
        self.storeTruthKeys = storeTruthKeys

        self.WPs = ['tight', 'medium', 'loose']

        self.outputName_dict = OrderedDict()
        for id in id_type:
            self.outputName_dict[id] = OrderedDict()
            for wp in self.WPs:
                self.outputName_dict[id][wp] = wp+'_'+id+'_'+'Electrons'

        self.mva_id_wp = {
            'tight': {
                '2016': 'mvaFall17V2Iso_WP80', '2016preVFP': 'mvaFall17V2Iso_WP80', 
                '2017': 'mvaFall17V2Iso_WP80', '2018': 'mvaFall17V2Iso_WP80',
                '2022': 'mvaIso_WP80', '2022EE': 'mvaIso_WP80'
            },
            'medium': {
                '2016': 'mvaFall17V2Iso_WP90', '2016preVFP': 'mvaFall17V2Iso_WP90', 
                '2017': 'mvaFall17V2Iso_WP90', '2018': 'mvaFall17V2Iso_WP90',
                '2022': 'mvaIso_WP90', '2022EE': 'mvaIso_WP90' 
            },
            'loose': {
                '2016': 'mvaFall17V2Iso_WPL', '2016preVFP': 'mvaFall17V2Iso_WPL', 
                '2017': 'mvaFall17V2Iso_WPL', '2018': 'mvaFall17V2Iso_WPL',
                '2022': '', '2022EE': '' 
            }
            # no MVA ID loose WP for Nano v12 --> Run 3 
            #https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv12/2022/2023/doc_DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2.html#Electron
        }


    def triggerMatched(self, electron, trigger_object):
        #return 2 arguments: 
        #   -1st: flag-->True, if electron matches a trigger object (deltaR<0.3); False, otherwise
        #   -2nd: idx of the matched trigger object
        min_deltaR, matchedTrgObj_id = 1., None
        if self.triggerMatch:
            #trig_deltaR = math.pi
            for itrgObj, trig_obj in enumerate(trigger_object):
                if abs(trig_obj.id) != 11:
                    continue
                deltaR_trgObj_lep = deltaR(electron,trig_obj)
                if deltaR_trgObj_lep < min_deltaR:
                    min_deltaR = deltaR_trgObj_lep
                    matchedTrgObj_id = itrgObj
            if min_deltaR < 0.3:
                return True, matchedTrgObj_id
            else:
                return False, None
        else:
            return True, None 

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch('nElectron', 'I')
        self.out.branch("electron_cutBasedID","F",lenVar="nElectron")
        for wp in ['tight', 'medium', 'loose']:
            self.out.branch("electron_MVA_"+wp+"ID","F",lenVar="nElectron")

        for id_type in self.id_type:
            for wp in self.outputName_dict[id_type].keys():
                self.out.branch("n"+self.outputName_dict[id_type][wp], "I")

                for variable in self.storeKinematics:
                    self.out.branch(self.outputName_dict[id_type][wp]+"_"+variable,"F",lenVar="n"+self.outputName_dict[id_type][wp])
                if not Module.globalOptions["isData"]:
                    for variable in self.storeTruthKeys:
                        self.out.branch(self.outputName_dict[id_type][wp]+"_"+variable,"F",lenVar="n"+self.outputName_dict[id_type][wp])
                self.out.branch(self.outputName_dict[id_type][wp]+"_trigger_matching", "I", lenVar="n"+self.outputName_dict[id_type][wp])

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = self.inputCollection(event)
        muons = Collection(event, "tightRelIso_looseID_Muons")
        triggerObjects = self.triggerObjectCollection(event)

        selectedElectrons = OrderedDict(
            [("MVA", OrderedDict(
                [("tight", []), ("medium",[]), ("loose",[])])), 
             ("cutBased", OrderedDict(
                 [("tight", []), ("medium",[]), ("loose",[])])) ])   
        unselectedElectrons = OrderedDict([("MVA", []), ("cutBased", []) ])
        
        matched_trgObj_id_list = []
        electronCutBasedID = []
        electronMVAID = {'tight': [], 'medium': [], 'loose': []}
        nElectron = 0

        for electron in electrons:
            trigger_matching = self.triggerMatched(electron, triggerObjects)
            matched_trgObj_id_list.append(trigger_matching[1])

            # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            if electron.pt > self.electronMinPt and math.fabs(electron.eta) < self.electronMaxEta:
                # and \
                # trigger_matching[0] and matched_trgObj_id_list.count(trigger_matching[1])==1:

                setattr(electron, "trigger_matching", False)
                if trigger_matching[0] and matched_trgObj_id_list.count(trigger_matching[1])==1:
                    setattr(electron, "trigger_matching", True)
                
                if not Module.globalOptions["isData"]:
                    setattr(electron, 'genPartIdx', electron.genPartIdx) 

                # dxy = math.fabs(electron.dxy)
                # dz = math.fabs(electron.dz)

                # if math.fabs(electron.eta) < 1.479 and (dxy>0.05 or dz>0.10):
                    # unselectedElectrons.append(electron)
                    # continue
                # elif dxy>0.10 or dz>0.20:
                    # unselectedElectrons.append(electron)
                    # continue

                # gap between barrel and endcap
                # if math.fabs(electron.eta) > 1.4442 and math.fabs(electron.eta) < 1.5660:
                #     continue

                #reject electron if close-by muon
                if len(muons)>0:
                    mindr = min(map(lambda muon: deltaR(muon, electron), muons))
                    if mindr < 0.05:
                        #unselectedElectrons.append(electron)
                        continue

                #saving relIso, cutBased Id 
                nElectron+=1
                electronCutBasedID.append(electron.cutBased)

                for wp in self.WPs:
                    mva_wp = self.mva_id_wp[wp][Module.globalOptions["year"]]
                    if mva_wp == '': continue
                    electronMVAID[wp].append(getattr(electron, mva_wp))

                for id_type in self.id_type:
                    if id_type == 'MVA':
                        tight_mva = self.mva_id_wp['tight'][Module.globalOptions["year"]]
                        medium_mva = self.mva_id_wp['medium'][Module.globalOptions["year"]]
                        loose_mva = self.mva_id_wp['loose'][Module.globalOptions["year"]]
                        if getattr(electron, tight_mva) == 1: 
                            selectedElectrons[id_type]['tight'].append(electron)
                            selectedElectrons[id_type]['medium'].append(electron)
                            selectedElectrons[id_type]['loose'].append(electron)
                        elif getattr(electron, medium_mva) == 1: 
                            selectedElectrons[id_type]['medium'].append(electron)
                            selectedElectrons[id_type]['loose'].append(electron)
                        elif loose_mva != '' and getattr(electron, loose_mva) == 1: 
                            selectedElectrons[id_type]['loose'].append(electron)
                        else:
                            unselectedElectrons[id_type].append(electron)
                    elif id_type == 'cutBased':
                        if electron.cutBased == 4: 
                            selectedElectrons[id_type]['tight'].append(electron)
                            selectedElectrons[id_type]['medium'].append(electron)
                            selectedElectrons[id_type]['loose'].append(electron)
                        elif electron.cutBased == 3: 
                            selectedElectrons[id_type]['medium'].append(electron)
                            selectedElectrons[id_type]['loose'].append(electron)
                        elif electron.cutBased == 2: 
                            selectedElectrons[id_type]['loose'].append(electron)
                        else:
                            unselectedElectrons[id_type].append(electron)  
            else:
                continue
                #unselectedElectrons[id_type].append(electron)

        self.out.fillBranch("nElectron", nElectron) 
        self.out.fillBranch("electron_cutBasedID", map(lambda cutBased_id: cutBased_id, electronCutBasedID))
        for wp in self.WPs:
            self.out.fillBranch("electron_MVA_"+wp+"ID", map(lambda id: id, electronMVAID[wp]))

        for id_type in self.id_type:
            for wp in self.outputName_dict[id_type].keys():
                setattr(event, self.outputName_dict[id_type][wp], selectedElectrons[id_type][wp])
                self.out.fillBranch("n"+self.outputName_dict[id_type][wp], len(selectedElectrons[id_type][wp]))

                for variable in self.storeKinematics:
                    self.out.fillBranch(self.outputName_dict[id_type][wp]+"_"+variable, map(lambda electron: getattr(electron,variable), selectedElectrons[id_type][wp]))
                self.out.fillBranch(self.outputName_dict[id_type][wp]+"_trigger_matching",map(lambda electron: getattr(electron,"trigger_matching"),selectedElectrons[id_type][wp]))
            
                if not Module.globalOptions["isData"]:
                    for variable in self.storeTruthKeys:
                        # if len(selectedElectrons[id_type][wp])>0: print(dir(selectedElectrons[id_type][wp][0]))
                        self.out.fillBranch(self.outputName_dict[id_type][wp]+"_"+variable, map(lambda electron: getattr(electron,variable), selectedElectrons[id_type][wp]))
                
            setattr(event,"unselectedElectrons", unselectedElectrons[id_type])
        setattr(event, "nElectron", nElectron)

        return True
