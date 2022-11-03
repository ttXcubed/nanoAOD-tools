import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getGraph, getHist, combineHist2D, getSFXY, deltaR

class MuonSelection(Module):
    VERYTIGHT = 1
    TIGHT = 1
    MEDIUM = 2
    LOOSE = 3
    NONE = 4
    INV = 5

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Muon"),
        outputName_list=["tightMuons","mediumMuons","looseMuons"],
        triggerMatch=False,
        #muonID=TIGHT,
        #muonIso=TIGHT,
        muonMinPt=25.,
        muonMaxEta=2.4,
        storeKinematics= [],#['pt','eta'],
    ):
        
        self.inputCollection = inputCollection
        self.outputName_list = outputName_list
        self.muonMinPt = muonMinPt
        self.muonMaxEta = muonMaxEta
        self.storeKinematics = storeKinematics
        self.triggerMatch = triggerMatch
        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []
        
	
    def triggerMatched(self, muon, trigger_object): 
        #return 2 arguments: 
        #   -1st: flag-->True, if muon matches a trigger object (deltaR<0.3); False, otherwise
        #   -2nd: idx of the matched trigger object
        min_deltaR, matchedTrgObj_id = 1., None
        if self.triggerMatch:
            #trig_deltaR = math.pi
            for itrgObj, trig_obj in enumerate(trigger_object):
                if abs(trig_obj.id) != 13:
                    continue
                #trig_deltaR = min(trig_deltaR, deltaR(trig_obj, muon))
                else:
                    deltaR_trgObj_lep = deltaR(muon,trig_obj)
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
        
        for var in ['PFRelIso04', 'tightID', 'mediumID', 'looseID']:        
            self.out.branch("muon_"+var,"F",lenVar="nMuon")

        for outputName in self.outputName_list:
            self.out.branch("n"+outputName, "I")

            for variable in self.storeKinematics:
                self.out.branch(outputName+"_"+variable,"F",lenVar="n"+outputName)
                if not Module.globalOptions["isData"]:
                        self.out.branch(outputName+"_genPartFlav","F",lenVar="n"+outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = self.inputCollection(event)
        triggerObjects = self.triggerObjectCollection(event)

        selectedMuons = {'tight': [], 'medium': [], 'loose': []}
        unselectedMuons = []
        matched_trgObj_id_list = []
        muonRelIso = []
        muonID = {'tight': [], 'medium': [], 'loose': []}
        nMuon = 0
        
        #https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
        for muon in muons:
            matched_trgObj_id_list.append(self.triggerMatched(muon, triggerObjects)[1])
            
            #baseline selection (pT, eta, trigger matching requirement)
            if muon.pt>self.muonMinPt and math.fabs(muon.eta)<self.muonMaxEta and \
                self.triggerMatched(muon, triggerObjects)[0] and matched_trgObj_id_list.count(self.triggerMatched(muon, triggerObjects)[1])==1:

                #saving relIso, cutBased Id 
                nMuon+=1
                muonRelIso.append(muon.pfRelIso04_all)
                muonID['tight'].append(muon.tightId)
                muonID['medium'].append(muon.mediumId)
                muonID['loose'].append(muon.looseId)

                #selectedMuons: tightIso, differentID --> https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2020/085 (4 top in dilepton final state)
                if muon.pfRelIso04_all<0.15:
                    if muon.tightId==1: 
                        selectedMuons['tight'].append(muon)
                        selectedMuons['medium'].append(muon)
                        selectedMuons['loose'].append(muon)
                    elif muon.mediumId==1:
                        selectedMuons['medium'].append(muon)
                        selectedMuons['loose'].append(muon)
                    elif muon.looseId==1:
                        selectedMuons['loose'].append(muon)
                    else:
                        unselectedMuons.append(muon)
                # elif muon.mediumId==1 and muon.pfRelIso04_all<0.20:
                #     selectedMuons['medium'].append(muon)
                #     selectedMuons['loose'].append(muon)
                # elif muon.looseId==1 and muon.pfRelIso04_all<0.25:
                #     print('muone buon')
                #     selectedMuons['loose'].append(muon)

            else:
                unselectedMuons.append(muon)

        self.out.fillBranch("nMuon",nMuon) 
        self.out.fillBranch("muon_PFRelIso04", map(lambda iso: iso, muonRelIso))
        for wp in muonID.keys():
            self.out.fillBranch("muon_"+wp+"ID", map(lambda id: id, muonID[wp]))

        for outputName, muon_ID in zip(self.outputName_list, ['tight','medium','loose']):
            self.out.fillBranch("n"+outputName,len(selectedMuons[muon_ID]))
		
            for variable in self.storeKinematics:
                self.out.fillBranch(outputName+"_"+variable,map(lambda muon: getattr(muon,variable),selectedMuons[muon_ID]))
                if not Module.globalOptions["isData"]:
                    self.out.fillBranch(outputName+"_genPartFlav",map(lambda muon: getattr(muon,'genPartFlav'),selectedMuons[muon_ID]))
            setattr(event,outputName,selectedMuons[muon_ID])
            #setattr(event, muon_ID+'Muons',selectedMuons[muon_ID])

        setattr(event,"unselectedMuons",unselectedMuons)
        setattr(event,'nMuon',nMuon)

        return True

