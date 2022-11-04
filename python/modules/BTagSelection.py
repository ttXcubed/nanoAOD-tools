import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi
from collections import OrderedDict

class BTagSelection(Module):
    #tight DeepFlav WP (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy)
    LOOSE=0
    MEDIUM=1
    TIGHT=2
    
    def __init__(
        self,
        btaggingWP = {},
        inputCollection=lambda event: Collection(event, "Jet"),
        outputName_list= [], #"btaggedJets",
        jetMinPt=30.,
        jetMaxEta=2.4,
        workingpoint = [], 
        storeKinematics=['pt', 'eta','phi','mass'],
        storeTruthKeys=[]
    ):
        self.btaggingWP = btaggingWP
        self.inputCollection = inputCollection
        #self.outputName = outputName
        self.outputName_list = outputName_list
        self.jetMinPt = jetMinPt
        self.jetMaxEta = jetMaxEta
        self.storeKinematics = storeKinematics
        self.storeTruthKeys = storeTruthKeys
        self.workingpoint = workingpoint
        
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        for outputName in self.outputName_list:
            self.out.branch("n"+outputName, "I")
            if not Module.globalOptions["isData"]:
                self.out.branch(outputName+"_genJetIdx", "I", lenVar="n"+outputName)
            self.out.branch(outputName+"_btagDeepFlavB", "F", lenVar="n"+outputName)
            for variable in self.storeKinematics:
                self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)
            if not Module.globalOptions["isData"]:
                for variable in self.storeTruthKeys:
                    self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        bJets = OrderedDict([('tight',[]), ('medium',[]), ('loose',[])])
        lJets = []


        for jet in jets:
            if jet.pt<self.jetMinPt:
                lJets.append(jet)
                continue
        
            if math.fabs(jet.eta) > self.jetMaxEta:
                lJets.append(jet)
                continue

            if jet.btagDeepFlavB>self.btaggingWP[2]:
                bJets['tight'].append(jet)
                setattr(jet,"is_tightBTagged",True)
            
            if jet.btagDeepFlavB>self.btaggingWP[1]:
                bJets['medium'].append(jet)
                setattr(jet,"is_mediumBTagged",True)
                
            if jet.btagDeepFlavB>self.btaggingWP[0]:
                bJets['loose'].append(jet)
                setattr(jet,"is_looseBTagged",True)
            
            #if not self.taggerFct(jet):
            #    lJets.append(jet)
            #    continue
            
            #bJets.append(jet)
            
        # for wp in bJets.keys():
        #     for jet in bJets[wp]:
                
        # for jet in lJets:
        #     setattr(jet,"is_not_b_tagged",False)

        for outputName, bJet_type in zip(self.outputName_list, ['tight','medium','loose']):

            self.out.fillBranch("n"+outputName, len(bJets[bJet_type]))

            if not Module.globalOptions["isData"]:
                self.out.fillBranch(outputName+"_genJetIdx", map(lambda jet: getattr(jet, 'genJetIdx'), bJets[bJet_type]))
            self.out.fillBranch(outputName+"_btagDeepFlavB", map(lambda jet: getattr(jet, 'btagDeepFlavB'), bJets[bJet_type]))
            #print("b jet idxtogen ", list(map(lambda jet: getattr(jet, 'genJetIdx'), bJets)))        	

            for variable in self.storeKinematics:
                self.out.fillBranch(outputName+"_"+variable, map(lambda jet: getattr(jet, variable), bJets[bJet_type]))

            if not Module.globalOptions["isData"]:
                for variable in self.storeTruthKeys:
                    self.out.fillBranch(outputName+"_"+variable, map(lambda jet: getattr(jet, variable), bJets[bJet_type]))


                setattr(event, outputName, bJets[bJet_type])
                #setattr(event, self.outputName+"_unselected", lJets)
                    
        return True

