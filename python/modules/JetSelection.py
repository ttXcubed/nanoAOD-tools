import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class JetSelection(Module):

    NONE = -1
    LOOSE = 0
    TIGHT = 1
    TIGHTLEPVETO = 2
    
    def __init__(
         self,
         inputCollection=lambda event: Collection(event, "Jet"),
         leptonCollectionDRCleaning=lambda event: [],
         outputName_list=["selectedJets","unselectedJets"],
         jetMinPt=30.,
         jetMaxEta=4.8,
         dRCleaning=0.4,
         storeKinematics=['pt', 'eta','phi','mass'],
         jetId=LOOSE,
         fatFlag=True,
         metInput = lambda event: Object(event, "MET"),
         storeTruthKeys=[]
     ):

        self.inputCollection = inputCollection
        self.leptonCollectionDRCleaning = leptonCollectionDRCleaning
        self.outputName_list = outputName_list
        self.jetMinPt = jetMinPt
        self.jetMaxEta = jetMaxEta
        self.dRCleaning = dRCleaning
        self.storeKinematics = storeKinematics
        self.jetId=jetId
        self.fatFlag = fatFlag
        self.metInput = metInput
        self.storeTruthKeys = storeTruthKeys

        #loose jet ID does not exists for UL -> accepting all jets  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
        if self.jetId==JetSelection.LOOSE and Module.globalOptions['year'] in ['2016','2016preVFP', '2017', '2018', '2022']:
            self.jetId = JetSelection.NONE

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("MET_energy", "F")
        for outputName in self.outputName_list:
            self.out.branch("n"+outputName, "I")

            for variable in self.storeKinematics:
                if variable=='_index': 
                    self.out.branch(outputName+variable, "I", lenVar="n"+outputName)
                    continue
                self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)
            if not Module.globalOptions["isData"]:
                for variable in self.storeTruthKeys:
                    self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)
                    

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)
        met = self.metInput(event)

        selectedJets = []
        unselectedJets = []

        leptonsForDRCleaning = self.leptonCollectionDRCleaning(event)
        for i, jet in enumerate(jets):
            # -- lepton cleaning
            jet_radius = self.dRCleaning
            # in case of HOTVR, the radius is be calculated as 600/jet pT
            if self.dRCleaning == None:
                jet_radius = 600./ jet.pt if 600./ jet.pt <= 1.5 else 1.5 

            minDeltaRSubtraction = 999.
            if len(leptonsForDRCleaning) > 0:
                # mindphi = min(map(lambda lepton: math.fabs(deltaPhi(lepton, jet)), leptonsForDRCleaning))
                mindphi = min(map(lambda lepton: deltaPhi(lepton, jet), leptonsForDRCleaning))
                mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsForDRCleaning))

                setattr(jet, "minDPhiClean", mindphi)
                setattr(jet, "minDRClean", mindr)

                if mindr < jet_radius:
                    unselectedJets.append(jet)
                    continue

            else:
                setattr(jet, "minDPhiClean", 100)
                setattr(jet, "minDRClean", 100)
            # --

            if jet.pt < self.jetMinPt:
                unselectedJets.append(jet)
                continue
        
            if math.fabs(jet.eta) > self.jetMaxEta:
                unselectedJets.append(jet)
                continue

            if self.jetId>=0 and ((jet.jetId & (1 << self.jetId)) == 0):
                unselectedJets.append(jet)
                continue

            selectedJets.append(jet)

        def metP4(obj):
            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM(obj.pt,0,obj.phi,0)
            return p4

        #print((metP4(met)).E())
        self.out.fillBranch("MET_energy", (metP4(met)).E())

        for outputName, jet_list in zip(self.outputName_list, [selectedJets, unselectedJets]):
            setattr(event, outputName, jet_list)
            self.out.fillBranch("n"+outputName, len(jet_list))
            for variable in self.storeKinematics:
                if Module.globalOptions["isData"] and variable=='genJetAK8Idx': continue
                if variable == '_index': 
                    self.out.fillBranch(outputName+variable, map(lambda jet: getattr(jet, variable), jet_list))
                    continue
                self.out.fillBranch(outputName+"_"+variable, map(lambda jet: getattr(jet, variable), jet_list))

            if not Module.globalOptions["isData"]:
                for variable in self.storeTruthKeys:
                    self.out.fillBranch(outputName+"_"+variable, map(lambda jet: getattr(jet, variable), jet_list))

        return True

