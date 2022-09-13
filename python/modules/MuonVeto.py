import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class MuonVeto(Module):

    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Muon"),
        outputName = "vetoMuons",
        muonMinPt = 10.,
        muonMaxEta = 2.4,
        storeKinematics=['pt','eta','charge','phi','mass'],
    ):
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.muonMinPt = muonMinPt
        self.muonMaxEta = muonMaxEta
        self.storeKinematics = storeKinematics
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        #self.out.branch("n"+self.outputName,"I")
        
	self.out.branch("n"+self.outputName, "I")
	self.out.branch("nunselectedMuons", "I")

	for variable in self.storeKinematics:
	    self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
	    self.out.branch("unselectedMuons"+"_"+variable,"F",lenVar="nunselectedMuons")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = self.inputCollection(event)
        
        vetoMuons = []
        vetounselectedMuons = []
        
        #https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
        for muon in muons:
            if muon.pt>self.muonMinPt and math.fabs(muon.eta)<self.muonMaxEta and muon.isPFcand==1 and (muon.pfRelIso04_all<0.25):
                vetoMuons.append(muon)
            else:
                vetounselectedMuons.append(muon)
  
	self.out.fillBranch("n"+self.outputName,len(vetoMuons))
	self.out.fillBranch("nunselectedMuons",len(vetounselectedMuons))
		
	for variable in self.storeKinematics:
	    self.out.fillBranch(self.outputName+"_"+variable,map(lambda muon: getattr(muon,variable),vetoMuons))
	    self.out.fillBranch("unselectedMuons""_"+variable,map(lambda muon: getattr(muon,variable),vetounselectedMuons))
        #self.out.fillBranch("n"+self.self.outputName,len(selectedMuons))
        
        setattr(event,self.outputName,vetoMuons)
        setattr(event,self.outputName+"_unselected",vetounselectedMuons)

        return True
        
