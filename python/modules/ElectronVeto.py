import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class ElectronVeto(Module):

    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Electron"),
        outputName = "vetoElectrons",
        electronMinPt = 10.,
        electronMaxEta = 2.4,
        storeKinematics=['pt','eta','charge','phi','mass'],
    ):
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.electronMinPt = electronMinPt
        self.electronMaxEta = electronMaxEta
        self.storeKinematics = storeKinematics
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        #self.out.branch("n"+self.outputName,"I")
        
        self.out.branch("n"+self.outputName, "I")
	self.out.branch("nunselectedElectrons", "I")

	for variable in self.storeKinematics:
	    self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
	    self.out.branch("unselectedElectrons"+"_"+variable,"F",lenVar="nunselectedElectrons")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = self.inputCollection(event)
        
        vetoElectrons = []
        vetounselectedElectrons = []
        
        #use electron cut-based veto ID
        for electron in electrons:
            if electron.pt>self.electronMinPt and math.fabs(electron.eta)<self.electronMaxEta and electron.cutBased>=1:
                vetoElectrons.append(electron)
            else:
                vetounselectedElectrons.append(electron)
  
        #self.out.fillBranch("n"+self.outputName,len(selectedElectrons))
        
        self.out.fillBranch("n"+self.outputName,len(vetoElectrons))
	self.out.fillBranch("nunselectedElectrons",len(vetounselectedElectrons))
		
	for variable in self.storeKinematics:
	    self.out.fillBranch(self.outputName+"_"+variable,map(lambda muon: getattr(muon,variable),vetoElectrons))
	    self.out.fillBranch("unselectedElectrons""_"+variable,map(lambda muon: getattr(muon,variable),vetounselectedElectrons))
        
        setattr(event,self.outputName,vetoElectrons)
        setattr(event,self.outputName+"_unselected",vetounselectedElectrons)

        return True
        
