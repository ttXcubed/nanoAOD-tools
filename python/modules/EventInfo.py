import os
import sys
import math
import json
import ROOT
import random
import yaml

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class EventInfo(Module):
    def __init__(
        self,
        storeVariables = []
    ):
        self.storeVariables = storeVariables
        
        
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.nGenWeights = 0	   
        self.out = wrappedOutputTree
        for variable in self.storeVariables:
            variable[0](self.out)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    	if not Module.globalOptions["isData"]:
            nGenWeight_parameter = ROOT.TParameter(float)("sumGenWeights", self.nGenWeights)
            outputFile.cd()
            nGenWeight_parameter.Write() 
        else:
            pass
        
    def analyze(self, event):
        #nGenWeights = 0 #sum of the genWeights for all the samples
        
        if not Module.globalOptions["isData"]:
            self.nGenWeights += event.Generator_weight
        
        else: 
            with open(os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/13TeV_UL_Era_runNumber.yaml") as yaml_f:
                try:
                    era_yaml = yaml.safe_load(yaml_f)
                except yaml.YAMLError as exc:
                    print(exc)
        
            for era in era_yaml[Module.globalOptions['year']]:
                if event.run>=era_yaml[Module.globalOptions['year']][era][0] and event.run<=era_yaml[Module.globalOptions['year']][era][1]:
                    Module.globalOptions['era']=era


        for variable in self.storeVariables:
            variable[1](self.out,event)
        return True
        
