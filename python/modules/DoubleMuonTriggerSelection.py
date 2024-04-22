import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from functools import reduce
from utils import getHist,combineHist2D,getSFXY,deltaR,deltaPt

class DoubleMuonTriggerSelection(Module):
    def __init__(
        self,
        outputName = "Trigger", 
    ):
        self.outputName = outputName
         
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_mumu_flag","I")
   
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        HLT_trigFlag = 0    
        #to be modify: different ERAS have different trigger flag --> https://gitlab.cern.ch/tdr/notes/AN-19-140/-/raw/miniAODv2/TrigEffAN_temp.pdf
        if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':
            if Module.globalOptions["era"]=='H':
                HLT_trigFlag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ or event.HLT_IsoMu24 or event.HLT_IsoTkMu24
            else: 
                HLT_trigFlag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL or event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL or event.HLT_IsoMu24 or event.HLT_IsoTkMu24

        if Module.globalOptions["year"] == '2017':
            if Module.globalOptions["era"]=='B':
                HLT_trigFlag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or event.HLT_IsoMu27
            else:
                HLT_trigFlag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or event.HLT_IsoMu27

        if Module.globalOptions["year"] == '2018':
            HLT_trigFlag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or event.HLT_IsoMu24

        if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
            HLT_trigFlag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8

        self.out.fillBranch(self.outputName+"_mumu_flag", HLT_trigFlag)
        setattr(event, self.outputName+"_mumu_flag", HLT_trigFlag)	

        return True
    