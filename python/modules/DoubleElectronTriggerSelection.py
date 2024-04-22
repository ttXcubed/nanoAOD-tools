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

class DoubleElectronTriggerSelection(Module):
    def __init__(
        self,
        outputName = "DiElectronTrigger",
    ):
        self.outputName = outputName
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_ee_flag","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
       
        HLT_trigFlag = 0    
        
        if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':
            HLT_trigFlag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ

        if Module.globalOptions["year"] == '2017':
            HLT_trigFlag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL>0
            	
        if Module.globalOptions["year"] == '2018':
            HLT_trigFlag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL>0

        if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
            HLT_trigFlag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL

        self.out.fillBranch(self.outputName+"_ee_flag", HLT_trigFlag)
        setattr(event, self.outputName+"_ee_flag", HLT_trigFlag)	

        return True      
