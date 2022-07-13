import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class NNLeptonicTopInputs():
    def __init__(self,outputName):
        self.outputName=outputName
        
        #note: order matters and should not be changed in other functions!
        self.variableList = [
            ['lepton_pt', lambda event,lepton,jet: lepton.pt],
            ['jet_pt', lambda event,lepton,jet: jet.pt],
            ['jet_eta', lambda event,lepton,jet: jet.eta]
        ]
    
    def bookBranches(self,out):
        out.branch('n'+self.outputName, "I")
        for varName, varFct in self.variableList:
            out.branch(self.outputName+'_'+varName, "F", lenVar='n'+self.outputName)
        
    def queryValues(self,event,leptons,jets):
        # [ncombinations, nvariables]
        values = np.zeros(
            (len(leptons)*len(jets),len(self.variableList)),
            dtype=np.float32
        )  

        for ilepton, lepton in enumerate(leptons):
            for ijet, jet in enumerate(jets):
                for ivar, (varName, varFct) in enumerate(self.variableList):
                    values[ilepton*len(jets)+ijet,ivar] = varFct(event,lepton,jet)
        return values
        
    def writeBranches(self,out,values):
        out.fillBranch('n'+self.outputName, values.shape[0])
        for ivar, (varName, varFct) in enumerate(self.variableList):
            out.fillBranch(self.outputName+'_'+varName,values[:,ivar])
'''
class NNGenLeptonicTop():
    def __init__(self):
        pass
        
        #check gen matched top; resolved/boosted
        #neutrino pt,kt,phi in top restframe wrt jet
        #gen top pt
        
    def bookBranches(self,out):
        for varName in self.variableDict.keys():
            self.out.branch(self.outputName+'_'+varName, "F")
        
    def queryValues(event,leptons,jets):
        valueDict = {}
        for varName, varFct in self.variableDict.items():
            valueDict[varName] = varFct(event,lepton,jet)
        return valueDict
        
    def writeBranches(self,out,valueDict):
        for varName in self.variableDict.keys():
            self.out.fillBranch(self.outputName+'_'+varName,valueDict[varName])
'''


class TopNNRecoInputs(Module):
    def __init__(self,
        leptonCollection,
        jetCollection,
        outputName,
    ):
        self.leptonCollection = leptonCollection
        self.jetCollection = jetCollection
        self.outputName=outputName
        
        self.leptonicTopInputs = NNLeptonicTopInputs('leptonicTopNN')

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.leptonicTopInputs.bookBranches(self.out)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
    
        leptons = self.leptonCollection(event)
        jets = self.jetCollection(event)
        
        values = self.leptonicTopInputs.queryValues(event,leptons,jets)
        #print values
        
        self.leptonicTopInputs.writeBranches(self.out,values)
        return True
        
        
