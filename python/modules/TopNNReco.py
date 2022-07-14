import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from gen_helper import isHardProcess, isLastCopy, isPrompt, fromHardProcess

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

class GenTops():
    def __init__(self,outputName):
        self.outputName=outputName
        
        #check gen matched top; resolved/boosted
        #neutrino pt,kt,phi in top restframe wrt jet
        #gen top pt
        
    def bookBranches(self,out,skipLenBranch=True):
        if not skipLenBranch:
            out.branch('n'+self.outputName, "I")
        for varName, varFct in self.variableList:
            out.branch(self.outputName+'_'+varName, "F", lenVar='n'+self.outputName)
        
    def queryValues(self,event,leptons,jets):
        genParticles = Collection(event,"GenPart")
        
        def findTopIdx(p):
            motherIdx = p.genPartIdxMother
            while (motherIdx>=0):
                if abs(genParticles[motherIdx].pdgId)==6:
                    return motherIdx
                motherIdx = genParticles[motherIdx].genPartIdxMother
            return -1
            
        topDict = {}
        
        #NB: following will ignore hadronic taus
        for genParticle in genParticles:
            if abs(genParticle.pdgId)==5 and isLastCopy(genParticle) and fromHardProcess(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0:
                    if topIdx not in topDict.keys():
                        topDict[topIdx] = {'bquark':[], 'lepton': [], 'neutrino':[], 'quarks': []}
                    topDict[topIdx]['bquark'].append(genParticle)
            elif abs(genParticle.pdgId)<5 and abs(genParticle.pdgId)>0 and isLastCopy(genParticle) and fromHardProcess(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0:
                    if topIdx not in topDict.keys():
                        topDict[topIdx] = {'bquark':[], 'lepton': [], 'neutrino':[], 'quarks': []}
                    topDict[topIdx]['quarks'].append(genParticle)
                       
            elif abs(genParticle.pdgId) in [11,13] and isPrompt(genParticle) and isLastCopy(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0:
                    if topIdx not in topDict.keys():
                        topDict[topIdx] = {'bquark':[], 'lepton': [], 'neutrino':[], 'quarks': []}
                    topDict[topIdx]['lepton'].append(genParticle)
                    
            elif abs(genParticle.pdgId) in [12,14] and isPrompt(genParticle) and isLastCopy(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0:
                    if topIdx not in topDict.keys():
                        topDict[topIdx] = {'bquark':[], 'lepton': [], 'neutrino':[], 'quarks': []}
                    topDict[topIdx]['neutrino'].append(genParticle)
            
        for idx in topDict.keys():
            top = topDict[idx]
            if len(top['quarks'])==2 and len(top['lepton'])==0 and len(top['neutrino'])==0 and len(top['bquark'])==1:
                top['leptonic']=False
            elif len(top['quarks'])==0 and len(top['lepton'])==1 and len(top['neutrino'])==1 and len(top['bquark'])==1:
                top['leptonic']=True
            else:
                del topDict[idx]
        print len(topDict)
        '''
        for ilepton, lepton in enumerate(leptons):
            for ijet, jet in enumerate(jets):
                for ivar, (varName, varFct) in enumerate(self.variableList):
                    values[ilepton*len(jets)+ijet,ivar] = varFct(event,lepton,jet)
        return values
        '''
        
    def writeBranches(self,out,values):
        out.fillBranch('n'+self.outputName, values.shape[0])
        for ivar, (varName, varFct) in enumerate(self.variableList):
            out.fillBranch(self.outputName+'_'+varName,values[:,ivar])



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
        self.genTops = GenTops('genTop')
        

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
        
        self.genTops.queryValues(event,leptons,jets)
        
        self.leptonicTopInputs.writeBranches(self.out,values)
        return True
        
        
