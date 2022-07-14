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
from utils import deltaR

class NNLeptonicTopInputs():
    def __init__(self,outputName):
        self.outputName=outputName
        
        self.variableDict = {
            'lepton_pt': lambda event,lepton,jet: lepton.pt,
            'jet_pt': lambda event,lepton,jet: jet.pt,
            'jet_eta': lambda event,lepton,jet: jet.eta,
        }
    
    def bookBranches(self,out):
        out.branch(self.getLenVarName(), "I")
        for varName in self.variableDict.keys():
            out.branch(self.outputName+'_'+varName, "F", lenVar=self.getLenVarName())
            
            
    def getLenVarName(self):
        return 'n'+self.outputName
        
    def queryValues(self,event,leptonJetCombinations):
        valueDict = {
            varName: np.zeros(len(leptonJetCombinations),dtype=np.float32) for varName in self.variableDict.keys()
        }
        
        for icomb, (lepton, jet) in enumerate(leptonJetCombinations):
            for varName, varFct in self.variableDict.items():
                valueDict[varName][icomb] = varFct(event,lepton,jet)
                
        return valueDict
        
    def writeBranches(self,out,valueDict):
        out.fillBranch('n'+self.outputName, len(valueDict.values()[0]))
        for varName in self.variableDict.keys():
            out.fillBranch(self.outputName+'_'+varName,valueDict[varName])


class TopNNRecoInputs(Module):
    def __init__(self,
        leptonCollection,
        jetCollection,
        outputName,
    ):
        self.leptonCollection = leptonCollection
        self.jetCollection = jetCollection
        self.outputName=outputName
        
        self.leptonicTopInputs = NNLeptonicTopInputs(self.outputName)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.leptonicTopInputs.bookBranches(self.out)
        
        self.out.branch(self.outputName+'_gen_dRMatch', "F", lenVar=self.leptonicTopInputs.getLenVarName())
        self.out.branch(self.outputName+'_gen_leptonPtMatch', "F", lenVar=self.leptonicTopInputs.getLenVarName())
        self.out.branch(self.outputName+'_gen_jetPtMatch', "F", lenVar=self.leptonicTopInputs.getLenVarName())
        
        self.out.branch(self.outputName+'_gen_topPt', "F", lenVar=self.leptonicTopInputs.getLenVarName())

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def getGenLevelTopDecays(self,event):
        genParticles = Collection(event,"GenPart")
        
        def findTopIdx(p):
            motherIdx = p.genPartIdxMother
            while (motherIdx>=0):
                if abs(genParticles[motherIdx].pdgId)==6:
                    return motherIdx
                motherIdx = genParticles[motherIdx].genPartIdxMother
            return -1
            
        topDict = {}
        
        #NB: following will not handle hadronic taus cases
        for topIdx,genParticle in enumerate(genParticles):
            if abs(genParticle.pdgId)==6:
                topDict[topIdx] = {'top': genParticle, 'bquark':[], 'lepton': [], 'neutrino':[], 'quarks': []}
        
        for genParticle in genParticles:
            if abs(genParticle.pdgId)==5 and isLastCopy(genParticle) and fromHardProcess(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0 and topIdx in topDict.keys():
                    topDict[topIdx]['bquark'].append(genParticle)
                    
            elif abs(genParticle.pdgId)<5 and abs(genParticle.pdgId)>0 and isLastCopy(genParticle) and fromHardProcess(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0 and topIdx in topDict.keys():
                    topDict[topIdx]['quarks'].append(genParticle)
                       
            elif abs(genParticle.pdgId) in [11,13] and isPrompt(genParticle) and isLastCopy(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0 and topIdx in topDict.keys():
                    topDict[topIdx]['lepton'].append(genParticle)
                    
            elif abs(genParticle.pdgId) in [12,14] and isPrompt(genParticle) and isLastCopy(genParticle):
                topIdx = findTopIdx(genParticle)
                if topIdx>=0 and topIdx in topDict.keys():
                    topDict[topIdx]['neutrino'].append(genParticle)
            
        hadronicTops = []
        leptonicTops = []
        otherTops = []
            
        for idx in topDict.keys():
            top = topDict[idx]
            if len(top['quarks'])==2 and len(top['lepton'])==0 and len(top['neutrino'])==0 and len(top['bquark'])==1:
                hadronicTops.append(top)
            elif len(top['quarks'])==0 and len(top['lepton'])==1 and len(top['neutrino'])==1 and len(top['bquark'])==1:
                leptonicTops.append(top)
            else:
                # rare more complex decays (e.g. involving taus or FSR)
                otherTops.append(top)
                
        return {"hadronic": hadronicTops, "leptonic": leptonicTops, "other": otherTops}

    def makeLeptonJetCombinations(self, leptons, jets):
        combinations = []
        for lepton in leptons:
            for jet in jets:
                #NB: can do some filtering here if needed; ie. min/max mT
                combinations.append([lepton,jet])
        return combinations

    def analyze(self, event):
    
        leptons = self.leptonCollection(event)
        jets = self.jetCollection(event)
        
        leptonJetCombinations = self.makeLeptonJetCombinations(leptons,jets)
        
        valueDict = self.leptonicTopInputs.queryValues(event,leptonJetCombinations)
        
        genLeptonicTopDecays = self.getGenLevelTopDecays(event)["leptonic"]
        
        gen_dRMatch = 10.*np.ones(len(leptonJetCombinations),dtype=np.float32)
        gen_leptonPtMatch = np.zeros(len(leptonJetCombinations),dtype=np.float32)
        gen_jetPtMatch = np.zeros(len(leptonJetCombinations),dtype=np.float32)
        gen_topPt = np.zeros(len(leptonJetCombinations),dtype=np.float32)
        
        for icomb, (lepton, jet) in enumerate(leptonJetCombinations):
            for genLeptonicTopDecay in genLeptonicTopDecays:
                genBquark = genLeptonicTopDecay['bquark'][0]
                genLepton = genLeptonicTopDecay['lepton'][0]
                #genNeutrino = genLeptonicTopDecay['neutrino'][0]
                genTopQuark =  genLeptonicTopDecay['top']
                
                dRLepton = deltaR(lepton,genLepton)
                dRJet = deltaR(jet,genBquark)

                #bias matching towards lepton since jets are less precise
                dRMatch = math.sqrt((dRLepton/0.1)**2+(dRJet/0.4)**2)
                
                if dRMatch<gen_dRMatch[icomb]:
                    gen_dRMatch[icomb] = dRMatch
                    gen_leptonPtMatch[icomb] = lepton.pt/genLepton.pt
                    gen_jetPtMatch[icomb] = jet.pt/genBquark.pt
                    gen_topPt[icomb] = genTopQuark.pt
                    
        self.leptonicTopInputs.writeBranches(self.out,valueDict)
        
        self.out.fillBranch(self.outputName+'_gen_dRMatch', gen_dRMatch)
        self.out.fillBranch(self.outputName+'_gen_leptonPtMatch', gen_leptonPtMatch)
        self.out.fillBranch(self.outputName+'_gen_jetPtMatch', gen_jetPtMatch)
        self.out.fillBranch(self.outputName+'_gen_topPt', gen_topPt)
        
        return True
        
        
