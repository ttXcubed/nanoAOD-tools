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
from utils import deltaR,deltaPhi

class NNLeptonicTopInputs():
    def __init__(self,outputName):
        self.outputName=outputName
        
        self.variableDict = {
            #event vars
            'rho': lambda event,lepton,jet: event.fixedGridRhoFastjetAll,
            'npv': lambda event,lepton,jet: event.PV_npvsGood,
            
            
            #jet vars
            'jet_area': lambda event,lepton,jet: jet.area,
            'jet_btagDeepFlavB': lambda event,lepton,jet: jet.btagDeepFlavB,
            'jet_qgl': lambda event,lepton,jet: jet.qgl,
            'jet_nConstituents': lambda event,lepton,jet: jet.nConstituents,
            
            'jet_chEmEF': lambda event,lepton,jet: jet.chEmEF,
            'jet_chFPV0EF': lambda event,lepton,jet: jet.chFPV0EF,
            
            'jet_chHEF': lambda event,lepton,jet: jet.chHEF,
            
            'jet_neEmEF': lambda event,lepton,jet: jet.neEmEF,
            'jet_neHEF': lambda event,lepton,jet: jet.neHEF,
            
            'jet_pt': lambda event,lepton,jet: jet.pt,
            'jet_abseta': lambda event,lepton,jet: math.fabs(jet.eta),
            'jet_mass': lambda event,lepton,jet: jet.mass,
            
            
            #lepton vars
            'lepton_pt': lambda event,lepton,jet: lepton.pt,
            'lepton_abseta': lambda event,lepton,jet: math.fabs(lepton.pt),
            
            'lepton_dxy': lambda event,lepton,jet:lepton.dxy,
            'lepton_dxysig': lambda event,lepton,jet:lepton.dxy/(1e-4+math.fabs(lepton.dxyErr)),
            'lepton_dz': lambda event,lepton,jet:lepton.dxy,
            'lepton_dzsig': lambda event,lepton,jet:lepton.dz/(1e-4+math.fabs(lepton.dzErr)),
            
            'lepton_isMuon': lambda event,lepton,jet:1. if lepton.isMuon else 0.,
            
            
            #electron specific vars
            'ele_deltaEtaSC': lambda event,lepton,jet:math.log(max(1e-5,math.fabs(lepton.deltaEtaSC))) if lepton.isElectron else 0.,
            'ele_dr03EcalRecHitSumEt': lambda event,lepton,jet:lepton.dr03EcalRecHitSumEt if lepton.isElectron else 0.,
            'ele_dr03HcalDepth1TowerSumEt': lambda event,lepton,jet:lepton.dr03HcalDepth1TowerSumEt if lepton.isElectron else 0.,
            'ele_dr03TkSumPt': lambda event,lepton,jet:lepton.dr03TkSumPt if lepton.isElectron else 0.,
            'ele_dr03TkSumPtHEEP': lambda event,lepton,jet:lepton.dr03TkSumPtHEEP if lepton.isElectron else 0.,
            
            'ele_eInvMinusPInv': lambda event,lepton,jet:lepton.eInvMinusPInv if lepton.isElectron else 0.,
            'ele_hoe': lambda event,lepton,jet:lepton.hoe if lepton.isElectron else 0.,
            'ele_ip3d': lambda event,lepton,jet:lepton.ip3d if lepton.isElectron else 0.,
            'ele_jetPtRelv2': lambda event,lepton,jet: math.log(max(1e-5,lepton.jetPtRelv2)) if lepton.isElectron else 0.,
            'ele_jetRelIso': lambda event,lepton,jet:math.log(max(1e-5,lepton.jetRelIso)) if lepton.isElectron else 0.,
            'ele_lostHits': lambda event,lepton,jet:lepton.lostHits if lepton.isElectron else 0.,
            'ele_r9': lambda event,lepton,jet:lepton.r9 if lepton.isElectron else 0.,
            'ele_scEtOverPt': lambda event,lepton,jet:lepton.scEtOverPt if lepton.isElectron else 0.,
            'ele_seedGain': lambda event,lepton,jet:lepton.seedGain if lepton.isElectron else 0.,
            'ele_sieie': lambda event,lepton,jet:lepton.sieie if lepton.isElectron else 0.,
            'ele_sip3d': lambda event,lepton,jet:lepton.sip3d if lepton.isElectron else 0.,
            
            'ele_pfRelIso03_all': lambda event,lepton,jet: math.log(max(1e-5,lepton.pfRelIso03_all)) if lepton.isElectron else 0.,
            'ele_pfRelIso03_chg': lambda event,lepton,jet: math.log(max(1e-5,lepton.pfRelIso03_chg)) if lepton.isElectron else 0.,
        
        
            #muon specific vars
            'mu_ip3d': lambda event,lepton,jet:lepton.ip3d if lepton.isMuon else 0.,
            'mu_jetPtRelv2': lambda event,lepton,jet: math.log(max(1e-5,lepton.jetPtRelv2)) if lepton.isMuon else 0.,
            'mu_jetRelIso': lambda event,lepton,jet:math.log(max(1e-5,lepton.jetRelIso)) if lepton.isMuon else 0.,
            'mu_miniPFRelIso_all': lambda event,lepton,jet: math.log(max(1e-5,lepton.miniPFRelIso_all)) if lepton.isMuon else 0.,
            'mu_miniPFRelIso_chg': lambda event,lepton,jet: math.log(max(1e-5,lepton.miniPFRelIso_chg)) if lepton.isMuon else 0.,
            'mu_nStations': lambda event,lepton,jet:lepton.nStations if lepton.isMuon else 0.,
            'mu_nTrackerLayers': lambda event,lepton,jet:lepton.nTrackerLayers if lepton.isMuon else 0.,
            'mu_pfRelIso03_all': lambda event,lepton,jet: math.log(max(1e-5,lepton.pfRelIso03_all)) if lepton.isMuon else 0.,
            'mu_pfRelIso03_chg': lambda event,lepton,jet: math.log(max(1e-5,lepton.pfRelIso03_chg)) if lepton.isMuon else 0.,
            'mu_pfRelIso04_all': lambda event,lepton,jet: math.log(max(1e-5,lepton.pfRelIso04_all)) if lepton.isMuon else 0.,
            'mu_segmentComp': lambda event,lepton,jet:lepton.segmentComp if lepton.isMuon else 0.,
            'mu_sip3d': lambda event,lepton,jet:lepton.sip3d if lepton.isMuon else 0.,
            'mu_tkRelIso': lambda event,lepton,jet: math.log(max(1e-5,lepton.tkRelIso)) if lepton.isMuon else 0.,
            
            
            #correlations
            'leptonJet_pt': lambda event,lepton,jet: (lepton.p4()+jet.p4()).Pt(),
            'leptonJet_ptratio': lambda event,lepton,jet: jet.pt/lepton.pt,
            'leptonJet_mass': lambda event,lepton,jet: (lepton.p4()+jet.p4()).M(),
            'leptonJet_mt': lambda event,lepton,jet: (lepton.p4()+jet.p4()).Mt(),
            'leptonJet_abseta': lambda event,lepton,jet: math.fabs((lepton.p4()+jet.p4()).Eta()),
            
            'leptonJet_dR': lambda event,lepton,jet: deltaR(lepton,jet),
            'leptonJet_dPhi': lambda event,lepton,jet: deltaPhi(lepton,jet),
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
        muonCollection,
        electronCollection,
        jetCollection,
        outputName,
    ):
        self.muonCollection = muonCollection
        self.electronCollection = electronCollection
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
        muons = self.muonCollection(event)
        for muon in muons:
            setattr(muon,'isMuon',True)
            setattr(muon,'isElectron',False)
        electrons = self.electronCollection(event)
        for electron in electrons:
            setattr(electron,'isMuon',False)
            setattr(electron,'isElectron',True)
            
        leptons = muons+electrons
        
       
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
        
        
