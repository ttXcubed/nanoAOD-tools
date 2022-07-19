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
            'lepton_abseta': lambda event,lepton,jet: math.fabs(lepton.eta),
            
            'lepton_dxy': lambda event,lepton,jet:lepton.dxy,
            'lepton_dxysig': lambda event,lepton,jet:lepton.dxy/(1e-4+math.fabs(lepton.dxyErr)),
            'lepton_dz': lambda event,lepton,jet:lepton.dz,
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


class NNHadronicTopInputs():
    def __init__(self,outputName):
        self.outputName=outputName
        
        def minAngleInBJJ(bjet, jet1, jet2):
            bjjP4 = bjet.p4()+jet1.p4()+jet2.p4()
            bjetP4 = bjet.p4()
            bjetP4.Boost(-bjjP4.BoostVector())
            jet1P4 = jet1.p4()
            jet1P4.Boost(-bjjP4.BoostVector())
            jet2P4 = jet2.p4()
            jet2P4.Boost(-bjjP4.BoostVector())
            
            return min([bjetP4.Angle(jet1P4.Vect()),bjetP4.Angle(jet2P4.Vect()),jet1P4.Angle(jet1P4.Vect())])
            
        def maxAngleInBJJ(bjet, jet1, jet2):
            bjjP4 = bjet.p4()+jet1.p4()+jet2.p4()
            bjetP4 = bjet.p4()
            bjetP4.Boost(-bjjP4.BoostVector())
            jet1P4 = jet1.p4()
            jet1P4.Boost(-bjjP4.BoostVector())
            jet2P4 = jet2.p4()
            jet2P4.Boost(-bjjP4.BoostVector())
            
            return max([bjetP4.Angle(jet1P4.Vect()),bjetP4.Angle(jet2P4.Vect()),jet1P4.Angle(jet1P4.Vect())])
        
        self.variableDict = {
            #event vars
            'rho': lambda event,bjet,jet1,jet2: event.fixedGridRhoFastjetAll,
            'npv': lambda event,bjet,jet1,jet2: event.PV_npvsGood,
            
            
            #bjet vars
            'bjet_area': lambda event,bjet,jet1,jet2: bjet.area,
            'bjet_btagDeepFlavB': lambda event,bjet,jet1,jet2: bjet.btagDeepFlavB,
            #'bjet_btagDeepFlavCvB':lambda event,bjet,jet1,jet2: bjet.btagDeepFlavCvB,
            #'bjet_btagDeepFlavCvL':lambda event,bjet,jet1,jet2: bjet.btagDeepFlavCvL,
            #'bjet_btagDeepFlavQG':lambda event,bjet,jet1,jet2: bjet.btagDeepFlavQG,
            'bjet_qgl': lambda event,bjet,jet1,jet2: bjet.qgl,
            'bjet_nConstituents': lambda event,bjet,jet1,jet2: bjet.nConstituents,
            
            'bjet_chEmEF': lambda event,bjet,jet1,jet2: bjet.chEmEF,
            'bjet_chFPV0EF': lambda event,bjet,jet1,jet2: bjet.chFPV0EF,
            
            'bjet_chHEF': lambda event,bjet,jet1,jet2: bjet.chHEF,
            
            'bjet_neEmEF': lambda event,bjet,jet1,jet2: bjet.neEmEF,
            'bjet_neHEF': lambda event,bjet,jet1,jet2: bjet.neHEF,
            
            
            'bjet_pt': lambda event,bjet,jet1,jet2: bjet.pt,
            'bjet_abseta': lambda event,bjet,jet1,jet2: math.fabs(bjet.eta),
            'bjet_mass': lambda event,bjet,jet1,jet2: bjet.mass,
            
            
            #jet vars
            'jet_area_max': lambda event,bjet,jet1,jet2: max(jet1.area, jet2.area),
            'jet_btagDeepFlavB_max': lambda event,bjet,jet1,jet2: max(jet1.btagDeepFlavB, jet2.btagDeepFlavB),
            #'jet_btagDeepFlavCvB_max':lambda event,bjet,jet1,jet2: max(jet1.btagDeepFlavCvB,jet2.btagDeepFlavCvB),
            #'jet_btagDeepFlavCvL_max':lambda event,bjet,jet1,jet2: max(jet1.btagDeepFlavCvL,jet2.btagDeepFlavCvL),
            #'jet_btagDeepFlavQG_max':lambda event,bjet,jet1,jet2: bjet.btagDeepFlavQG,
            
            
            'jet_qgl_max': lambda event,bjet,jet1,jet2: max(jet1.qgl, jet2.qgl),
            'jet_qgl_min': lambda event,bjet,jet1,jet2: min(jet1.qgl, jet2.qgl),
            'jet_nConstituents_max': lambda event,bjet,jet1,jet2: max(jet1.nConstituents, jet2.nConstituents),
            'jet_nConstituents_min': lambda event,bjet,jet1,jet2: min(jet1.nConstituents, jet2.nConstituents),
            
            'jet_chEmEF_max': lambda event,bjet,jet1,jet2: max(jet1.chEmEF, jet2.chEmEF),
            'jet_chFPV0EF_max': lambda event,bjet,jet1,jet2: max(jet1.chFPV0EF, jet2.chFPV0EF),
            
            'jet_chHEF_max': lambda event,bjet,jet1,jet2: max(jet1.chHEF,jet2.chHEF),
            
            'jet_neEmEF_max': lambda event,bjet,jet1,jet2: max(jet1.neEmEF,jet2.neEmEF),
            'jet_neHEF_max': lambda event,bjet,jet1,jet2: max(jet1.neHEF,jet2.neEmEF),
            
            'jet1_pt': lambda event,bjet,jet1,jet2: jet1.pt,
            'jet1_abseta': lambda event,bjet,jet1,jet2: math.fabs(jet1.eta),
            'jet1_mass': lambda event,bjet,jet1,jet2: jet1.mass,
            
            'jet2_pt': lambda event,bjet,jet1,jet2: jet2.pt,
            'jet2_abseta': lambda event,bjet,jet1,jet2: math.fabs(jet2.eta),
            'jet2_mass': lambda event,bjet,jet1,jet2: jet2.mass,
            
            
            #jj correlations
            'jj_ptrel': lambda event,bjet,jet1,jet2: jet1.pt/jet2.pt,
            
            'jj_dR': lambda event,bjet,jet1,jet2: deltaR(jet1,jet2),
            'jj_dPhi': lambda event,bjet,jet1,jet2: deltaPhi(jet1,jet2),
            'jj_dEta': lambda event,bjet,jet1,jet2: math.fabs(jet1.eta-jet2.eta),
            
            'jj_mass': lambda event,bjet,jet1,jet2: (jet1.p4()+jet2.p4()).M(),
            'jj_pt': lambda event,bjet,jet1,jet2: (jet1.p4()+jet2.p4()).Pt(),
            
            
            #bj correlations
            'bj_ptrel_min': lambda event,bjet,jet1,jet2: min(jet1.pt/bjet.pt,jet2.pt/bjet.pt),
            'bj_ptrel_max': lambda event,bjet,jet1,jet2: max(jet1.pt/bjet.pt,jet2.pt/bjet.pt),
            
            'bj_dR_min': lambda event,bjet,jet1,jet2: min(deltaR(bjet,jet1),deltaR(bjet,jet2)),
            'bj_dPhi_min': lambda event,bjet,jet1,jet2: min(deltaPhi(bjet,jet1),deltaPhi(bjet,jet2)),
            'bj_dEta_min': lambda event,bjet,jet1,jet2: min(math.fabs(bjet.eta-jet1.eta),math.fabs(bjet.eta-jet2.eta)),
            
            'bj_dR_max': lambda event,bjet,jet1,jet2: max(deltaR(bjet,jet1),deltaR(bjet,jet2)),
            'bj_dPhi_max': lambda event,bjet,jet1,jet2: max(deltaPhi(bjet,jet1),deltaPhi(bjet,jet2)),
            'bj_dEta_max': lambda event,bjet,jet1,jet2: max(math.fabs(bjet.eta-jet1.eta),math.fabs(bjet.eta-jet2.eta)),
            
            'bj_mass_min': lambda event,bjet,jet1,jet2: min((bjet.p4()+jet1.p4()).M(),(bjet.p4()+jet2.p4()).M()),
            'bj_mass_max': lambda event,bjet,jet1,jet2: max((bjet.p4()+jet1.p4()).M(),(bjet.p4()+jet2.p4()).M()),
            
            'bj_pt_min': lambda event,bjet,jet1,jet2: min((bjet.p4()+jet1.p4()).Pt(),(bjet.p4()+jet2.p4()).Pt()),
            'bj_pt_max': lambda event,bjet,jet1,jet2: max((bjet.p4()+jet1.p4()).Pt(),(bjet.p4()+jet2.p4()).Pt()),
            
            
            #bjj correlations
            'bjj_mass': lambda event,bjet,jet1,jet2: (bjet.p4()+jet1.p4()+jet2.p4()).M(),
            'bjj_pt': lambda event,bjet,jet1,jet2: (bjet.p4()+jet1.p4()+jet2.p4()).Pt(),
            'bjj_min_angle': lambda event,bjet,jet1,jet2: minAngleInBJJ(bjet,jet1,jet2),
            'bjj_max_angle': lambda event,bjet,jet1,jet2: maxAngleInBJJ(bjet,jet1,jet2)
        }
    
    def bookBranches(self,out):
        out.branch(self.getLenVarName(), "I")
        for varName in self.variableDict.keys():
            out.branch(self.outputName+'_'+varName, "F", lenVar=self.getLenVarName())
            
            
    def getLenVarName(self):
        return 'n'+self.outputName
        
    def queryValues(self,event,jetCombinations):
        valueDict = {
            varName: np.zeros(len(jetCombinations),dtype=np.float32) for varName in self.variableDict.keys()
        }
        
        for icomb, (bjet, jet1, jet2) in enumerate(jetCombinations):
            for varName, varFct in self.variableDict.items():
                valueDict[varName][icomb] = varFct(event, bjet, jet1, jet2)
                
        return valueDict
        
    def writeBranches(self,out,valueDict):
        out.fillBranch('n'+self.outputName, len(valueDict.values()[0]))
        for varName in self.variableDict.keys():
            out.fillBranch(self.outputName+'_'+varName,valueDict[varName])


class TopNNRecoInputs(Module):
    def __init__(self,
        bjetCollection,
        jetCollection,
        outputName,
    ):
        self.bjetCollection = bjetCollection
        self.jetCollection = jetCollection
        self.outputName=outputName
        
        self.hadronicTopInputs = NNHadronicTopInputs(self.outputName)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.hadronicTopInputs.bookBranches(self.out)
        
        self.out.branch(self.outputName+'_gen_dRMatch', "F", lenVar=self.hadronicTopInputs.getLenVarName())
        self.out.branch(self.outputName+'_gen_dRTopMatch', "F", lenVar=self.hadronicTopInputs.getLenVarName())
        self.out.branch(self.outputName+'_gen_dRBJetMatch', "F", lenVar=self.hadronicTopInputs.getLenVarName())
        self.out.branch(self.outputName+'_gen_dRJet1Match', "F", lenVar=self.hadronicTopInputs.getLenVarName())
        self.out.branch(self.outputName+'_gen_dRJet2Match', "F", lenVar=self.hadronicTopInputs.getLenVarName())
        
        self.out.branch(self.outputName+'_gen_bjetPtMatch', "F", lenVar=self.hadronicTopInputs.getLenVarName())
        
        self.out.branch(self.outputName+'_gen_topPt', "F", lenVar=self.hadronicTopInputs.getLenVarName())

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
                top['quarks'] = sorted(top['quarks'],key=lambda x: x.pt, reverse=True) #sort by pT like jets
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
        
    def makeJetCombinations(self, bjets, jets):
        combinations = []
        
        # sort by pT so that pt(jet1)>pt(jet2)
        jets = sorted(jets,key=lambda x: x.pt,reverse=True)
        
        for bjet in bjets:
            for i,jet1 in enumerate(jets):
                for j,jet2 in enumerate(jets[i+1:]):
                
                    # remove potential double counting if bjets/jets are not mutually exclusive
                    if deltaR(bjet,jet1)<0.4:
                        continue
                    if deltaR(bjet,jet2)<0.4:
                        continue
                        
                    bjjMass = (bjet.p4()+jet1.p4()+jet2.p4()).M()
                    if bjjMass<100. or bjjMass>300.:
                        continue
                        
                    combinations.append([bjet,jet1,jet2])
                    
        #optional only keep up to 10 closest combinations
        #combinations = sorted(combinations,key = lambda x: math.fabs(172.5-(x[0].p4()+x[1].p4()+x[2].p4()).M()))
        #combinations = combinations[:10]
        
        return combinations

    def analyze(self, event):
        
        bjets = self.bjetCollection(event)
        jets = self.jetCollection(event)
        
        jetCombinations = self.makeJetCombinations(bjets,jets)
        
        
        valueDict = self.hadronicTopInputs.queryValues(event,jetCombinations)
        
        genHadronicTopDecays = self.getGenLevelTopDecays(event)["hadronic"]
        
        gen_dRMatch = 10.*np.ones(len(jetCombinations),dtype=np.float32)
        
        gen_dRTopMatch = 10.*np.ones(len(jetCombinations),dtype=np.float32)
        gen_dRBJetMatch = 10.*np.ones(len(jetCombinations),dtype=np.float32)
        gen_dRJet1Match = 10.*np.ones(len(jetCombinations),dtype=np.float32)
        gen_dRJet2Match = 10.*np.ones(len(jetCombinations),dtype=np.float32)
        
        
        gen_bjetPtMatch = np.zeros(len(jetCombinations),dtype=np.float32)
        gen_topPt = np.zeros(len(jetCombinations),dtype=np.float32)
        
        for icomb, (bjet, jet1, jet2) in enumerate(jetCombinations):
            for genHadronicTopDecay in genHadronicTopDecays:
                genBquark = genHadronicTopDecay['bquark'][0]
                genQuark1 = genHadronicTopDecay['quarks'][0]
                genQuark2 = genHadronicTopDecay['quarks'][1]
                genTopQuark =  genHadronicTopDecay['top']
                
                dRbjet = deltaR(bjet,genBquark)
                if deltaR(jet1,genQuark1)<deltaR(jet1,genQuark2):
                    dRjet1 = deltaR(jet1,genQuark1)
                    dRjet2 = deltaR(jet2,genQuark2)
                else:
                    dRjet1 = deltaR(jet1,genQuark2)
                    dRjet2 = deltaR(jet2,genQuark1)
                    

                dRMatch = math.sqrt((dRbjet)**2+(dRjet1)**2+(dRjet2)**2)
                
                if dRMatch<gen_dRMatch[icomb]:
                    gen_dRMatch[icomb] = dRMatch
                    gen_dRBJetMatch[icomb] = dRbjet
                    gen_dRJet1Match[icomb] = dRjet1
                    gen_dRJet2Match[icomb] = dRjet2
                    
                    gen_bjetPtMatch[icomb] = bjet.pt/genBquark.pt
                    gen_topPt[icomb] = genTopQuark.pt
                    
                    topP4 = (bjet.p4()+jet1.p4()+jet2.p4())
                    gen_dRTopMatch[icomb] = topP4.DeltaR(genTopQuark.p4())
                    
        self.hadronicTopInputs.writeBranches(self.out,valueDict)
        
        self.out.fillBranch(self.outputName+'_gen_dRMatch', gen_dRMatch)
        self.out.fillBranch(self.outputName+'_gen_dRTopMatch', gen_dRTopMatch)
        
        self.out.fillBranch(self.outputName+'_gen_dRBJetMatch', gen_dRBJetMatch)
        self.out.fillBranch(self.outputName+'_gen_dRJet1Match', gen_dRJet1Match)
        self.out.fillBranch(self.outputName+'_gen_dRJet2Match', gen_dRJet2Match)
        
        self.out.fillBranch(self.outputName+'_gen_bjetPtMatch', gen_bjetPtMatch)
        self.out.fillBranch(self.outputName+'_gen_topPt', gen_topPt)
        
        return True
        
        
