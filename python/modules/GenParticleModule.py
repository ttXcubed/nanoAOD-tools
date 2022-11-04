import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from collections import OrderedDict

from gen_helper import *
from utils import deltaR, deltaPhi

class GenParticleModule(Module):

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "GenPart"),
        inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
        inputGenJetCollection=lambda event: Collection(event, "GenJet"),
        inputFatJetCollection=lambda event: Collection(event, "Jet"),
        inputJetCollection=lambda event: Collection(event, "Jet"),
        inputMuonCollection=lambda event: Collection(event, "Muon"),
        inputElectronCollection=lambda event: Collection(event, "Electron"),
        outputName="genPart",
        storeKinematics= ['pt','eta','phi','mass'],
    ):
        
        self.inputCollection = inputCollection
        self.inputFatGenJetCollection = inputFatGenJetCollection
        self.inputGenJetCollection = inputGenJetCollection
        self.inputFatJetCollection = inputFatJetCollection
        self.inputJetCollection = inputJetCollection
        self.inputMuonCollection = inputMuonCollection
        self.inputElectronCollection = inputElectronCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics
		
    def genP_genJet_matching(self, ptGenP, etaGenP, phiGenP, gJets, minDeltaR, minRelDeltaPt):
    
    	matched_gJet_idx = -99    
    	
    	for i_gJet, gJet in enumerate(gJets):
    		min_deltaR =  min(minDeltaR, math.sqrt((gJet.eta - etaGenP)**2 + (self.deltaPhi(gJet.phi,phiGenP))**2))
		min_relDeltaPt = min(minRelDeltaPt, self.relDeltaPt(gJet.pt, ptGenP))
		
		if min_deltaR < minDeltaR and min_relDeltaPt < minRelDeltaPt: 
    			matched_gJet_idx = i_gJet 
    		else: continue
    	
    	return matched_gJet_idx
    	
    def deltaPhi(self,phi1, phi2):
	res = phi1-phi2
	while (res > math.pi):
		res -= 2 * math.pi
	while (res <= -math.pi):
		res += 2 * math.pi

	return res

    def relDeltaPt(self,pt1, pt2):
	if pt2 != 0.0:
		return abs((pt1-pt2)/pt2)
 	return 1.
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("n"+self.outputName+"_selectedJets", "I")
	self.out.branch("n"+self.outputName+"_selectedFatJet", "I")
        
        self.out.branch("selectedFatJets_flavor","I",lenVar="n"+self.outputName+"_selectedFatJet") 
        self.out.branch("selectedJets_flavor","I",lenVar="n"+self.outputName+"_selectedJets") 	
	

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        genParticles = self.inputCollection(event)
        fatGenJets = self.inputFatGenJetCollection(event)
        genJets = self.inputGenJetCollection(event)
        fatGenJets = sorted(fatGenJets,key=lambda x: x.pt, reverse=True)
        genJets = sorted(genJets,key=lambda x: x.pt, reverse=True)
        
        muons = self.inputMuonCollection(event)
        electrons = self.inputElectronCollection(event)
        
        fatjets = self.inputFatJetCollection(event)
        jets = self.inputJetCollection(event)
 
	genP_firstcopy_pdgid, genP_firstcopy_id,genP_firstcopy_mother_idx = [], [], []
	

 	genJetIdx_list, genFatJetIdx_list = [], []
 	genP_matchedToJet_pdgId_list, genP_matchedToFatJet_pdgId_list = [], []
 	
 	for igenP, genParticle in enumerate(genParticles):
 		if isFirstCopy(genParticle) and isPrompt(genParticle) and isHardProcess(genParticle) and fromHardProcess(genParticle):
 			#print("genpart pdg ", genParticle.pdgId)
 			genP_firstcopy_id.append(igenP)
 			genP_firstcopy_pdgid.append(genParticle.pdgId)
 			genP_firstcopy_mother_idx.append(genParticle.genPartIdxMother)
 			
 			gFJets_idx = self.genP_genJet_matching(genParticle.pt, genParticle.eta, genParticle.phi, fatGenJets, 0.8, 0.5)
 			#print("gfjet idx ", gFJets_idx)
			if gFJets_idx != -99 and gFJets_idx not in genFatJetIdx_list: 
				genFatJetIdx_list.append(gFJets_idx)
				#genFatJetIdx_list.append(igenP)
				genP_matchedToFatJet_pdgId_list.append(genParticle.pdgId)
						
			gJets_idx = self.genP_genJet_matching(genParticle.pt, genParticle.eta, genParticle.phi, genJets, 0.4, 0.5) 
			if gJets_idx != -99 and gJets_idx not in genJetIdx_list: 
				genJetIdx_list.append(gJets_idx)	
				genP_matchedToJet_pdgId_list.append(genParticle.pdgId)			
 		else:
 			continue

	#print("jet ", genJetIdx_list)
	#print("fj ", genFatJetIdx_list)
	
	#print("gen P fj ", genP_matchedToFatJet_pdgId_list)
	#print("gen P j ", genP_matchedToJet_pdgId_list)

	'''
	for mu in muons:
		#print('mu charge :', mu.charge)
		#print('mu genPartFlav:', mu.genPartFlav)
		for ig, gp in enumerate(genP_firstcopy_id): 
			if mu.genPartIdx==gp:
				#print('mu genPart asso: ', genParticles[mu.genPartIdx].pdgId)
				#print('madre mu genPart asso: ', genParticles[genP_firstcopy_mother_idx[ig]].pdgId)
				#print('madre idx mu genPart asso: ', genParticles[genP_firstcopy_mother_idx[ig]])
				for i, ak8 in enumerate(fatjets):
					if ak8.genJetAK8Idx in genFatJetIdx_list:
						#print(genP_matchedToFatJet_pdgId_list[i])
						#print(math.sqrt((mu.eta - fatjets[i].eta)**2 + (self.deltaPhi(mu.phi,fatjets[i].phi))**2))
					else: continue

	for e in electrons:
		#print('e charge :', e.charge)
		#print('e genPartFlav:', e.genPartFlav)
		#print('e genPart asso: ', genParticles[e.genPartIdx].pdgId) 
		for i, ak8 in enumerate(fatjets):
					if ak8.genJetAK8Idx in genFatJetIdx_list:
						#print(genP_matchedToFatJet_pdgId_list[i])
						#print(math.sqrt((e.eta - fatjets[i].eta)**2 + (self.deltaPhi(e.phi,fatjets[i].phi))**2))
					else: continue
	'''

	jet_flavor_list, fatjet_flavor_list = [],[]

	for iak8, ak8 in enumerate(fatjets):
		
		if ak8.genJetAK8Idx in genFatJetIdx_list:
			#print('gen fjet pt {} eta {} phi {}'.format(fatGenJets[ak8.genJetAK8Idx].pt, fatGenJets[ak8.genJetAK8Idx].eta, fatGenJets[ak8.genJetAK8Idx].phi))
			#print('fjet pt {} eta {} phi {}'.format(ak8.pt, ak8.eta, ak8.phi))
			#print('fjet tvsqcd {} '.format(ak8.particleNet_TvsQCD))
			fatjet_flavor_list.append(genP_matchedToFatJet_pdgId_list[genFatJetIdx_list.index(ak8.genJetAK8Idx)])
		else: 
			fatjet_flavor_list.append(-99)
	
	for ak4 in jets:
		if ak4.genJetIdx in genJetIdx_list:
			jet_flavor_list.append(genP_matchedToJet_pdgId_list[genJetIdx_list.index(ak4.genJetIdx)])
		else:
			jet_flavor_list.append(-99)
	
	self.out.fillBranch("n"+self.outputName+"_selectedFatJet", len(fatjets))
	self.out.fillBranch("n"+self.outputName+"_selectedJets", len(jets))

	self.out.fillBranch("selectedFatJets_flavor", fatjet_flavor_list)
	self.out.fillBranch("selectedJets_flavor", jet_flavor_list)
	
		          
        return True

