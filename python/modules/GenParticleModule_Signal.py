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

class GenParticleModule_Signal(Module):

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "GenPart"),
        inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
        inputGenJetCollection=lambda event: Collection(event, "GenJet"),
        inputFatJetCollection=lambda event: Collection(event, "Jet"),
        inputJetCollection=lambda event: Collection(event, "Jet"),
        #eventReco_flags = [],
        outputName="genPart",
        storeKinematics= ['pt','eta','phi','mass'],
    ):
        
        self.inputCollection = inputCollection
        self.inputFatGenJetCollection = inputFatGenJetCollection
        self.inputGenJetCollection = inputGenJetCollection
        self.inputFatJetCollection = inputFatJetCollection
        self.inputJetCollection = inputJetCollection
        #self.eventReco_flags = eventReco_flags
        self.outputName = outputName
        self.storeKinematics = storeKinematics
     
    def getFirstCopy(self,motherIdx,pdg,index):
    	if old_pdg_list[motherIdx] == pdg and motherIdx != -1: 
		index = motherIdx
		return self.getFirstCopy(old_motherIdx_list[motherIdx],pdg,index)
    	elif old_pdg_list[motherIdx] != pdg and motherIdx != -1:
		return index 
		
    def genP_genJet_matching(self, genP_idx, kinematicsGenP, gJets, minDeltaR, minRelDeltaPt):
    
    	matched_gJet_idx = -99    
    	
    	for i_gJet, gJet in enumerate(gJets):
    		min_deltaR =  min(minDeltaR, math.sqrt((gJet.eta - kinematicsGenP[1][genP_idx])**2 + (self.deltaPhi(gJet.phi,kinematicsGenP[2][genP_idx]))**2))
		min_relDeltaPt = min(minRelDeltaPt, self.relDeltaPt(gJet.pt, kinematicsGenP[0][genP_idx]))
		
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
	return abs((pt1-pt2)/pt2)
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")
        self.out.branch("n"+self.outputName+"_selectedJets", "I")
	self.out.branch("n"+self.outputName+"_selectedFatJet", "I")
	
        #self.out.branch("n"+self.outputName+"_diHadronic_cleanedAK4s", "I")
	#self.out.branch("n"+self.outputName+"_diHadronic_FatJet", "I")
	
        self.out.branch(self.outputName+"_pdgId","I",lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_IdxMother","I",lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_status","I",lenVar="n"+self.outputName)  #1=stable

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
        
        self.out.branch("selectedFatJet_spectatorFlag_topMatching","I",lenVar="n"+self.outputName+"_selectedFatJet") 
	self.out.branch("selectedFatJet_resonantFlag_topMatching","I",lenVar="n"+self.outputName+"_selectedFatJet")
	self.out.branch("selectedFatJet_spectatorFlag_wMatching","I",lenVar="n"+self.outputName+"_selectedFatJet") 
	self.out.branch("selectedFatJet_resonantFlag_wMatching","I",lenVar="n"+self.outputName+"_selectedFatJet")
	
	self.out.branch("selectedJets_spectatorFlag_topMatching","I",lenVar="n"+self.outputName+"_selectedJets")
	self.out.branch("selectedJets_resonantFlag_topMatching","I",lenVar="n"+self.outputName+"_selectedJets") 
	self.out.branch("selectedJets_spectatorFlag_bMatching","I",lenVar="n"+self.outputName+"_selectedJets")
	self.out.branch("selectedJets_resonantFlag_bMatching","I",lenVar="n"+self.outputName+"_selectedJets") 
        
	#self.out.branch("diHadronic_FatJet_spectatorFlag","I",lenVar="n"+self.outputName+"_diHadronic_FatJet") 
	#self.out.branch("diHadronic_FatJet_resonantFlag","I",lenVar="n"+self.outputName+"_diHadronic_FatJet")
	
	#self.out.branch("diHadronic_cleanedAK4s_spectatorFlag_topMatching","I",lenVar="n"+self.outputName+"_diHadronic_cleanedAK4s")
	#self.out.branch("diHadronic_cleanedAK4s_resonantFlag_topMatching","I",lenVar="n"+self.outputName+"_diHadronic_cleanedAK4s") 
	#self.out.branch("diHadronic_cleanedAK4s_spectatorFlag_bMatching","I",lenVar="n"+self.outputName+"_diHadronic_cleanedAK4s")
	#self.out.branch("diHadronic_cleanedAK4s_resonantFlag_bMatching","I",lenVar="n"+self.outputName+"_diHadronic_cleanedAK4s") 
	

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        genParticles = self.inputCollection(event)
        fatGenJets = self.inputFatGenJetCollection(event)
        genJets = self.inputGenJetCollection(event)
        
        fatjets = self.inputFatJetCollection(event)
        #fatJets_dict = self.inputFatJetCollection(event)
        #print("gen  ", len(fatJets_dict[0:2])
        jets = self.inputJetCollection(event)
              
        #eventReco_flag_dict = self.eventReco_flags(event)
 
        fatGenJets = sorted(fatGenJets,key=lambda x: x.pt, reverse=True)
        genJets = sorted(genJets,key=lambda x: x.pt, reverse=True)
       
        oldToNew = OrderedDict()
        
        pdg_list = []
        motherIdx_list = []
        
        status_list = []
        statusFlag_list = []
        variable_list = []
        
        global old_motherIdx_list
        old_motherIdx_list = []
        global old_pdg_list
        old_pdg_list = []
        
        selected_pdg_list = [1,2,3,4,5,6,11,12,13,14,15,16,24,6000055]
        
        i=0
        while i < len(self.storeKinematics):
        	variable_list.append([])
        	i+=1	
        
	for igenP, genParticle in enumerate(genParticles):
		old_motherIdx_list.append(genParticle.genPartIdxMother)
		old_pdg_list.append(genParticle.pdgId)

		
		if genParticle.genPartIdxMother>igenP :  #CHECK NECESSARY!!!!!  igenP==0 or igenP==1  (genParticle.genPartIdxMother !=-1) 
			break	
	
		else:
			firstCopy_idx = self.getFirstCopy(genParticle.genPartIdxMother,genParticle.pdgId,igenP)			
			if firstCopy_idx in oldToNew.keys() or firstCopy_idx is None:
        			#print("index already in the dict")
        			continue
        		else:
        			particleFirstCopy = genParticles[firstCopy_idx]
        			partPdg = particleFirstCopy.pdgId
				if abs(partPdg) in selected_pdg_list:
				
				    if abs(partPdg) == 6 or abs(partPdg) == 6000055: 
					#print("X or t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    							
				    if abs(partPdg) == 5 and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 6 and isPrompt(particleFirstCopy) and isHardProcess(particleFirstCopy) and fromHardProcess(particleFirstCopy) and isFirstCopy(particleFirstCopy):	    
				    
				    	# and particleFirstCopy.statusFlags!=8193 and particleFirstCopy.statusFlags!=257 and particleFirstCopy.statusFlags!=1 and particleFirstCopy.statusFlags!=20481 and particleFirstCopy.statusFlags!=8449 and particleFirstCopy.statusFlags!=16385 and particleFirstCopy.statusFlags!=12289: 
					#print("b from t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    			
				    if (abs(partPdg) == 1 or abs(partPdg) == 2 or abs(partPdg) == 3 or abs(partPdg) == 4 or abs(partPdg) == 5 or abs(partPdg) == 11 or abs(partPdg) == 12 or abs(partPdg) == 13 or abs(partPdg) == 14 or abs(partPdg) == 15 or abs(partPdg) == 16) and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 24 and abs(old_pdg_list[old_motherIdx_list[self.getFirstCopy(old_motherIdx_list[firstCopy_idx],old_pdg_list[old_motherIdx_list[firstCopy_idx]],firstCopy_idx)]]) == 6 and isPrompt(particleFirstCopy) and isHardProcess(particleFirstCopy) and fromHardProcess(particleFirstCopy) and isFirstCopy(particleFirstCopy):
	
				    	#and particleFirstCopy.statusFlags!=12289 and particleFirstCopy.statusFlags!=8193  and particleFirstCopy.statusFlags!=12291 and particleFirstCopy.statusFlags!=4097 and particleFirstCopy.statusFlags!=20481:  
				    	pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    			
				    if abs(partPdg) == 24 and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 6:
					#print("W from t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
				    else: continue
				else:
				    continue
	
				
	for oldIdx in oldToNew.keys():
    		oldMotherIdx = self.getFirstCopy(old_motherIdx_list[oldIdx],old_pdg_list[old_motherIdx_list[oldIdx]],oldIdx)
    		if oldMotherIdx not in oldToNew.keys():
        		newMotherIdx = -99
        		motherIdx_list.append(newMotherIdx)
    		else:
        		newMotherIdx = oldToNew[oldMotherIdx]
        		motherIdx_list.append(newMotherIdx)
	
	
	genFatJetIdx_top_list_resonant, genFatJetIdx_top_list_spectator = [], []
	genFatJetIdx_w_list_resonant, genFatJetIdx_w_list_spectator = [], [] 
	genJetIdx_top_list_resonant, genJetIdx_top_list_spectator = [], [] #### idx of gen ak4 linked to top part
	genJetIdx_b_list_resonant, genJetIdx_b_list_spectator = [], [] #### idx of gen ak4 linked to b part
	
	genMatching_jet_top_spectator_flag, genMatching_jet_top_resonant_flag = [], [] 
	genMatching_jet_b_spectator_flag, genMatching_jet_b_resonant_flag = [], [] 
	genMatching_fatjet_top_spectator_flag, genMatching_fatjet_top_resonant_flag = [], []
	genMatching_fatjet_w_spectator_flag, genMatching_fatjet_w_resonant_flag = [], []

	for iP, particle in enumerate(pdg_list):
				
		if abs(particle) == 6:
			gFJets_idx = self.genP_genJet_matching(iP, variable_list, fatGenJets, 0.8, 0.5)
			gJets_idx = self.genP_genJet_matching(iP, variable_list, genJets, 0.4, 0.5)
		
			if gJets_idx != -99:	
				if motherIdx_list[iP] == -99:		
					genJetIdx_top_list_spectator.append(gJets_idx)
				elif abs(pdg_list[motherIdx_list[iP]]) == 6000055:
					genJetIdx_top_list_resonant.append(gJets_idx)
										
			if gFJets_idx != -99: 
				if motherIdx_list[iP] == -99:		
					genFatJetIdx_top_list_spectator.append(gFJets_idx)
				elif abs(pdg_list[motherIdx_list[iP]]) == 6000055:
					genFatJetIdx_top_list_resonant.append(gFJets_idx)
																
		elif abs(particle) == 24:
			gJets_idx = self.genP_genJet_matching(iP, variable_list, fatGenJets, 0.8, 0.5)
			if gJets_idx != -99:	
				if motherIdx_list[motherIdx_list[iP]] == -99:
					genFatJetIdx_w_list_spectator.append(gJets_idx)		
				elif abs(pdg_list[motherIdx_list[motherIdx_list[iP]]]) == 6000055:
					genFatJetIdx_w_list_resonant.append(gJets_idx)
				
		elif abs(particle) == 5:
			gJets_idx = self.genP_genJet_matching(iP, variable_list, genJets, 0.4, 0.5)
			if gJets_idx != -99:	
				if motherIdx_list[motherIdx_list[iP]] == -99:	
					genJetIdx_b_list_spectator.append(gJets_idx)	
				elif abs(pdg_list[motherIdx_list[motherIdx_list[iP]]]) == 6000055:
			 		genJetIdx_b_list_resonant.append(gJets_idx)
		else: continue					
			
	#print("EVENTO")
	
	#print("LEN FJ ", len(fatjets))
	#print("fj reso {} spec {} top".format(genFatJetIdx_top_list_resonant, genFatJetIdx_top_list_spectator))
	#print("fj reso {} spec {} w".format(genFatJetIdx_w_list_resonant, genFatJetIdx_w_list_spectator))
	#print("j reso {} spect {} top".format(genJetIdx_top_list_resonant, genJetIdx_top_list_spectator))
	#print("j reso {} spect {} b".format(genJetIdx_b_list_resonant, genJetIdx_b_list_spectator))
	
			
	for ak8 in fatjets: #fatJets_dict[eventReco_flag][0:2]:	
		#print("gen jet idx ", ak8.genJetAK8Idx)
		if ak8.genJetAK8Idx in genFatJetIdx_top_list_spectator:
			genMatching_fatjet_top_spectator_flag.append(True)
			genMatching_fatjet_top_resonant_flag.append(False)
		elif ak8.genJetAK8Idx in genFatJetIdx_top_list_resonant:	
			genMatching_fatjet_top_resonant_flag.append(True)
			genMatching_fatjet_top_spectator_flag.append(False)
		else: 
			genMatching_fatjet_top_spectator_flag.append(False)
			genMatching_fatjet_top_resonant_flag.append(False)

	for ak8 in fatjets:
		if ak8.genJetAK8Idx in genFatJetIdx_w_list_spectator:
			genMatching_fatjet_w_spectator_flag.append(True)
			genMatching_fatjet_w_resonant_flag.append(False)
		elif ak8.genJetAK8Idx in genFatJetIdx_w_list_resonant:	
			genMatching_fatjet_w_resonant_flag.append(True)
			genMatching_fatjet_w_spectator_flag.append(False)
		else: 
			genMatching_fatjet_w_spectator_flag.append(False)
			genMatching_fatjet_w_resonant_flag.append(False)
	
	for ak4 in jets:
		if ak4.genJetIdx in genJetIdx_top_list_spectator:
			genMatching_jet_top_spectator_flag.append(True)
			genMatching_jet_top_resonant_flag.append(False)				
		elif ak4.genJetIdx in genJetIdx_top_list_resonant:
			genMatching_jet_top_resonant_flag.append(True)
			genMatching_jet_top_spectator_flag.append(False)
		else: 
			genMatching_jet_top_spectator_flag.append(False)
			genMatching_jet_top_resonant_flag.append(False)
			
	for ak4 in jets:
		#print("gen jet idx ", ak4.genJetIdx)	
		if ak4.genJetIdx in genJetIdx_b_list_spectator:
			genMatching_jet_b_spectator_flag.append(True)
			genMatching_jet_b_resonant_flag.append(False)				
		elif ak4.genJetIdx in genJetIdx_b_list_resonant:
			genMatching_jet_b_resonant_flag.append(True)
			genMatching_jet_b_spectator_flag.append(False)
		else: 
			genMatching_jet_b_spectator_flag.append(False)
			genMatching_jet_b_resonant_flag.append(False)	
		

	#print("fj reso {} spce {} FLAG TOP".format(genMatching_fatjet_top_resonant_flag, genMatching_fatjet_top_spectator_flag))
	#print("fj reso {} spce {} FLAG W".format(genMatching_fatjet_w_resonant_flag, genMatching_fatjet_w_spectator_flag))
	#print("j reso {} spce {} FLAG TOP".format(genMatching_jet_top_resonant_flag, genMatching_jet_top_spectator_flag))
	#print("j reso {} spce {} FLAG B".format(genMatching_jet_b_resonant_flag, genMatching_jet_b_spectator_flag))
	
	self.out.fillBranch("n"+self.outputName+"_selectedFatJet", len(fatjets))
	self.out.fillBranch("n"+self.outputName+"_selectedJets", len(jets))
	self.out.fillBranch("n"+self.outputName,len(pdg_list))
			
	for i,variable in enumerate(self.storeKinematics):
		self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) 
	
	self.out.fillBranch(self.outputName+"_pdgId", pdg_list) 
	self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) 
	self.out.fillBranch(self.outputName+"_status", status_list) 
	

	self.out.fillBranch("selectedFatJet_spectatorFlag_topMatching", genMatching_fatjet_top_spectator_flag)
	self.out.fillBranch("selectedFatJet_resonantFlag_topMatching", genMatching_fatjet_top_resonant_flag)
	self.out.fillBranch("selectedFatJet_spectatorFlag_wMatching", genMatching_fatjet_w_spectator_flag)
	self.out.fillBranch("selectedFatJet_resonantFlag_wMatching", genMatching_fatjet_w_resonant_flag)
	
	self.out.fillBranch("selectedJets_spectatorFlag_topMatching", genMatching_jet_top_spectator_flag)
	self.out.fillBranch("selectedJets_resonantFlag_topMatching", genMatching_jet_top_resonant_flag)
	self.out.fillBranch("selectedJets_spectatorFlag_bMatching", genMatching_jet_b_spectator_flag)
	self.out.fillBranch("selectedJets_resonantFlag_bMatching", genMatching_jet_b_resonant_flag)	

	'''
	for eventReco_flag in eventReco_flag_dict.keys():
		genMatching_jet_top_spectator_flag, genMatching_jet_top_resonant_flag = [], [] 
		genMatching_jet_b_spectator_flag, genMatching_jet_b_resonant_flag = [], [] 
		genMatching_fatjet_top_spectator_flag, genMatching_fatjet_top_resonant_flag = [], []
		genMatching_fatjet_w_spectator_flag, genMatching_fatjet_w_resonant_flag = [], []

		if eventReco_flag_dict[eventReco_flag]:	
			
			for iP, particle in enumerate(pdg_list):
				
					if abs(particle) == 6:
						gFJets_idx = self.genP_genJet_matching(genP_idx, variable_list, fatGenJets, 1.6, 0.5)
						gJets_idx = self.genP_genJet_matching(genP_idx, variable_list, genJets, 0.8, 0.5)
					
						if gJets_index != -99:	
							if motherIdx_list[motherIdx_list[iP]] == -99:		
								genJetIdx_top_list_spectator.append(gJets_idx)
							elif abs(pdg_list[motherIdx_list[motherIdx_list[iP]]) == 6000055:
								genFatJetIdx_top_list_resonant.append(gJets_idx)
													
						if gFJets_index != -99: 
							if motherIdx_list[motherIdx_list[iP]] == -99:		
								genFatJetIdx_top_list_spectator.append(gFJets_idx)
							elif abs(pdg_list[motherIdx_list[motherIdx_list[iP]]) == 6000055:
								genFatJetIdx_top_list_resonant.append(gFJets_idx)
																			
					elif abs(particle) == 24:
						gJets_idx = self.genP_genJet_matching(genP_idx, variable_list, fatGenJets, 1.6, 0.5)
						if gJets_index != -99:	
							if motherIdx_list[motherIdx_list[iP]] == -99:
								genFatJetIdx_w_list_spectator.append(gJets_idx)		
							elif abs(pdg_list[motherIdx_list[motherIdx_list[iP]]) == 6000055:
								genFatJetIdx_w_list_resonant.append(gJets_idx)
							
					elif abs(particle) == 5:
						gJets_idx = self.genP_genJet_matching(genP_idx, variable_list, genJets, 0.8, 0.5)
						if gJets_index != -99:	
							if motherIdx_list[motherIdx_list[iP]] == -99:	
								genJetIdx_b_list_resonant.append(gJets_idx)	
							elif abs(pdg_list[motherIdx_list[motherIdx_list[iP]]]) == 6000055:
						 		genJetIdx_b_list_spectator.append(gJets_idx)
					else: continue					
			
			
			for ak8 in fatJets_dict[eventReco_flag][0:2]:	
				if ak8.genJetAK8Idx in genFatJetIdx_list_spectator:
					genMatching_fatjet_spectator_flag.append(True)
					genMatching_fatjet_resonant_flag.append(False)				
				elif ak8.genJetAK8Idx in genFatJetIdx_list_resonant:
					genMatching_fatjet_resonant_flag.append(True)
					genMatching_fatjet_spectator_flag.append(False)
				else: 
					genMatching_fatjet_spectator_flag.append(False)
					genMatching_fatjet_resonant_flag.append(False)
			
			for ak4 in jets:
				if ak4.genJetIdx in genJetIdx_top_list_spectator:
					genMatching_jet_top_spectator_flag.append(True)
					genMatching_jet_top_resonant_flag.append(False)				
				elif ak4.genJetIdx in genJetIdx_top_list_resonant:
					genMatching_jet_top_resonant_flag.append(True)
					genMatching_jet_top_spectator_flag.append(False)
				else: 
					genMatching_jet_top_spectator_flag.append(False)
					genMatching_jet_top_resonant_flag.append(False)
				
				if ak4.genJetIdx in genJetIdx_b_list_spectator:
					genMatching_jet_b_spectator_flag.append(True)
					genMatching_jet_b_resonant_flag.append(False)				
				elif ak4.genJetIdx in genJetIdx_b_list_resonant:
					genMatching_jet_b_resonant_flag.append(True)
					genMatching_jet_b_spectator_flag.append(False)
				else: 
					genMatching_jet_b_spectator_flag.append(False)
					genMatching_jet_b_resonant_flag.append(False)	
			

			self.out.fillBranch("n"+self.outputName+"_diHadronic_FatJet", len(fatJets_dict[eventReco_flag][0:2]))
			self.out.fillBranch("n"+self.outputName+"_diHadronic_cleanedAK4s", len(jets))
			self.out.fillBranch("n"+self.outputName,len(pdg_list))
					
			for i,variable in enumerate(self.storeKinematics):
				self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) 
			
			self.out.fillBranch(self.outputName+"_pdgId", pdg_list) 
			self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) 
			self.out.fillBranch(self.outputName+"_status", status_list) 
			

			self.out.fillBranch("diHadronic_FatJet_spectatorFlag", genMatching_fatjet_spectator_flag)
			self.out.fillBranch("diHadronic_FatJet_resonantFlag", genMatching_fatjet_resonant_flag)
			
			self.out.fillBranch("diHadronic_cleanedAK4s_spectatorFlag_topMatching", genMatching_jet_top_spectator_flag)
			self.out.fillBranch("diHadronic_cleanedAK4s_resonantFlag_topMatching", genMatching_jet_top_resonant_flag)
			
			self.out.fillBranch("diHadronic_cleanedAK4s_spectatorFlag_bMatching", genMatching_jet_b_spectator_flag)
			self.out.fillBranch("diHadronic_cleanedAK4s_resonantFlag_bMatching", genMatching_jet_b_resonant_flag)
	'''
				
	setattr(event,"n"+self.outputName,len(pdg_list))
		          
        return True

