import os
import sys
import math
import json
import ROOT
import random
import heapq
import copy

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getGraph, getHist, combineHist2D, getSFXY, deltaR, deltaPhi
from collections import OrderedDict

class EventReconstruction(Module):

    def __init__(
        self,
        inputMuonCollection=[], #lambda event: Collection(event, "Muon"),
        inputElectronCollection=[], #lambda event: Collection(event, "Electron"),
        inputJetCollection=lambda event: Collection(event, "Jet"),
        inputFatJetCollection=lambda event: Collection(event, "Jet"),
        inputBJetCollection= [],  #ALL THE BJETS ARE AK4 JETS
        inputMETCollection=lambda event: Object(event, "MET"),
        outputName="diHadronic",
        storeKinematics_jets= ['pt','eta','phi','mass'],
        storeKinematics_leptons= ['pt','eta','phi','mass','charge','leptonFlavour'],
    ):
        
        self.inputMuonCollection = inputMuonCollection
        self.inputElectronCollection = inputElectronCollection
        self.inputJetCollection = inputJetCollection
        self.inputFatJetCollection = inputFatJetCollection
        self.inputBJetCollection = inputBJetCollection
        self.inputMETCollection = inputMETCollection
        self.outputName = outputName
        self.storeKinematics_jets = storeKinematics_jets 
 	self.storeKinematics_leptons = storeKinematics_leptons
 
    def in_Z_mass_or_in_low_resonance_window(self, lep1, lep2):
    	#print(str(lep1)[1:-4]) 
    	#print(str(lep2)[1:-4])
    	if str(lep1)[1:-4]==str(lep2)[1:-4]:

    		lep1_TLorentz, lep2_TLorentz = lep1.p4(), lep2.p4()#ROOT.TLorentzVector(), ROOT.TLorentzVector()
    		#lep1_TLorentz.SetPtEtaPhiM(lep1.pt,lep1.eta,lep1.phi,lep1.mass)
    		#lep2_TLorentz.SetPtEtaPhiM(lep2.pt,lep2.eta,lep2.phi,lep2.mass)
    		dilep_TLorentz = lep1_TLorentz+lep2_TLorentz

    		if dilep_TLorentz.M() > 80 and 100 > dilep_TLorentz.M() and dilep_TLorentz.M() > 20: return True
    		else: return False    		
    	else: return False 
    	
    def ak4_vs_ak8_cleaning(self, ak8_list, ak4_list):
    	
    	ak4_dict = {"cleaned": [], "overlapped": []}
    	for ak4 in ak4_list:
		if deltaR(ak8_list[0], ak4) >= 1.2 and deltaR(ak8_list[1], ak4) >= 1.2:
			ak4_dict["cleaned"].append(ak4)
		else: 	ak4_dict["overlapped"].append(ak4)
	
	return ak4_dict
	
    def OS_dilepton_system(self, muons_list, electrons_list): 
    
    	muons_list = sorted(muons_list,key=lambda x: x.pt, reverse=True)
	for muon in muons_list: 
		setattr(muon,"leptonFlavour",13.)	
	electrons_list = sorted(electrons_list,key=lambda x: x.pt, reverse=True)
	for electron in electrons_list:	
		setattr(electron,"leptonFlavour",11.)
		
	leptons_list = muons_list + electrons_list
	leptons_list = set(leptons_list)
	leptons_list = sorted(leptons_list,key=lambda x: x.pt, reverse=True)

	if leptons_list[0].charge != leptons_list[1].charge and self.in_Z_mass_or_in_low_resonance_window(leptons_list[0],leptons_list[1])==False and leptons_list[0].pt>=25. : ###### OPPOSITE SIGN AND NOT IN Z WINDOW
		return True, leptons_list[0:2]
	else: return False, leptons_list[0:2]
	
    def jet_objects_selection(self, ak8_list, bjets_tight_list, bjets_medium_list, bjets_loose_list, ak4_list):
    	
    	objects_selected_dict = {"cleanedBjets_tight": [], "cleanedBjets_medium": [], "cleanedBjets_loose": [], "cleanedAK4s": [], "uncleanedBjets_tight": [], "uncleanedBjets_medium": [], "uncleanedBjets_loose": [], "uncleanedAK4s": []}

	objects_selected_dict["cleanedBjets_tight"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],bjets_tight_list)["cleaned"])
	objects_selected_dict["uncleanedBjets_tight"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],bjets_tight_list)["overlapped"])
	
	objects_selected_dict["cleanedBjets_medium"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],bjets_medium_list)["cleaned"])
	objects_selected_dict["uncleanedBjets_medium"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],bjets_medium_list)["overlapped"])
	
	objects_selected_dict["cleanedBjets_loose"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],bjets_loose_list)["cleaned"])
	objects_selected_dict["uncleanedBjets_loose"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],bjets_loose_list)["overlapped"])
	
	objects_selected_dict["cleanedAK4s"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],ak4_list)["cleaned"])
	objects_selected_dict["uncleanedAK4s"].extend(self.ak4_vs_ak8_cleaning(ak8_list[0:2],ak4_list)["overlapped"])
	
	return objects_selected_dict 
	
    def leptonID_flag_selection(self, m_tight, e_tight, m_medium, e_medium, m_loose, e_loose, fj_tagger):
    	
    	leptonsID_selection_flag_dict = {}
    	OS_selected_dict = {}

    			    	
	for lepton_id in ["2tight","1tight1medium","2medium","1tight1loose","1medium1loose", "2loose"]:
	    	for tagger in ["2top", "2w", "1top1w", "low_pt_2w"]:
   			leptonsID_selection_flag_dict["OS_"+lepton_id+"_"+tagger+"_tagged"] = False
   			OS_selected_dict["OS_"+lepton_id+"_"+tagger+"_tagged"] = []	
 	
 	
    	if len(m_tight)+len(e_tight)==2:
		if self.OS_dilepton_system(m_tight,e_tight)[0]: 
			leptonsID_selection_flag_dict["OS_2tight_"+fj_tagger+"_tagged"] = True
			OS_selected_dict["OS_2tight_"+fj_tagger+"_tagged"].extend(self.OS_dilepton_system(m_tight,e_tight)[1])
	
	if len(m_medium)+len(e_medium)==2:
		if self.OS_dilepton_system(m_medium,e_medium)[0]: 
			leptonsID_selection_flag_dict["OS_2medium_"+fj_tagger+"_tagged"] = True
			OS_selected_dict["OS_2medium_"+fj_tagger+"_tagged"].extend(self.OS_dilepton_system(m_medium,e_medium)[1])
	
	if len(m_loose)+len(e_loose)==2:
		if self.OS_dilepton_system(m_loose,e_loose)[0]: 
			leptonsID_selection_flag_dict["OS_2loose_"+fj_tagger+"_tagged"] = True
			OS_selected_dict["OS_2loose_"+fj_tagger+"_tagged"].extend(self.OS_dilepton_system(m_loose,e_loose)[1])
	
	if len(m_medium)+len(e_medium)==2 and len(m_tight)+len(e_tight)>=1:
		m_selected_tightmedium = m_tight + m_medium
		e_selected_tightmedium = e_tight + e_medium

		if self.OS_dilepton_system(m_selected_tightmedium,e_selected_tightmedium)[0]: 
			leptonsID_selection_flag_dict["OS_1tight1medium_"+fj_tagger+"_tagged"] = True
			OS_selected_dict["OS_1tight1medium_"+fj_tagger+"_tagged"].extend(self.OS_dilepton_system(m_selected_tightmedium,e_selected_tightmedium)[1])
	
	if len(m_loose)+len(e_loose)==2 and len(m_tight)+len(e_tight)>=1:
		m_selected_tightloose = m_loose + m_tight
		e_selected_tightloose = e_loose + e_tight
	
		if self.OS_dilepton_system(m_selected_tightloose,e_selected_tightloose)[0]: 
			leptonsID_selection_flag_dict["OS_1tight1loose_"+fj_tagger+"_tagged"] = True
			OS_selected_dict["OS_1tight1loose_"+fj_tagger+"_tagged"].extend(self.OS_dilepton_system(m_selected_tightloose,e_selected_tightloose)[1])		

	if len(m_loose)+len(e_loose)==2 and len(m_medium)+len(e_medium)>=1:
		m_selected_mediumloose = m_loose + m_medium
		e_selected_mediumloose = e_loose + e_medium

		if self.OS_dilepton_system(m_selected_mediumloose,e_selected_mediumloose)[0]: 
			leptonsID_selection_flag_dict["OS_1medium1loose_"+fj_tagger+"_tagged"] = True
			OS_selected_dict["OS_1medium1loose_"+fj_tagger+"_tagged"].extend(self.OS_dilepton_system(m_selected_mediumloose,e_selected_mediumloose)[1])
	
 	return leptonsID_selection_flag_dict, OS_selected_dict
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
       	 
        
        for event_flag in ["OS_2tight_2top_tagged", "OS_1tight1medium_2top_tagged", "OS_2medium_2top_tagged" , "OS_1medium1loose_2top_tagged", "OS_1tight1loose_2top_tagged", "OS_2loose_2top_tagged", "OS_2tight_1top1w_tagged", "OS_1tight1medium_1top1w_tagged", "OS_2medium_1top1w_tagged" , "OS_1medium1loose_1top1w_tagged", "OS_2loose_1top1w_tagged", "OS_1tight1loose_1top1w_tagged", "OS_2tight_2w_tagged", "OS_1tight1medium_2w_tagged", "OS_2medium_2w_tagged" , "OS_1medium1loose_2w_tagged", "OS_1tight1loose_2w_tagged", "OS_2loose_2w_tagged", "OS_2tight_low_pt_2w_tagged", "OS_1tight1medium_low_pt_2w_tagged", "OS_2medium_low_pt_2w_tagged" , "OS_1medium1loose_low_pt_2w_tagged", "OS_2loose_low_pt_2w_tagged", "OS_1tight1loose_low_pt_2w_tagged"]:
        	self.out.branch(self.outputName+"_event_selection_"+event_flag, "I")	
        
        self.out.branch("n"+self.outputName+"_diFatJets", "I")
        self.out.branch("n"+self.outputName+"_diLepton", "I")
        
	for ak4_type in ["cleanedBjets_tight", "cleanedBjets_medium", "cleanedBjets_loose", "cleanedAK4s", "uncleanedBjets_tight", "uncleanedBjets_medium", "uncleanedBjets_loose", "uncleanedAK4s"]:
		self.out.branch("n"+self.outputName+"_"+ak4_type, "I")
		self.out.branch(self.outputName+"_"+ak4_type+"_btagDeepFlavB","F", lenVar="n"+self.outputName+"_"+ak4_type)
		for variable in self.storeKinematics_jets:
			self.out.branch(self.outputName+"_"+ak4_type+"_"+variable,"F", lenVar="n"+self.outputName+"_"+ak4_type)	    
	
	self.out.branch(self.outputName+"_diFatJet_deltaPhi", "F")
	
	self.out.branch(self.outputName+"_FatJet_particleNet_mass", "F",lenVar="n"+self.outputName+"_diFatJets")
	self.out.branch(self.outputName+"_FatJet_particleNet_TvsQCD", "F",lenVar="n"+self.outputName+"_diFatJets")
	self.out.branch(self.outputName+"_FatJet_particleNet_WvsQCD", "F",lenVar="n"+self.outputName+"_diFatJets")
	self.out.branch(self.outputName+"_FatJet_deepTag_WvsQCD", "F",lenVar="n"+self.outputName+"_diFatJets")
	
	for variable in self.storeKinematics_jets:
	    self.out.branch(self.outputName+"_FatJet_"+variable,"F",lenVar="n"+self.outputName+"_diFatJets")
	    self.out.branch(self.outputName+"_diFatJet_"+variable, "F")
	    
	for variable in self.storeKinematics_leptons:   #['charge']:#,'leptonFlavour']: 
  	    self.out.branch(self.outputName+"_Lepton_"+variable, "F",lenVar="n"+self.outputName+"_diLepton")
  	
  	self.out.branch(self.outputName+"_MET_pt", "F")
  	self.out.branch(self.outputName+"_MET_phi", "F")
  	
	self.out.branch(self.outputName+"_HT", "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        muons_tight = self.inputMuonCollection(event)[0]
	electrons_tight = self.inputElectronCollection(event)[0]
	muons_medium = self.inputMuonCollection(event)[1]
	electrons_medium = self.inputElectronCollection(event)[1]
	muons_loose = self.inputMuonCollection(event)[2]
	electrons_loose = self.inputElectronCollection(event)[2]
	
	jetsAK4 = self.inputJetCollection(event)
	jetsAK8 = self.inputFatJetCollection(event)
      	bjets_tight = self.inputBJetCollection(event)[0]
      	bjets_medium = self.inputBJetCollection(event)[1]
      	bjets_loose = self.inputBJetCollection(event)[2]
      	met = self.inputMETCollection(event)

	selectedAK8s_top_tagged, selectedAK8s_w_tagged, selectedAK8s_low_pt_w_tagged, unselectedAK8s = [], [], [], []

        selectedOS = {}
	jet_collection_dict =  {"cleanedBjets_tight": [], "cleanedBjets_medium": [], "cleanedBjets_loose": [], "cleanedAK4s": [], "uncleanedBjets_tight": [], "uncleanedBjets_medium": [], "uncleanedBjets_loose": [], "uncleanedAK4s": []}
        diFatJets_TLorentz_var_dict = OrderedDict()
        event_selection_flag_dict = {}
        event_selection_diFatJets_system_dict = {}
        
	for lepton_id in ["2tight","1tight1medium","2medium","1tight1loose","1medium1loose", "2loose"]:
    		for tagger in ["2top", "2w", "1top1w", "low_pt_2w"]:
			event_selection_flag_dict["OS_"+lepton_id+"_"+tagger+"_tagged"] = False
			event_selection_diFatJets_system_dict["OS_"+lepton_id+"_"+tagger+"_tagged"] = []
			selectedOS["OS_"+lepton_id+"_"+tagger+"_tagged"] = []	
		
 	
        n_diFatJet, n_diLepton = 0,0

        if len(jetsAK8)>=2:

		#print("EVENT")

        	for ak8 in jetsAK8:
        		if ak8.pt>=400:
				if ak8.particleNet_TvsQCD>=0.58: 
        				selectedAK8s_top_tagged.append(ak8)
        				
        			elif ak8.particleNet_WvsQCD>=0.71: 
        				selectedAK8s_w_tagged.append(ak8)
        			else:
        				unselectedAK8s.append(ak8)
        		else: 
        			if ak8.particleNet_WvsQCD>=0.71: 
        				selectedAK8s_low_pt_w_tagged.append(ak8)
        			else:
        				unselectedAK8s.append(ak8)
        				
		if len(selectedAK8s_top_tagged) >= 2:

			#print("2 top: ", len( bjets_medium))

			ak8_TLorentz_lead, ak8_TLorentz_sublead = selectedAK8s_top_tagged[0].p4(), selectedAK8s_top_tagged[1].p4()       
			diFatJets_TLorentz = ak8_TLorentz_lead + ak8_TLorentz_sublead

			event_selection_flag_dict = copy.deepcopy(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2top")[0])
			for flag in self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2top")[0].keys():
								
				if self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2top")[0][flag]: 
					selectedOS[flag].extend(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2top")[1][flag])

					
					for jet_collection_key in self.jet_objects_selection(selectedAK8s_top_tagged, bjets_tight, bjets_medium, bjets_loose, jetsAK4).keys():
						jet_collection_dict[jet_collection_key] = self.jet_objects_selection(selectedAK8s_top_tagged, bjets_tight, bjets_medium, bjets_loose, jetsAK4)[jet_collection_key]
					
					event_selection_diFatJets_system_dict[flag].extend([selectedAK8s_top_tagged[0],selectedAK8s_top_tagged[1], diFatJets_TLorentz])	
				else: 
					continue
					
				
		elif len(selectedAK8s_top_tagged)==1 and len(selectedAK8s_w_tagged) >= 1:
				
			selectedAK8s = copy.deepcopy(selectedAK8s_top_tagged)
			for w_tagged in selectedAK8s_w_tagged:
			#	if deltaR(selectedAK8s[0], w_tagged) >= 0.8:
				selectedAK8s.append(w_tagged)
			#	else: continue
			selectedAK8s = set(selectedAK8s)
			selectedAK8s = sorted(selectedAK8s,key=lambda x: x.pt, reverse=True)
			
			if len(selectedAK8s)<2: return
			
			ak8_TLorentz_lead, ak8_TLorentz_sublead = selectedAK8s[0].p4(), selectedAK8s[1].p4()       
			diFatJets_TLorentz = ak8_TLorentz_lead + ak8_TLorentz_sublead

			event_selection_flag_dict =  copy.deepcopy(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"1top1w")[0])
			for flag in self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"1top1w")[0].keys():


				if self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"1top1w")[0][flag]: 
					
					selectedOS[flag].extend(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"1top1w")[1][flag])	
					for jet_collection_key in self.jet_objects_selection(selectedAK8s, bjets_tight, bjets_medium, bjets_loose, jetsAK4).keys():
						jet_collection_dict[jet_collection_key] = self.jet_objects_selection(selectedAK8s, bjets_tight, bjets_medium, bjets_loose, jetsAK4)[jet_collection_key]
					
					event_selection_diFatJets_system_dict[flag].extend([selectedAK8s[0],selectedAK8s[1], diFatJets_TLorentz])	
					#print("lep sys ", selectedOS[flag])
					#print("jet col ", jet_collection_dict)
				else: 
					continue
			
		elif len(selectedAK8s_w_tagged) >= 2:
		
			print("2 w: ", len( bjets_medium))
			
			ak8_TLorentz_lead, ak8_TLorentz_sublead = selectedAK8s_w_tagged[0].p4(), selectedAK8s_w_tagged[1].p4()       
			diFatJets_TLorentz = ak8_TLorentz_lead + ak8_TLorentz_sublead
		
			event_selection_flag_dict = copy.deepcopy(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2w")[0])		
			for flag in self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2w")[0].keys():

			
				if self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2w")[0][flag]: 
					
					selectedOS[flag].extend(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"2w")[1][flag])	
					for jet_collection_key in self.jet_objects_selection(selectedAK8s_w_tagged, bjets_tight, bjets_medium, bjets_loose, jetsAK4).keys():
						jet_collection_dict[jet_collection_key] = self.jet_objects_selection(selectedAK8s_w_tagged, bjets_tight, bjets_medium, bjets_loose, jetsAK4)[jet_collection_key]
					
					event_selection_diFatJets_system_dict[flag].extend([selectedAK8s_w_tagged[0],selectedAK8s_w_tagged[1], diFatJets_TLorentz])	
					#print("lep sys ", selectedOS[flag])
					#print("jet col ", jet_collection_dict)
				else: 
					continue
			
			
		elif len(selectedAK8s_low_pt_w_tagged) >= 2:
		
		
			ak8_TLorentz_lead, ak8_TLorentz_sublead = selectedAK8s_low_pt_w_tagged[0].p4(), selectedAK8s_low_pt_w_tagged[1].p4()       
			diFatJets_TLorentz = ak8_TLorentz_lead + ak8_TLorentz_sublead
		
			event_selection_flag_dict = copy.deepcopy(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"low_pt_2w")[0])
			for flag in self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"low_pt_2w")[0].keys():
				
			
				if self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"low_pt_2w")[0][flag]: 
					
					selectedOS[flag].extend(self.leptonID_flag_selection(muons_tight, electrons_tight, muons_medium, electrons_medium, muons_loose, electrons_loose,"low_pt_2w")[1][flag])	
					for jet_collection_key in self.jet_objects_selection(selectedAK8s_low_pt_w_tagged, bjets_tight, bjets_medium, bjets_loose, jetsAK4).keys():
						jet_collection_dict[jet_collection_key] = self.jet_objects_selection(selectedAK8s_low_pt_w_tagged, bjets_tight, bjets_medium, bjets_loose, jetsAK4)[jet_collection_key]
					
					event_selection_diFatJets_system_dict[flag].extend([selectedAK8s_low_pt_w_tagged[0],selectedAK8s_low_pt_w_tagged[1], diFatJets_TLorentz])	
					#print("lep sys ", selectedOS[flag])
					#print("jet col ", jet_collection_dict)
				else: 
					continue

	#print("capa ", event_selection_flag_dict)
        for event_flag in ["OS_2tight_2top_tagged", "OS_1tight1medium_2top_tagged", "OS_2medium_2top_tagged" , "OS_1medium1loose_2top_tagged", "OS_1tight1loose_2top_tagged", "OS_2loose_2top_tagged", "OS_2tight_1top1w_tagged", "OS_1tight1medium_1top1w_tagged", "OS_2medium_1top1w_tagged" , "OS_1medium1loose_1top1w_tagged", "OS_2loose_1top1w_tagged", "OS_1tight1loose_1top1w_tagged", "OS_2tight_2w_tagged", "OS_1tight1medium_2w_tagged", "OS_2medium_2w_tagged" , "OS_1medium1loose_2w_tagged", "OS_1tight1loose_2w_tagged", "OS_2loose_2w_tagged", "OS_2tight_low_pt_2w_tagged", "OS_1tight1medium_low_pt_2w_tagged", "OS_2medium_low_pt_2w_tagged" , "OS_1medium1loose_low_pt_2w_tagged", "OS_2loose_low_pt_2w_tagged", "OS_1tight1loose_low_pt_2w_tagged"]:	
        	
		self.out.fillBranch(self.outputName+"_event_selection_"+event_flag, event_selection_flag_dict[event_flag])	

		if event_selection_flag_dict[event_flag]:

			diFatJets_TLorentz_var_dict['pt'] =  event_selection_diFatJets_system_dict[event_flag][2].Pt()
			diFatJets_TLorentz_var_dict['eta'] = event_selection_diFatJets_system_dict[event_flag][2].Eta()
			diFatJets_TLorentz_var_dict['phi'] = event_selection_diFatJets_system_dict[event_flag][2].Phi()
			diFatJets_TLorentz_var_dict['mass'] = event_selection_diFatJets_system_dict[event_flag][2].M()
			
			self.out.fillBranch("n"+self.outputName+"_diFatJets", len(event_selection_diFatJets_system_dict[event_flag][0:2]))		
			
			self.out.fillBranch(self.outputName+"_diFatJet_deltaPhi", deltaPhi(event_selection_diFatJets_system_dict[event_flag][0], event_selection_diFatJets_system_dict[event_flag][1]))

			for variable in ['particleNet_mass','particleNet_TvsQCD', 'particleNet_WvsQCD', 'deepTag_WvsQCD']:
				self.out.fillBranch(self.outputName+"_FatJet_"+variable, map(lambda ak8: getattr(ak8,variable), event_selection_diFatJets_system_dict[event_flag][0:2]))

			for variable in self.storeKinematics_jets:
				self.out.fillBranch(self.outputName+"_FatJet_"+variable, map(lambda ak8: getattr(ak8,variable), event_selection_diFatJets_system_dict[event_flag][0:2]))
				self.out.fillBranch(self.outputName+"_diFatJet_"+variable, diFatJets_TLorentz_var_dict[variable])
			
			#print(event_flag)
			for jet_collection_key in jet_collection_dict.keys():
				
				self.out.fillBranch("n"+self.outputName+"_"+jet_collection_key, len(jet_collection_dict[jet_collection_key]))
		
				for variable in self.storeKinematics_jets:
					self.out.fillBranch(self.outputName+"_"+jet_collection_key+"_"+variable, map(lambda ak4: getattr(ak4,variable), jet_collection_dict[jet_collection_key]))

				self.out.fillBranch(self.outputName+"_"+jet_collection_key+"_btagDeepFlavB", map(lambda ak4: getattr(ak4,'btagDeepFlavB'), jet_collection_dict[jet_collection_key]))
			
			
			self.out.fillBranch("n"+self.outputName+"_diLepton", len(selectedOS[event_flag]))
			
			for variable in self.storeKinematics_leptons:
		    		self.out.fillBranch(self.outputName+"_Lepton_"+variable, map(lambda lep: getattr(lep,variable), selectedOS[event_flag]))	
				
			self.out.fillBranch(self.outputName+"_MET_pt", getattr(met,'pt'))
			self.out.fillBranch(self.outputName+"_MET_phi", getattr(met,'phi'))

					
		'''
		else:
			self.out.fillBranch("n"+self.outputName+"_diFatJets", 0)	
			self.out.fillBranch("n"+self.outputName+"_diLepton", 0)

			self.out.fillBranch(self.outputName+"_diFatJet_deltaPhi", 0.)

			self.out.fillBranch(self.outputName+"_FatJet_particleNet_mass", [-2., -2.])
			self.out.fillBranch(self.outputName+"_FatJet_particleNet_TvsQCD", [0., 0.])
			self.out.fillBranch(self.outputName+"_FatJet_particleNet_WvsQCD", [0., 0.])
			self.out.fillBranch(self.outputName+"_FatJet_deepTag_WvsQCD", [0., 0.])

			for variable in self.storeKinematics_jets:
				if variable == 'pt' or variable == 'mass':
					self.out.fillBranch(self.outputName+"_FatJet_"+variable, [-2., -2.])
					self.out.fillBranch(self.outputName+"_diFatJet_"+variable, -2.)

				else: 
					self.out.fillBranch(self.outputName+"_FatJet_"+variable, [0., 0.])
					self.out.fillBranch(self.outputName+"_diFatJet_"+variable, 0.)


			for jet_collection_key in jet_collection_dict.keys():
				
				self.out.fillBranch("n"+self.outputName+"_"+jet_collection_key, 0)
		
				for variable in self.storeKinematics_jets:
					if variable == 'pt' or variable == 'mass':
						self.out.fillBranch(self.outputName+"_"+jet_collection_key+"_"+variable, [-2., -2.])
					else: 
						self.out.fillBranch(self.outputName+"_"+jet_collection_key+"_"+variable, [0., 0.])

				self.out.fillBranch(self.outputName+"_"+jet_collection_key+"_btagDeepFlavB", [0., 0.])
			
					

			for variable in self.storeKinematics_leptons:
				if variable == 'pt' or variable == 'mass': 
					self.out.fillBranch(self.outputName+"_Lepton_"+variable, [-2.,-2.])
				else: 
					self.out.fillBranch(self.outputName+"_Lepton_"+variable, [0.,0.])

			self.out.fillBranch(self.outputName+"_MET_pt", -2.)
			self.out.fillBranch(self.outputName+"_MET_phi", 0.)
		'''
		

	setattr(event, self.outputName+"_selectedAK8s", event_selection_diFatJets_system_dict)
	
	setattr(event, self.outputName+"_event_selection_flags", event_selection_flag_dict)

	setattr(event,self.outputName+"_cleanedAK4s", jet_collection_dict["cleanedAK4s"])
	setattr(event,self.outputName+"_uncleanedAK4s", jet_collection_dict["uncleanedAK4s"])	
	

        return True

