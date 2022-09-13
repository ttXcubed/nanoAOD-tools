import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class JetSelection(Module):

    NONE = -1
    LOOSE = 0
    TIGHT = 1
    TIGHTLEPVETO = 2
    
    def __init__(
         self,
         inputCollection=lambda event: Collection(event, "Jet"),
         leptonCollectionDRCleaning=lambda event: [],
         outputName_list=["selectedJets","unselectedJets"],
         jetMinPt=30.,
         jetMaxEta=4.8,
         dRCleaning=0.4,
         storeKinematics=['pt', 'eta','phi','mass'],
         jetId=LOOSE,
         fatFlag=True,
         metInput = lambda event: Object(event, "MET"),
     ):

        self.inputCollection = inputCollection
        self.leptonCollectionDRCleaning = leptonCollectionDRCleaning
        self.outputName_list = outputName_list
        self.jetMinPt = jetMinPt
        self.jetMaxEta = jetMaxEta
        self.dRCleaning = dRCleaning
        self.storeKinematics = storeKinematics
        self.jetId=jetId
        self.fatFlag = fatFlag
        self.metInput = metInput
        
        #loose jet ID does not exists for UL 2017 or 2018 -> accepting all jets
        if self.jetId==JetSelection.LOOSE and Module.globalOptions['year'] in ['2017','2018']:
            self.jetId = JetSelection.NONE
            

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("MET_energy", "F")
        for outputName in self.outputName_list:
		self.out.branch("n"+outputName, "I")

		if self.fatFlag:
	 	        self.out.branch(outputName+"_genJetAK8Idx", "I", lenVar="n"+outputName)
		        self.out.branch(outputName+"_deepTag_TvsQCD", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_deepTag_WvsQCD", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_particleNet_TvsQCD", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_particleNet_WvsQCD", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_particleNet_QCD", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_particleNet_mass", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_btagDeepB", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_tau2", "F", lenVar="n"+outputName)
		        self.out.branch(outputName+"_tau3", "F", lenVar="n"+outputName)
		else:
			self.out.branch(outputName+"_genJetIdx", "I", lenVar="n"+outputName)
			#self.out.branch(outputName+"_HT", "F")
			self.out.branch(outputName+"_btagDeepFlavB", "F", lenVar="n"+outputName)	
		for variable in self.storeKinematics:
		    self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)
		    
		
		

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)
	met = self.metInput(event)

        selectedJets = []
        unselectedJets = []

        leptonsForDRCleaning = self.leptonCollectionDRCleaning(event)

        for jet in jets:


            if jet.pt<self.jetMinPt:
                unselectedJets.append(jet)
                continue
        
            if math.fabs(jet.eta) > self.jetMaxEta:
                unselectedJets.append(jet)
                continue

            if self.jetId>=0 and ((jet.jetId & (1 << self.jetId)) == 0):
                unselectedJets.append(jet)
                continue

            minDeltaRSubtraction = 999.

            if len(leptonsForDRCleaning) > 0:
                mindphi = min(map(lambda lepton: math.fabs(deltaPhi(lepton, jet)), leptonsForDRCleaning))
                mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsForDRCleaning))
                
                if mindr < self.dRCleaning:
                    unselectedJets.append(jet)
                    continue
                    
                setattr(jet,"minDPhiClean",mindphi)
                setattr(jet,"minDRClean",mindr)
            else:
                setattr(jet,"minDPhiClean",100)
                setattr(jet,"minDRClean",100)
            
            selectedJets.append(jet)
            
	def metP4(obj):
                p4 = ROOT.TLorentzVector()
                p4.SetPtEtaPhiM(obj.pt,0,obj.phi,0)
                return p4
             
             
        #print((metP4(met)).E())
        self.out.fillBranch("MET_energy", (metP4(met)).E())
	
	for outputName, jet_list in zip(self.outputName_list, [selectedJets,unselectedJets]):
		self.out.fillBranch("n"+outputName, len(jet_list))
	   
		if self.fatFlag:
			self.out.fillBranch(outputName+"_genJetAK8Idx", map(lambda jet: jet.genJetAK8Idx, jet_list))
			self.out.fillBranch(outputName+"_deepTag_TvsQCD", map(lambda jet: jet.deepTag_TvsQCD, jet_list))
			self.out.fillBranch(outputName+"_deepTag_WvsQCD", map(lambda jet: jet.deepTag_WvsQCD, jet_list))
			self.out.fillBranch(outputName+"_particleNet_TvsQCD", map(lambda jet: jet.particleNet_TvsQCD, jet_list))
			self.out.fillBranch(outputName+"_particleNet_WvsQCD", map(lambda jet: jet.particleNet_WvsQCD, jet_list))
			self.out.fillBranch(outputName+"_particleNet_QCD", map(lambda jet: jet.particleNet_WvsQCD, jet_list))
			self.out.fillBranch(outputName+"_particleNet_mass", map(lambda jet: jet.particleNet_mass, jet_list))
			self.out.fillBranch(outputName+"_btagDeepB", map(lambda jet: jet.btagDeepB, jet_list)) 
			self.out.fillBranch(outputName+"_tau2", map(lambda jet: jet.tau2, jet_list)) 
			self.out.fillBranch(outputName+"_tau3", map(lambda jet: jet.tau3, jet_list)) 
			#print("jet fat idxtogen ", list(map(lambda jet: getattr(jet, 'genJetAK8Idx'), jet_list)))
		else: 
			self.out.fillBranch(outputName+"_genJetIdx", map(lambda jet: jet.genJetIdx, jet_list))
			#print("jet idxtogen ", list(map(lambda jet: getattr(jet, 'genJetIdx'), jet_list)))        	
	  		#self.out.fillBranch(outputName+"_HT", sum(list(map(lambda jet: jet.pt, jet_list))))
	  		self.out.fillBranch(outputName+"_btagDeepFlavB", map(lambda jet: jet.btagDeepFlavB, jet_list)) 
		
		for variable in self.storeKinematics:
		    self.out.fillBranch(outputName+"_"+variable, map(lambda jet: getattr(jet, variable), jet_list))

		setattr(event, outputName, jet_list)

        return True

