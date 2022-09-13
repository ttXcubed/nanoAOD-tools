import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class BTagSelection(Module):
    #tight DeepFlav WP (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy)
    LOOSE=0
    MEDIUM=1
    TIGHT=2
    
    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Jet"),
        flagName = "isBTagged",
        outputName_list= [], #"btaggedJets",
        jetMinPt=30.,
        jetMaxEta=2.4,
        workingpoint = [], #TIGHT,
        storeKinematics=['pt', 'eta','phi','mass'],
        storeTruthKeys=[]
    ):

        self.inputCollection = inputCollection
        self.flagName = flagName
        #self.outputName = outputName
        self.outputName_list = outputName_list
        self.jetMinPt = jetMinPt
        self.jetMaxEta = jetMaxEta
        self.storeKinematics = storeKinematics
        self.storeTruthKeys = storeTruthKeys
        self.workingpoint = workingpoint
        
        #DONE - but also available in files
        global wpValues
        wpValues = {
            '2016preVFP': [0.0614, 0.3093, 0.7221], #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
            '2016': [0.0480, 0.2489, 0.6377], #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16postVFP
            '2017': [0.0532, 0.3040, 0.7476], #https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
            '2018': [0.0490, 0.2783, 0.7100] #https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        }
        
        '''
        if workingpoint==BTagSelection.TIGHT:
            self.taggerFct = lambda jet: jet.btagDeepFlavB>wpValues[Module.globalOptions['year']][2]
        elif workingpoint==BTagSelection.MEDIUM:
            self.taggerFct = lambda jet: jet.btagDeepFlavB>wpValues[Module.globalOptions['year']][1]
        elif workingpoint==BTagSelection.LOOSE:
            self.taggerFct = lambda jet: jet.btagDeepFlavB>wpValues[Module.globalOptions['year']][0]
        else:
            raise Exception("Btagging workingpoint not defined")
        '''
         
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        for outputName in self.outputName_list:
		self.out.branch("n"+outputName, "I")
		self.out.branch(outputName+"_genJetIdx", "I", lenVar="n"+outputName)
		self.out.branch(outputName+"_btagDeepFlavB", "F", lenVar="n"+outputName)
		for variable in self.storeKinematics:
		    self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)
		for variable in self.storeTruthKeys:
		    self.out.branch(outputName+"_"+variable, "F", lenVar="n"+outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)


        #bJets_tight = []
        #bJets_medium = []
        bJets = {'tight': [], 'medium': [], 'loose': []}
        lJets = []


        for jet in jets:
            if jet.pt<self.jetMinPt:
                lJets.append(jet)
                continue
        
            if math.fabs(jet.eta) > self.jetMaxEta:
                lJets.append(jet)
                continue

	    if jet.btagDeepFlavB>wpValues[Module.globalOptions['year']][2]:
	    	#bJets_tight.append(jet)
	    	bJets['tight'].append(jet)
	    
	    if jet.btagDeepFlavB>wpValues[Module.globalOptions['year']][1]:
	    	#bJets_medium.append(jet)
	    	bJets['medium'].append(jet)
	    	
	    if jet.btagDeepFlavB>wpValues[Module.globalOptions['year']][0]:
	    	#bJets_medium.append(jet)
	    	bJets['loose'].append(jet)
	    	
            #if not self.taggerFct(jet):
            #    lJets.append(jet)
            #    continue
            
            #bJets.append(jet)
            
        #for jet in bJets:
        #    setattr(jet,self.flagName,True)
        #for jet in lJets:
        #    setattr(jet,self.flagName,False)


	for outputName, bJet_type in zip(self.outputName_list, ['tight','medium','loose']):
	
		self.out.fillBranch("n"+outputName, len(bJets[bJet_type]))
		
		
		self.out.fillBranch(outputName+"_genJetIdx", map(lambda jet: getattr(jet, 'genJetIdx'), bJets[bJet_type]))
		self.out.fillBranch(outputName+"_btagDeepFlavB", map(lambda jet: getattr(jet, 'btagDeepFlavB'), bJets[bJet_type]))
		#print("b jet idxtogen ", list(map(lambda jet: getattr(jet, 'genJetIdx'), bJets)))        	
		
		for variable in self.storeKinematics:
		    self.out.fillBranch(outputName+"_"+variable,
		                        map(lambda jet: getattr(jet, variable), bJets[bJet_type]))

		for variable in self.storeTruthKeys:
		    self.out.fillBranch(outputName+"_"+variable,
		                        map(lambda jet: getattr(jet, variable), bJets[bJet_type]))


        	setattr(event, outputName, bJets[bJet_type])
        	#setattr(event, self.outputName+"_unselected", lJets)
                
        return True

