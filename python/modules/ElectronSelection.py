import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getHist, getSFXY, deltaR

class ElectronSelection(Module):
    WP80 = 1
    WP90 = 2
    INV = 3
    NONE = 4

    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Electron"),
        outputName_list = ["tightElectrons","mediumElectrons","looseElectrons"],
        triggerMatch=False,
        #electronID = WP80,
        electronMinPt = 29.,
        electronMaxEta = 2.4,
        storeKinematics= [] ,#['pt','eta'],
        storeWeights=False,
        flagMC=True
    ):

        self.inputCollection = inputCollection
        self.outputName_list = outputName_list
        self.electronMinPt = electronMinPt
        self.electronMaxEta = electronMaxEta
        self.storeKinematics = storeKinematics
        self.storeWeights = storeWeights
        self.triggerMatch = triggerMatch
        self.flagMC = flagMC
        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []

	'''
        if electronID == ElectronSelection.WP90:
            self.electronID = lambda electron: electron.mvaFall17V2Iso_WP90==1
        elif electronID == ElectronSelection.WP80:
            self.electronID = lambda electron: electron.mvaFall17V2Iso_WP80==1
        elif electronID == ElectronSelection.INV: 
            self.electronID = lambda electron: electron.mvaFall17V2Iso_WP80==0
        elif electronID == ElectronSelection.NONE:
            self.storeWeights = False
            self.electronID = lambda electron: True
        else:
            raise Exception("Electron ID undefined")
	'''

        #TODO: save the reco/ID weights if storeWeights==True
        

    def triggerMatched(self, electron, trigger_object):
        if self.triggerMatch:
            trig_deltaR = math.pi
            for trig_obj in trigger_object:
                if abs(trig_obj.id) != 11:
                    continue
                trig_deltaR = min(trig_deltaR, deltaR(trig_obj, muon))
            if trig_deltaR < 0.3:
                return True
            else:
                return False
        else:
            return True  

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        for outputName in self.outputName_list:
		self.out.branch("n"+outputName, "I")
		if not Module.globalOptions["isData"] and self.storeWeights:
		    self.out.branch(outputName+"_weight_reco_nominal","F")
		    self.out.branch(outputName+"_weight_reco_up","F")
		    self.out.branch(outputName+"_weight_reco_down","F")

		    self.out.branch(outputName+"_weight_id_nominal","F")
		    self.out.branch(outputName+"_weight_id_up","F")
		    self.out.branch(outputName+"_weight_id_down","F")

		for variable in self.storeKinematics:
		    if variable=='genPartFlav':
		    	if not Module.globalOptions["isData"]:
		    		self.out.branch(outputName+"_"+variable,"F",lenVar="n"+outputName)
			else: continue
		    else: self.out.branch(outputName+"_"+variable,"F",lenVar="n"+outputName)
			
			
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = self.inputCollection(event)
        muons = Collection(event, "Muon")
        triggerObjects = self.triggerObjectCollection(event)

        selectedElectrons = {"WP80": [], "WP90": [] , "WPL": [] }
        unselectedElectrons = []
        
        weight_reco_nominal = 1.
        weight_reco_up = 1.
        weight_reco_down = 1.

        weight_id_nominal = 1.
        weight_id_up = 1.
        weight_id_down = 1.

        for electron in electrons:
            # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            if electron.pt>self.electronMinPt \
            and math.fabs(electron.eta)<self.electronMaxEta:
            #and self.electronID(electron)\
            #and self.triggerMatched(electron, triggerObjects):
	

                dxy = math.fabs(electron.dxy)
                dz = math.fabs(electron.dz)
                
                if math.fabs(electron.eta) < 1.479 and (dxy>0.05 or dz>0.10):
                    unselectedElectrons.append(electron)
                    continue
                elif dxy>0.10 or dz>0.20:
                    unselectedElectrons.append(electron)
                    continue

                #reject electron if close-by muon
                if len(muons)>0:
                    mindr = min(map(lambda muon: deltaR(muon, electron), muons))
                    if mindr < 0.05:
                        unselectedElectrons.append(electron)
                        continue

		if electron.mvaFall17V2Iso_WP80==1: 
			selectedElectrons["WP80"].append(electron)
			selectedElectrons["WP90"].append(electron)
			selectedElectrons["WPL"].append(electron)
		elif electron.mvaFall17V2Iso_WP90==1: 
			selectedElectrons["WP90"].append(electron)
			selectedElectrons["WPL"].append(electron)
		elif electron.mvaFall17V2Iso_WPL==1: 
			selectedElectrons["WPL"].append(electron)
		
					
		'''
       	 if electronID == ElectronSelection.WP90:
            self.electronID = lambda electron: electron.mvaFall17V2Iso_WP90==1
        	elif electronID == ElectronSelection.WP80:
            self.electronID = lambda electron: electron.mvaFall17V2Iso_WP80==1
        	elif electronID == ElectronSelection.INV: 
            self.electronID = lambda electron: electron.mvaFall17V2Iso_WP80==0
        	elif electronID == ElectronSelection.NONE:
            self.storeWeights = False
            self.electronID = lambda electron: True
        	else:
            raise Exception("Electron ID undefined")
		'''
		
               #selectedElectrons.append(electron)
                
                #TODO: electron reco/ID SFs
                
                
            else:
                unselectedElectrons.append(electron)

        for outputName, electron_mvaID in zip(self.outputName_list, ['WP80','WP90','WPL']):
	
		self.out.fillBranch("n"+outputName, len(selectedElectrons[electron_mvaID]))
		
		if not Module.globalOptions["isData"] and self.storeWeights:
            
		    self.out.fillBranch(outputName+"_weight_reco_nominal", weight_reco_nominal)
		    self.out.fillBranch(outputName+"_weight_reco_up", weight_reco_up)
		    self.out.fillBranch(outputName+"_weight_reco_down", weight_reco_down)

		    self.out.fillBranch(outputName+"_weight_id_nominal",weight_id_nominal)
		    self.out.fillBranch(outputName+"_weight_id_up",weight_id_up)
		    self.out.fillBranch(outputName+"_weight_id_down",weight_id_down)


		for variable in self.storeKinematics:
		    if variable=='genPartFlav':
		    	if not Module.globalOptions["isData"]:
		    		self.out.fillBranch(outputName+"_"+variable,map(lambda electron: getattr(electron,variable), selectedElectrons[electron_mvaID]))
		    	else:
		    		continue
	    	    else:
		    	self.out.fillBranch(outputName+"_"+variable,map(lambda electron: getattr(electron,variable), selectedElectrons[electron_mvaID]))		    	
		    
		    #print list(map(lambda electron: getattr(electron,variable), selectedElectrons))

		setattr(event,outputName,selectedElectrons[electron_mvaID])
	setattr(event,"unselectedElectrons",unselectedElectrons)
		
        return True
