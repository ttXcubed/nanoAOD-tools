import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from functools import reduce
from utils import getHist,combineHist2D,getSFXY,deltaR,deltaPt

class DoubleLeptonTriggerSelection(Module):
    def __init__(
        self,
        dileptonTriggerSFFileName,
        inputMuonCollection = lambda event: getattr(event,"looseMuons"),
        inputElectronCollection = lambda event: getattr(event,"looseElectrons"),
        storeWeights=False,
        outputName = "trigger",
        thresholdPt=15.,  #preliminary value  
    ):
        self.dileptonTriggerSFFileName = dileptonTriggerSFFileName
        self.inputMuonCollection = inputMuonCollection
        self.inputElectronCollection = inputElectronCollection
        self.outputName = outputName
        self.storeWeights = storeWeights
        self.thresholdPt = thresholdPt
                    
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_mumu_flag","I")
        self.out.branch(self.outputName+"_ee_flag","I")
        self.out.branch(self.outputName+"_emu_flag","I")
        self.out.branch(self.outputName+"_general_flag","I")
        
        if not self.globalOptions["isData"] and self.storeWeights:
            self.out.branch(self.outputName+"_weight_nominal","F")
            self.out.branch(self.outputName+"_weight_up","F")
            self.out.branch(self.outputName+"_weight_down","F")
            
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        muons = self.inputMuonCollection(event)
        electrons = self.inputElectronCollection(event)

        mumu_flag = 0
        ee_flag = 0
        emu_flag = 0
        general_flag = 0
        
        weight_trigger_nominal = 1.
        weight_trigger_up = 1.
        weight_trigger_down = 1.

        # trigger flags recomended by TOP PAG: https://gitlab.cern.ch/tdr/notes/AN-19-140/-/blob/miniAODv2/TrigEffAN_temp.pdf, https://indico.cern.ch/event/995560/contributions/4189577/attachments/2174069/3671077/Dilepton_TriggerSF_TOPPAG.pdf
        emu_flag, mumu_flag, ee_flag = False, False, False
        if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':
            emu_flag = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
            mumu_flag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL or event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL or event.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ 
            ee_flag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        
        if Module.globalOptions["year"] == '2017':
            emu_flag = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
            mumu_flag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            ee_flag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL>0
        
        if Module.globalOptions["year"] == '2018':
            emu_flag = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
            mumu_flag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 
            ee_flag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL>0
        
        if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
            emu_flag = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
            mumu_flag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            ee_flag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL

        #SFs are provided in range of leading/subleading lepton pT therefore are taken from h2D (leading vs subleading pT)
        triggerSFHist = None
        if (not Module.globalOptions["isData"]) and self.storeWeights: 
            if emu_flag and len(muons)>0 and len(electrons)>0:
                triggerSFHist = getHist(
                        self.dileptonTriggerSFFileName,
                        'h2D_SF_emu_lepABpt_FullError'
                )
                weight_trigger,weight_trigger_err = getSFXY(triggerSFHist, electrons[0].pt, muons[0].pt)
                weight_trigger_nominal *= weight_trigger
                weight_trigger_up *= (weight_trigger+weight_trigger_err)
                weight_trigger_down *= (weight_trigger-weight_trigger_err)
            elif mumu_flag and len(muons)>1:
                triggerSFHist = getHist(
                        self.dileptonTriggerSFFileName,
                        'h2D_SF_mumu_lepABpt_FullError'
                )
                weight_trigger,weight_trigger_err = getSFXY(triggerSFHist, muons[0].pt, muons[1].pt)
                weight_trigger_nominal *= weight_trigger
                weight_trigger_up *= (weight_trigger+weight_trigger_err)
                weight_trigger_down *= (weight_trigger-weight_trigger_err)
            elif ee_flag and len(electrons)>1:
                triggerSFHist = getHist(
                        self.dileptonTriggerSFFileName,
                        'h2D_SF_ee_lepABpt_FullError'
                )
                weight_trigger,weight_trigger_err = getSFXY(triggerSFHist, electrons[0].pt, electrons[1].pt)
                weight_trigger_nominal *= weight_trigger
                weight_trigger_up *= (weight_trigger+weight_trigger_err)
                weight_trigger_down *= (weight_trigger-weight_trigger_err)
       
        general_flag = mumu_flag or emu_flag or ee_flag

        self.out.fillBranch(self.outputName+"_mumu_flag", mumu_flag)
        self.out.fillBranch(self.outputName+"_ee_flag", ee_flag)
        self.out.fillBranch(self.outputName+"_emu_flag", emu_flag)
        self.out.fillBranch(self.outputName+"_general_flag", general_flag)
            
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_nominal",weight_trigger_nominal)
            self.out.fillBranch(self.outputName+"_weight_up",weight_trigger_up)
            self.out.fillBranch(self.outputName+"_weight_down",weight_trigger_down)

        setattr(event, self.outputName+"_general_flag", general_flag)
        for os, os_flag in zip(['ee','emu','mumu'],[ee_flag, emu_flag, mumu_flag]):
            setattr(event, self.outputName+"_"+os+"_flag", os_flag)

        return True
        
