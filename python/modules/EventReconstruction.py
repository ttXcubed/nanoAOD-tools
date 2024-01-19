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

from utils import getGraph, getHist, combineHist2D, getSFXY, deltaR, deltaPhi, PhysicsObject
from collections import OrderedDict

class EventReconstruction(Module):

    def __init__(
        self,
        inputTriggersCollection = {},
        inputMuonCollection = {},
        inputElectronCollection = {},
        inputJetCollection = {},
        inputBJetCollection = {},
        inputFatJetCollection = {},
        inputMETCollection = [],
        inputHOTVRJetCollection = {},
        inputHOTVRSubJetCollection = {},
        inputGenTopCollection = {}, 
    ):
        
        self.inputTriggersCollection = inputTriggersCollection
        self.inputMuonCollection = inputMuonCollection
        self.inputElectronCollection = inputElectronCollection
        self.inputJetCollection = inputJetCollection
        self.inputFatJetCollection = inputFatJetCollection
        self.inputBJetCollection = inputBJetCollection
        self.inputMETCollection = inputMETCollection
        self.inputHOTVRJetCollection = inputHOTVRJetCollection
        self.inputHOTVRSubJetCollection = inputHOTVRSubJetCollection
        if not Module.globalOptions['isData']:
            self.inputGenTopCollection = inputGenTopCollection
    
        self.print_out = False

        self.eventReconstructionKeys_ak4Jets = [
            'is_inside_hotvr','is_inside_ak8','min_deltaRVSak8Jet',
            'min_deltaRVShotvrJet', "rho_over_pt_closest_hotvr"
            ]

        self.hotvr_vars = [
            'nsubjets', 'tau3_over_tau2', 'fractional_subjet_pt', 
            'min_pairwise_subjets_mass', 'effective_radius', 
            'tau2_over_tau1', 
            ]
        self.jet_composition = ['has_genTop_inside', 'has_genTopHadronic_inside']
        if  Module.globalOptions['isSignal']: self.jet_composition.append('has_genTopFromResonance_inside')

    def minimum_pairwise_mass(self, subjets_in_hotvr):
        if len(subjets_in_hotvr)>2:
            min_pair_sum_mass = float('inf')

            for i, sjet_i in enumerate(subjets_in_hotvr[:2]):  # minimum pairwise mass of the three leading subjets
                for j, sjet_j in enumerate(subjets_in_hotvr[i+1:3]):
                    pair_sum = sjet_i.p4() + sjet_j.p4()
                    if pair_sum.M() < min_pair_sum_mass:
                        min_pair_sum_mass = pair_sum.M()
            # if self.print_out: print("The minimum sum pair mass is :", min_pair_sum_mass)
            return min_pair_sum_mass
        else: return -99.

    def is_inside_ak8(self, jet, ak8):
        if deltaR(jet, ak8)<0.8:
            return True
        else: return False

    def is_inside_hotvr(self, jet, hotvr):
        rho = 600 
        effective_radius = 600./ hotvr.pt if 600./ hotvr.pt <= 1.5 else 1.5   # parameter defined in the paper: https://arxiv.org/abs/1606.04961
        if deltaR(jet, hotvr)< effective_radius:
            return True
        else: return False

    def is_quark_lepton_gluon(self, genP):
        if genP.pdgId in self.quarks_list or genP.pdgId in self.leptons_list or genP.pdgId==21:
            return True
        else: return False

    def is_OS_lepton(self, leptons, **kwargs):
        if len(leptons) == 2:
            if kwargs['lepton_selection'] != 'emu':
                first_lep_p4, second_lep_p4 = ROOT.TLorentzVector(), ROOT.TLorentzVector()
                first_lep_p4.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, leptons[0].mass)
                second_lep_p4.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, leptons[1].mass)
                dilepton_system = first_lep_p4 + second_lep_p4

                if self.print_out:
                    print('Leading lepton pt {}; leptons charges {},{}; dilepton system invariant mass {}'.format(leptons[0].pt, leptons[0].charge, leptons[1].charge, dilepton_system.M()))
                return (leptons[0].charge*leptons[1].charge)<1 and leptons[0].pt>25 and dilepton_system.M() > 20. #and (dilepton_system.M() < 80. or 100. < dilepton_system.M())

            else: return (leptons[0].charge*leptons[1].charge)<1 and leptons[0].pt>25
        else: return False

    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for cut_selection in ['ee','emu','mumu']:
            self.out.branch('eventSelection_'+cut_selection+'_cut',"I")

        for jet_collection in ['selectedJets_nominal', 'preselectedHOTVRJets','selectedFatJets_nominal']:
            
            self.out.branch("n"+jet_collection,"I")
            # for largeRadiusJet in ['ak8','hotvr']:
            #     self.out.branch("n"+jet_collection+'_inside_'+largeRadiusJet,"I")
            #     self.out.branch("n"+jet_collection+'_outside_'+largeRadiusJet,"I")
            if jet_collection == 'selectedJets_nominal':
                for var in self.eventReconstructionKeys_ak4Jets:
                    if 'min' in var or 'rho' in var :
                        self.out.branch(jet_collection+"_"+var, "F", lenVar="n"+jet_collection)
                    else: self.out.branch(jet_collection+"_"+var, "I", lenVar="n"+jet_collection)

            if jet_collection == 'preselectedHOTVRJets':
                for var in self.hotvr_vars:
                    if var != 'nsubjets':
                        self.out.branch(jet_collection+"_"+var, "F", lenVar="n"+jet_collection)
                    else: self.out.branch(jet_collection+"_"+var, "I", lenVar="n"+jet_collection)

            if jet_collection != 'selectedJets_nominal':
                for jet_composition_flag in self.jet_composition:
                    self.out.branch(jet_collection+"_"+jet_composition_flag, "I", lenVar="n"+jet_collection)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        hotvrjets = self.inputHOTVRJetCollection(event)
        hotvrsubjets = self.inputHOTVRSubJetCollection(event)
        jets = self.inputJetCollection(event)
        bjets = self.inputBJetCollection(event)
        fatjets = self.inputFatJetCollection(event)
        muons = self.inputMuonCollection(event)
        electrons = self.inputElectronCollection(event)
        if not Module.globalOptions['isData']:
            gentops = self.inputGenTopCollection(event)

        # separation needed for data/MC:
        # data -> separate indipendent modules for ee, emu, mumu triggers
        # MC -> cumulative dilepton trigger (ee OR emu OR mumu)
        if Module.globalOptions['isData']:
            triggers = {'ee': False, 'emu': False, 'mumu': False}
            triggers[Module.globalOptions['trigger']] = self.inputTriggersCollection[Module.globalOptions['trigger']](event)
        else: triggers = {'ee': self.inputTriggersCollection['ee'](event), 'emu': self.inputTriggersCollection['emu'](event), 'mumu': self.inputTriggersCollection['mumu'](event)}

        # print('############# EVENT')
        # ---- AK4 jet cleaning 
        if self.print_out: 
            # print("####################### EVENT")
            print("---- Event Reconstruction Module ----")
            print('# of hotvr jets: {}, ak8 jets: {}, ak4: {}'.format(len(hotvrjets), len(fatjets), len(jets)))
        
        for ak4 in jets:
            setattr(ak4, "min_deltaRVSak8Jet", -99), setattr(ak4, "is_inside_ak8", False), setattr(ak4, "is_outside_ak8", False)
            if len(fatjets)>0:
                closest_ak8 = min(fatjets, key=lambda ak8: deltaR(ak8,ak4))
                inside_any_ak8 = any( self.is_inside_ak8(ak4, ak8) for ak8 in fatjets)
                setattr(ak4, "min_deltaRVSak8Jet", deltaR(closest_ak8,ak4))
                setattr(ak4, "is_inside_ak8", inside_any_ak8)

            setattr(ak4, "min_deltaRVShotvrJet", -99), setattr(ak4, "rho_over_pt_closest_hotvr", -99), setattr(ak4, "is_inside_hotvr", False), setattr(ak4, "is_outside_hotvr", False)
            if len(hotvrjets) > 0:
                closest_hotvr = min(hotvrjets, key=lambda hotvr: deltaR(hotvr,ak4))
                inside_any_hotvr = any( self.is_inside_hotvr(ak4, hotvr) for hotvr in hotvrjets)
                setattr(ak4, "min_deltaRVShotvrJet", deltaR(closest_hotvr,ak4))
                setattr(ak4, "rho_over_pt_closest_hotvr", 600./closest_hotvr.pt)
                setattr(ak4, "is_inside_hotvr", inside_any_hotvr)
        # ---- 

        # ---- OS CUT
        event_selection_OS_dilepton_cut = False
        event_selection_dilepton = {'ee': False, 'emu': False, 'mumu': False}
        if triggers['ee'] and len(muons) == 0:
            if self.print_out: print('Trigger {}'.format('ee'))
            event_selection_dilepton['ee'] = self.is_OS_lepton(electrons, lepton_selection='ee')
        if triggers['emu'] and len(electrons) == 1 and len(muons) == 1:
            if self.print_out: print('Trigger {}'.format('emu'))
            leptons = muons + electrons
            leptons.sort(key=lambda lep: lep.pt, reverse=True)
            event_selection_dilepton['emu'] = self.is_OS_lepton(leptons, lepton_selection='emu')
        if triggers['mumu'] and len(electrons) == 0:
            if self.print_out: print('Trigger {}'.format('mumu'))
            event_selection_dilepton['mumu'] = self.is_OS_lepton(muons, lepton_selection='mumu')

        event_selection_OS_dilepton_cut = event_selection_dilepton['ee'] or event_selection_dilepton['emu'] or event_selection_dilepton['mumu']
        # ----

        # ---- AK4, b CUT
        # // case 0 hotvr: check ak8 -> 1) no ak8; 2) 1 ak8 with reso flag; 3) 2 ak8 with reso flag
        # // case 1 hotvr with reso flag: check ak8 -> 1) no ak8; 2) 1 ak8 with reso flag (and no overlap); 
        # // case 2 hotvr with reso flag
        # ----

        # ---- HOTVR VARIABLES
        for ihotvr, hotvr in enumerate(hotvrjets):
            if self.print_out: print('HOTVR N. {}'.format(ihotvr+1))

            effective_radius = 600./ hotvr.pt if 600./ hotvr.pt <= 1.5 else 1.5
            setattr(hotvr, 'effective_radius', effective_radius)

            subjets_in_hotvr = []
            for hotvr_subjet in hotvrsubjets:
                if hotvr.subJetIdx1 == hotvr_subjet._index:
                    subjets_in_hotvr.insert(0, hotvr_subjet)
                if hotvr.subJetIdx2 == hotvr_subjet._index:  
                    subjets_in_hotvr.insert(1, hotvr_subjet)
                if hotvr.subJetIdx3 == hotvr_subjet._index:  
                    subjets_in_hotvr.insert(2, hotvr_subjet)                
            subjets_in_hotvr = sorted(subjets_in_hotvr,key=lambda x: x.pt, reverse=True)

            setattr(hotvr, 'nsubjets', len(subjets_in_hotvr))
            if hotvr.tau2 != 0.: 
                setattr(hotvr, 'tau3_over_tau2', hotvr.tau3/hotvr.tau2 )
            else: setattr(hotvr, 'tau3_over_tau2', -1 )
            if hotvr.tau1 != 0.:
                setattr(hotvr, 'tau2_over_tau1', hotvr.tau2/hotvr.tau1 )
            else: setattr(hotvr, 'tau2_over_tau1', -1 )
            if len(subjets_in_hotvr)!=0:
                setattr(hotvr, 'fractional_subjet_pt', subjets_in_hotvr[0].pt/hotvr.pt )
            else: 
                setattr(hotvr, 'fractional_subjet_pt', 0.) #-99 )
            setattr(hotvr, 'min_pairwise_subjets_mass', self.minimum_pairwise_mass(subjets_in_hotvr) )

            # gen-composition in top samples
            if not Module.globalOptions['isData']:
                setattr(hotvr, 'has_genTop_inside', False), setattr(hotvr, 'has_genTopFromResonance_inside', False), setattr(hotvr, 'has_genTopHadronic_inside', False)
                if len(gentops) != 0:
                    closest_gentop = min(gentops, key=lambda gentop: deltaR(gentop, hotvr))
                    if self.is_inside_hotvr(closest_gentop, hotvr):
                        setattr(hotvr, 'has_genTop_inside', True)
                        if closest_gentop.has_hadronically_decay: setattr(hotvr, 'has_genTopHadronic_inside', True)
                        if Module.globalOptions['isSignal'] and closest_gentop.from_resonance: setattr(hotvr, 'has_genTopFromResonance_inside', True)
        # ----
                
        # ---- AK8 VARIABLES
        for iak8, ak8 in enumerate(fatjets):
            if self.print_out: print('HOTVR N. {}'.format(iak8+1))

            # gen-composition in top samples
            if not Module.globalOptions['isData']:
                setattr(ak8, 'has_genTop_inside', False), setattr(ak8, 'has_genTopFromResonance_inside', False), setattr(ak8, 'has_genTopHadronic_inside', False)
                if len(gentops) != 0:
                    closest_gentop = min(gentops, key=lambda gentop: deltaR(gentop, ak8))
                    if self.is_inside_ak8(closest_gentop, ak8):
                        setattr(ak8, 'has_genTop_inside', True)
                        if closest_gentop.has_hadronically_decay: setattr(ak8, 'has_genTopHadronic_inside', True)
                        if Module.globalOptions['isSignal'] and closest_gentop.from_resonance: setattr(ak8, 'has_genTopFromResonance_inside', True)


        setattr(event, 'preselectedHOTVRJets', hotvrjets)
        setattr(event, 'selectedFatJets_nominal', fatjets)

        self.out.fillBranch("nselectedJets_nominal", len(jets))
        self.out.fillBranch("npreselectedHOTVRJets", len(hotvrjets))
        self.out.fillBranch("nselectedFatJets_nominal", len(fatjets))


        for var in self.eventReconstructionKeys_ak4Jets:
            self.out.fillBranch("selectedJets_nominal_"+var, map(lambda jet: getattr(jet,var), jets))
        for var in self.hotvr_vars:
            self.out.fillBranch("preselectedHOTVRJets_"+var, map(lambda hotvr: getattr(hotvr,var), hotvrjets))

        if not Module.globalOptions['isData']:
            for jet_composition_flag in self.jet_composition:
                self.out.fillBranch("preselectedHOTVRJets_"+jet_composition_flag, map(lambda hotvr: getattr(hotvr, jet_composition_flag), hotvrjets))
                self.out.fillBranch("selectedFatJets_nominal_"+jet_composition_flag, map(lambda ak8: getattr(ak8, jet_composition_flag), fatjets))

        for cut_selection in ['ee','emu','mumu']:
            self.out.fillBranch('eventSelection_'+cut_selection+'_cut', event_selection_dilepton[cut_selection])
        setattr(event, 'event_selection_OS_dilepton_cut', event_selection_OS_dilepton_cut)

        return True

