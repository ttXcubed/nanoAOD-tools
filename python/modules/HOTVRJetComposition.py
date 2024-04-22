import pandas as pd

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from gen_helper import *
from utils import deltaR

class HOTVRJetComposition(Module):
    def __init__(
            self,
            globalOptions={"isData":False}, 
            outputName=None,
            inputHOTVRJetCollection=lambda event: Collection(event, "preselectedHOTVRJets"),
            inputSubHOTVRJetCollection= lambda event: Collection(event, "preselectedHOTVRSubJets"),
            inputGenParticleCollections={},
            outputJetPrefix='selectedHOTVRJets_nominal'
        ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.inputHOTVRJetCollection = inputHOTVRJetCollection
        self.inputSubHOTVRJetCollection = inputSubHOTVRJetCollection
        self.inputGenParticleCollections = inputGenParticleCollections
        self.outputJetPrefix = outputJetPrefix
        
        self.quarks_list = [1, 2, 3, 4]
        self.leptons_list = [11, 13, 15]
        self.neutrinos_list = [12, 14, 16]

        # flags to distinguish multiple processes which the jet is originated from
        # distinguishing the cases if the top quark is present or not in the jet,
        # if its decay products are inside and if contains hadronization of 
        # light quarks/gluon not coming from top

        self.jetCompositions = ['has_gluon_or_quark_not_fromTop']
        self.jetCompositions.extend(['has_other','has_b_not_fromTop'])
        for flag_top_inside in ['topIsInside','topIsNotInside','topIsNotInside_and_has_gluon_or_quark_not_fromTop']:
            self.jetCompositions.append('has_hadronicTop_'+flag_top_inside)
            self.jetCompositions.append('has_other_'+flag_top_inside)
            self.jetCompositions.append('has_noTopDaughters_'+flag_top_inside)
            if flag_top_inside == 'topIsInside' or flag_top_inside == 'topIsNotInside_and_has_gluon_or_quark_not_fromTop':
                for top_label in ['fromTop']: #,'not_fromTop']:
                    self.jetCompositions.append('has_leptonicW_'+top_label+'_'+flag_top_inside)
                    self.jetCompositions.append('has_hadronicW_'+top_label+'_'+flag_top_inside)
                    self.jetCompositions.append('has_b_plus_quark_'+top_label+'_'+flag_top_inside)
                    self.jetCompositions.append('has_b_plus_lepton_'+top_label+'_'+flag_top_inside)
                    self.jetCompositions.append('has_b_'+top_label+'_'+flag_top_inside)
                    self.jetCompositions.append('has_quark_fromW_'+top_label+'_'+flag_top_inside)
        # --- special case, samples without top and samples like tW, ttH
        self.jetCompositions.append('has_leptonicW_not_fromTop')
        self.jetCompositions.append('has_hadronicW_not_fromTop')
        self.jetCompositions.append('has_b_plus_quark_not_fromTop')
        self.jetCompositions.append('has_b_plus_lepton_not_fromTop')
        self.jetCompositions.append('has_b_not_fromTop')
        self.jetCompositions.append('has_quark_fromW_not_fromTop')

        if Module.globalOptions['isSignal']: 
            self.jetCompositions.append('has_top_fromResonance')

        self.print_out = True

    def is_inside_hotvr(self, jet, hotvr):
        rho = 600 
        effective_radius = rho / hotvr.pt if rho / hotvr.pt <= 1.5 else 1.5   # parameter defined in the paper: https://arxiv.org/abs/1606.04961
        if deltaR(jet, hotvr)< effective_radius:
            return True
        else: return False

    def is_quark_lepton_gluon(self, genP):
        if genP.pdgId in self.quarks_list or genP.pdgId in self.leptons_list or genP.pdgId == 21:
            return True
        else: return False

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.outputName is not None:
            self.out.branch(self.outputName, "I")

        self.out.branch("n{}".format(self.outputJetPrefix), "I")
        for composition_flag in self.jetCompositions:
            self.out.branch("{}_{}".format(self.outputJetPrefix, composition_flag), "I", lenVar="n{}".format(self.outputJetPrefix))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        hotvrjets = self.inputHOTVRJetCollection(event)
        subhotvrjets = self.inputSubHOTVRJetCollection(event)
        gentops = self.inputGenParticleCollections['gentops'](event)
        genWs_not_from_top = self.inputGenParticleCollections['genWs_not_from_top'](event)
        genbs_not_from_top = self.inputGenParticleCollections['genbs_not_from_top'](event)
        genparticles_not_from_top = self.inputGenParticleCollections['genparticles_not_from_top'](event)

        if self.print_out: print('\n############# HOTVR_MVA module - ev. number {}'.format(event.run))
        
        for ihotvr, hotvr in enumerate(hotvrjets):
            if self.print_out: print('HOTVR N. {}'.format(ihotvr+1))
            effective_radius = 600./ hotvr.pt if 600./ hotvr.pt <= 1.5 else 1.5

            subjets_in_hotvr = []
            for hotvr_subjet in subhotvrjets:
                if hotvr.subJetIdx1 == hotvr_subjet._index:
                    subjets_in_hotvr.insert(0, hotvr_subjet)
                if hotvr.subJetIdx2 == hotvr_subjet._index:  
                    subjets_in_hotvr.insert(1, hotvr_subjet)
                if hotvr.subJetIdx3 == hotvr_subjet._index:  
                    subjets_in_hotvr.insert(2, hotvr_subjet)                
            subjets_in_hotvr = sorted(subjets_in_hotvr, key=lambda x: x.pt, reverse=True)

            # --- initializing the flags as attributes of the jet
            for flag in self.jetCompositions:
                setattr(hotvr, flag, False)
            # ---

            # --- check the origin of the jet (top, W or b not from top, pure QCD)
            # the functions to check top, W substructures are in the gen_helper.py

            if len(gentops) != 0:
                closest_gentop = min(gentops, key=lambda gentop: deltaR(gentop,hotvr))
                if self.print_out: print('Closest genTop idx.{}'.format(closest_gentop._index))

                top_daughters_inside_hotvr = []
                for top_daughter in closest_gentop.daughters:
                    if deltaR(top_daughter, hotvr) < effective_radius:
                        top_daughters_inside_hotvr.append(top_daughter)
                if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, top_daughters_inside_hotvr))))

                if deltaR(closest_gentop, hotvr) < effective_radius:
                    if Module.globalOptions['isSignal']: 
                        if closest_gentop.from_resonance: 
                            setattr(hotvr, 'has_top_fromResonance', True)
                    if self.print_out: print('Closest genTop [idx. {}] inside HOTVR'.format(closest_gentop._index))
                    substr_flag = gentop_substructures_check(closest_gentop, top_daughters_inside_hotvr, flag_is_top_inside='topIsInside')

                else: 
                    if self.print_out: 
                        print('Daugthers {} inside HOTVR of a genTop NOT inside HOTVR [deltaR {} > rho/pt {}]'.format(list(map(lambda daughter: daughter.pdgId, top_daughters_inside_hotvr)), deltaR(closest_gentop,hotvr), effective_radius))

                    # --- special case of tW processes -> the hotvr comes from the W not from the top
                    if len(genWs_not_from_top) != 0:
                        closest_genW = min(genWs_not_from_top, key=lambda genW: deltaR(genW,hotvr))
                        if deltaR(closest_genW, hotvr) < effective_radius:
                            if self.print_out: print('Closest genW not from top inside hotvr')
                            W_daughters_inside_hotvr = []
                            for W_daughter in closest_genW.daughters:
                                if deltaR(W_daughter, hotvr) < effective_radius:
                                    W_daughters_inside_hotvr.append(W_daughter)
                            if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, W_daughters_inside_hotvr))))
                            substr_flag = genW_substructures_check(W_daughters_inside_hotvr)
                            print(substr_flag)
                            import time
                            time.sleep(1)
                            
                            setattr(hotvr, substr_flag, True)
                        else: 
                            setattr(hotvr, 'has_other', True)
                            continue
                    # ---
                    # --- special case of ttH, H->bb -> the hotvr comes from the b from the Higgs
                    elif len(genbs_not_from_top) != 0:
                        closest_genb = min(genbs_not_from_top, key=lambda genb: deltaR(genb,hotvr))
                        if deltaR(closest_genb,hotvr) < effective_radius:
                            if self.print_out: print('Closest b not from inside hotvr')
                            setattr(hotvr, 'has_b_not_fromTop', True)
                            continue
                        else: 
                            setattr(hotvr, 'has_other', True)
                            continue
                    # ---
                    # --- special case: daughters are inside but also gluon radiation (which represents the biggest contribution)
                    elif any(self.is_inside_hotvr(genp, hotvr) for genp in list(filter(self.is_quark_lepton_gluon, genparticles_not_from_top))):
                        substr_flag = gentop_substructures_check(closest_gentop, top_daughters_inside_hotvr, flag_is_top_inside='topIsNotInside_and_has_gluon_or_quark_not_fromTop')
                    # --- very rare case 
                    else: 
                        substr_flag = gentop_substructures_check(closest_gentop, top_daughters_inside_hotvr, flag_is_top_inside='topIsNotInside')  
                setattr(hotvr, substr_flag, True)
            # ---

            # --- check on W, b, others not from top 
            elif len(genWs_not_from_top) != 0:
                closest_genW = min(genWs_not_from_top, key=lambda genW: deltaR(genW,hotvr))
                if deltaR(closest_genW, hotvr) < effective_radius:
                    if self.print_out: print('Closest genW not from top inside hotvr')
                    W_daughters_inside_hotvr = []
                    for W_daughter in closest_genW.daughters:
                        if deltaR(W_daughter, hotvr) < effective_radius:
                            W_daughters_inside_hotvr.append(W_daughter)
                    if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, W_daughters_inside_hotvr))))
                    substr_flag = genW_substructures_check(W_daughters_inside_hotvr)
                    setattr(hotvr, substr_flag, True)
                else: 
                    setattr(hotvr, 'has_other', True)

            elif len(genbs_not_from_top) != 0:
                closest_genb = min(genbs_not_from_top, key=lambda genb: deltaR(genb,hotvr))
                if deltaR(closest_genb,hotvr) < effective_radius:
                    if self.print_out: print('Closest b not from inside hotvr')
                    setattr(hotvr, 'has_b_not_fromTop', True)
                else: 
                    setattr(hotvr, 'has_other', True)

            elif any(self.is_inside_hotvr(genp, hotvr) for genp in list(filter( self.is_quark_lepton_gluon, genparticles_not_from_top))):
                if self.print_out: print('Quarks/Gluons not from top inside hotvr')
                setattr(hotvr, 'has_gluon_or_quark_not_fromTop', True)

            else: 
                setattr(hotvr, 'has_other', True)


        self.out.fillBranch("n{}".format(self.outputJetPrefix), len(hotvrjets))
        for composition_flag in self.jetCompositions:
            self.out.fillBranch("{}_{}".format(self.outputJetPrefix, composition_flag), map(lambda jet: getattr(jet, composition_flag), hotvrjets))

        return True
