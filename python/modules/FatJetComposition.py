import pandas as pd

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from gen_helper import *
from utils import deltaR

class FatJetComposition(Module):
    def __init__(
            self,
            globalOptions={"isData":False}, 
            outputName=None,
            inputFatJetCollection=lambda event: Collection(event, "selectedFatJets_nominal"),
            inputGenParticleCollections={},
            outputJetPrefix='selectedHOTVRJets_nominal'
        ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.inputFatJetCollection = inputFatJetCollection
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

        self.print_out = False

    def is_inside_ak8(self, jet, ak8):
        if deltaR(jet, ak8) < 0.8:
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

        fatjets = self.inputFatJetCollection(event)
        gentops = self.inputGenParticleCollections['gentops'](event)
        genWs_not_from_top = self.inputGenParticleCollections['genWs_not_from_top'](event)
        genbs_not_from_top = self.inputGenParticleCollections['genbs_not_from_top'](event)
        genparticles_not_from_top = self.inputGenParticleCollections['genparticles_not_from_top'](event)
        
        for iak8, ak8 in enumerate(fatjets):
            if self.print_out: print('AK8 N. {}'.format(iak8+1))

            # --- initializing the flags as attributes of the jet
            for flag in self.jetCompositions:
                setattr(ak8, flag, False)
            # ---

            # --- check the origin of the jet (top, W or b not from top, pure QCD)
            # the functions to check top, W substructures are in the gen_helper.py

            if len(gentops) != 0:
                closest_gentop = min(gentops, key=lambda gentop: deltaR(gentop,ak8))
                if self.print_out: print('Closest genTop idx.{}'.format(closest_gentop._index))

                top_daughters_inside_ak8 = []
                for top_daughter in closest_gentop.daughters:
                    if deltaR(top_daughter, ak8) < 0.8:
                        top_daughters_inside_ak8.append(top_daughter)
                if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, top_daughters_inside_ak8))))

                if deltaR(closest_gentop, ak8) < 0.8:
                    if Module.globalOptions['isSignal']: 
                        if closest_gentop.from_resonance: 
                            setattr(ak8, 'has_top_fromResonance', True)
                    if self.print_out: print('Closest genTop [idx. {}] inside AK8'.format(closest_gentop._index))
                    substr_flag = gentop_substructures_check(closest_gentop, top_daughters_inside_ak8, flag_is_top_inside='topIsInside')
                else: 
                    if self.print_out: 
                        print('Daugthers {} inside AK8 of a genTop NOT inside AK8 [deltaR {} > rho/pt {}]'.format(list(map(lambda daughter: daughter.pdgId, top_daughters_inside_ak8)), deltaR(closest_gentop,ak8), 0.8))
                    
                    # --- special case of tW processes -> the ak8 comes from the W not from the top
                    if len(genWs_not_from_top) != 0:
                        closest_genW = min(genWs_not_from_top, key=lambda genW: deltaR(genW,ak8))
                        if deltaR(closest_genW, ak8) < 0.8:
                            if self.print_out: print('Closest genW not from top inside ak8')
                            W_daughters_inside_ak8 = []
                            for W_daughter in closest_genW.daughters:
                                if deltaR(W_daughter, ak8) < 0.8:
                                    W_daughters_inside_ak8.append(W_daughter)
                            if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, W_daughters_inside_ak8))))
                            substr_flag = genW_substructures_check(W_daughters_inside_ak8)
                            setattr(ak8, substr_flag, True)
                        else: 
                            setattr(ak8, 'has_other', True)
                            continue
                    # ---
                    # --- special case of ttH, H->bb -> the ak8 comes from the b from the Higgs
                    elif len(genbs_not_from_top) != 0:
                        closest_genb = min(genbs_not_from_top, key=lambda genb: deltaR(genb,ak8))
                        if deltaR(closest_genb,ak8) < 0.8:
                            if self.print_out: print('Closest b not from inside ak8')
                            setattr(ak8, 'has_b_not_fromTop', True)
                            continue
                        else: 
                            setattr(ak8, 'has_other', True)
                            continue
                    # ---
                    # --- special case: daughters are inside but also gluon radiation (which represents the biggest contribution)
                    elif any(self.is_inside_ak8(genp, ak8) for genp in list(filter(self.is_quark_lepton_gluon, genparticles_not_from_top))):
                        substr_flag = gentop_substructures_check(closest_gentop, top_daughters_inside_ak8, flag_is_top_inside='topIsNotInside_and_has_gluon_or_quark_not_fromTop')
                    # --- very rare case 
                    else: 
                        substr_flag = gentop_substructures_check(closest_gentop, top_daughters_inside_ak8, flag_is_top_inside='topIsNotInside')
                setattr(ak8, substr_flag, True)
            # ---
                
            # --- check on W, b, others not from top 
            elif len(genWs_not_from_top) != 0:
                closest_genW = min(genWs_not_from_top, key=lambda genW: deltaR(genW,ak8))
                if deltaR(closest_genW, ak8) < 0.8:
                    if self.print_out: print('Closest genW not from top inside ak8')
                    W_daughters_inside_ak8 = []
                    for W_daughter in closest_genW.daughters:
                        if deltaR(W_daughter, ak8) < 0.8:
                            W_daughters_inside_ak8.append(W_daughter)
                    if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, W_daughters_inside_ak8))))
                    substr_flag = genW_substructures_check(W_daughters_inside_ak8)
                    setattr(ak8, substr_flag, True)
                else: 
                    setattr(ak8, 'has_other', True)

            elif len(genbs_not_from_top) != 0:
                closest_genb = min(genbs_not_from_top, key=lambda genb: deltaR(genb,ak8))
                if deltaR(closest_genb, ak8) < 0.8:
                    if self.print_out: print('Closest b not from inside ak8')
                    setattr(ak8, 'has_b_not_fromTop', True)
                else: 
                    setattr(ak8, 'has_other', True)

            elif any(self.is_inside_ak8(genp, ak8) for genp in list(filter( self.is_quark_lepton_gluon, genparticles_not_from_top))):
                if self.print_out: print('Quarks/Gluons not from top inside ak8')
                setattr(ak8, 'has_gluon_or_quark_not_fromTop', True)

            else: 
                setattr(ak8, 'has_other', True)

        self.out.fillBranch("n{}".format(self.outputJetPrefix), len(fatjets))
        for composition_flag in self.jetCompositions:
            self.out.fillBranch("{}_{}".format(self.outputJetPrefix, composition_flag), map(lambda jet: getattr(jet, composition_flag), fatjets))

        return True
