import pandas as pd

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from gen_helper import *
from utils import deltaR

class HOTVR_MVA(Module):
    def __init__(self,globalOptions={"isData":False}, outputName=None,
            inputHOTVRJetCollection=lambda event: Collection(event, "preselectedHOTVRJets"),
            inputSubHOTVRJetCollection= lambda event: Collection(event, "preselectedHOTVRSubJets"),
            inputGenParticleCollections={},
        ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.inputHOTVRJetCollection = inputHOTVRJetCollection
        self.inputSubHOTVRJetCollection = inputSubHOTVRJetCollection
        self.inputGenParticleCollections = inputGenParticleCollections
        
        self.hotvrjets_df_list = []

        self.hotvr_vars = [
            'pt', 'eta', 'phi', 'mass', 'nsubjets', 'tau3_over_tau2', 
            'fractional_subjet_pt', 'min_pairwise_subjets_mass',
            'tau2_over_tau1'
            ]

        self.quarks_list = [1, 2, 3, 4]
        self.leptons_list = [11, 13, 15]
        self.neutrinos_list = [12, 14, 16]

        self.MVA_keys = ['has_gluon_or_quark_not_fromTop']
        self.MVA_keys.extend(['has_other','has_b_not_fromTop'])
        for flag_top_inside in ['topIsInside','topIsNotInside','topIsNotInside_and_has_gluon_or_quark_not_fromTop']:
            self.MVA_keys.append('has_hadronicTop_'+flag_top_inside)
            self.MVA_keys.append('has_other_'+flag_top_inside)
            self.MVA_keys.append('has_noTopDaughters_'+flag_top_inside)
            for top_label in ['fromTop','not_fromTop']:
                self.MVA_keys.append('has_leptonicW_'+top_label+'_'+flag_top_inside)
                self.MVA_keys.append('has_hadronicW_'+top_label+'_'+flag_top_inside)
                self.MVA_keys.append('has_b_plus_quark_'+top_label+'_'+flag_top_inside)
                self.MVA_keys.append('has_b_plus_lepton_'+top_label+'_'+flag_top_inside)
                self.MVA_keys.append('has_b_'+top_label+'_'+flag_top_inside)
                self.MVA_keys.append('has_quark_fromW_'+top_label+'_'+flag_top_inside)
        for var in self.hotvr_vars:
            self.MVA_keys.append(var)
        self.MVA_keys.append('run_number')

        if Module.globalOptions['isSignal']: 
            self.MVA_keys.append('has_top_fromResonance')

        self.print_out = False

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

    def gen_top_substructures_check(self, hotvr_dict, gen_top_inside_jet, top_daughters_inside_jet, **keys):
        if len(top_daughters_inside_jet) == 3:
            if gen_top_inside_jet.has_hadronically_decay:
                if self.print_out: print('Has hadronic decay ...daugthers are {} {} {}'.format(gen_top_inside_jet.first_daughter, gen_top_inside_jet.second_daughter, gen_top_inside_jet.third_daughter))
                hotvr_dict['has_hadronicTop_'+keys['flag_is_top_inside']] = True
            elif not gen_top_inside_jet.has_hadronically_decay:
                if self.print_out: print('Has leptonic decay ...daugthers are {} {} {}'.format(gen_top_inside_jet.first_daughter, gen_top_inside_jet.second_daughter, gen_top_inside_jet.third_daughter))
                hotvr_dict['has_b_plus_lepton_fromTop_'+keys['flag_is_top_inside']] = True # scenario with neutrino inside is not considered as a separate category with respect to the one with just 2 objects inside the jet (b+lep)
            else: 
                if self.print_out: print('Has wierd decay ...daugthers are {} {} {}'.format(gen_top_inside_jet.first_daughter, gen_top_inside_jet.second_daughter, gen_top_inside_jet.third_daughter))
                hotvr_dict['has_other_'+keys['flag_is_top_inside']] = True

        elif len(top_daughters_inside_jet) == 2:
            daughters_pdg = list(map(lambda top_daughter: abs(top_daughter.pdgId), top_daughters_inside_jet))
            if gen_top_inside_jet.has_hadronically_decay and not any(daughter == 5 for daughter in daughters_pdg):
                if self.print_out: print('Has hadronic decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_hadronicW_fromTop_'+keys['flag_is_top_inside']] = True
            elif not gen_top_inside_jet.has_hadronically_decay and not any(daughter == 5 for daughter in daughters_pdg): 
                if self.print_out: print('Has leptonic decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_leptonicW_fromTop_'+keys['flag_is_top_inside']] = True
            elif any(daughter == 5 for daughter in daughters_pdg) and any(daughter in self.quarks_list for daughter in daughters_pdg): 
                if self.print_out: print('Has b+q(inside) decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_b_plus_quark_fromTop_'+keys['flag_is_top_inside']] = True
            elif any(daughter == 5 for daughter in daughters_pdg) and any(daughter in self.leptons_list for daughter in daughters_pdg):  
                if self.print_out: print('Has b+l(inside) decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_b_plus_lepton_fromTop_'+keys['flag_is_top_inside']] = True
            elif any(daughter == 5 for daughter in daughters_pdg) and any(daughter in self.neutrinos_list for daughter in daughters_pdg): 
                if self.print_out: print('Has b+neutrino(inside) decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_b_fromTop_'+keys['flag_is_top_inside']] = True
            else: 
                if self.print_out: print('Has wierd decay ...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_other_'+keys['flag_is_top_inside']] = True

        elif len(top_daughters_inside_jet) == 1:
            if self.print_out: print('Daugther is {}'.format(top_daughters_inside_jet[0].pdgId))
            if abs(top_daughters_inside_jet[0].pdgId) == 5:
                hotvr_dict['has_b_fromTop_'+keys['flag_is_top_inside']] = True
                if self.print_out: print('...has_b_fromTop_'+keys['flag_is_top_inside'])
            elif abs(top_daughters_inside_jet[0].pdgId) in self.quarks_list:
                hotvr_dict['has_quark_fromW_fromTop_'+keys['flag_is_top_inside']] = True
                if self.print_out: print('...has_quark_fromW_fromTop_'+keys['flag_is_top_inside'])
            elif abs(top_daughters_inside_jet[0].pdgId) in self.leptons_list:
                hotvr_dict['has_leptonicW_fromTop_'+keys['flag_is_top_inside']] = True
                if self.print_out: print('...has_leptonicW_fromW_fromTop_'+keys['flag_is_top_inside'])
            else: 
                hotvr_dict['has_other_'+keys['flag_is_top_inside']] = True
                if self.print_out: print('...has_other_'+keys['flag_is_top_inside'])
        else: 
            hotvr_dict['has_noTopDaughters_'+keys['flag_is_top_inside']] = True
            if self.print_out: print('...has_noTopDaughters_'+keys['flag_is_top_inside'])

    def gen_W_substructures_check(self, hotvr_dict, W_daughters_inside_jet):
        if len(W_daughters_inside_jet)==2:
            daughters_pdg = list(map(lambda top_daughter: abs(top_daughter.pdgId), W_daughters_inside_jet))
            are_W_leptonic = all(daughter in self.leptons_list for daughter in daughters_pdg)
            are_W_hadronic = all(daughter in self.quarks_list for daughter in daughters_pdg)

            if are_W_hadronic: 
                if self.print_out: print('Has hadronic decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_hadronicW_not_fromTop'] = True
            elif are_W_leptonic: 
                if self.print_out: print('Has leptonic decay...daugthers are {}'.format(daughters_pdg))
                hotvr_dict['has_leptonicW_not_fromTop'] = True
            else: 
                hotvr_dict['has_other'] = True
        
        elif len(W_daughters_inside_jet)==1:
            if abs(W_daughters_inside_jet[0].pdgId) in self.quarks_list:
                hotvr_dict['has_quark_fromW_not_fromTop'] = True
            elif abs(W_daughters_inside_jet[0].pdgId) in self.leptons_list:
                hotvr_dict['has_leptonicW_not_fromTop'] = True
            else: 
                hotvr_dict['has_other'] = True
        else:
            hotvr_dict['has_other'] = True
    

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.outputName is not None:
            self.out.branch(self.outputName, "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if len(self.hotvrjets_df_list)!=0:
            tot_hotvrjets_df = pd.concat(self.hotvrjets_df_list, ignore_index=True, sort=False)
            tot_hotvrjets_df.to_hdf(outputFile.GetName().replace('.root','.h5'), key='hotvrjets_dataset', mode='w')

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

            # --- initializing the dictionary
            hotvr_dict = dict()
            for key in self.MVA_keys:
                hotvr_dict[key] = False
            # ---

            for var in self.hotvr_vars:
                hotvr_dict[var] = getattr(hotvr, var)
            
            # saving run number so the training will be performed on jets with odd event number
            # while the evaluation on jets with even number to avoid overtraining 
            hotvr_dict['run_number'] = event.run

            # --- check on top
            if len(gentops) != 0:
                closest_gentop = min(gentops, key=lambda gentop: deltaR(gentop,hotvr))
                if self.print_out: print('Closest genTop idx.{}'.format(closest_gentop._index))
                if Module.globalOptions['isSignal']: 
                    if closest_gentop.from_resonance: hotvr_dict['has_top_fromResonance'] = True
                
                top_daughters_inside_hotvr = []
                for top_daughter in closest_gentop.daughters:
                    if deltaR(top_daughter, hotvr) < effective_radius:
                        top_daughters_inside_hotvr.append(top_daughter)
                if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, top_daughters_inside_hotvr))))

                if deltaR(closest_gentop, hotvr) < effective_radius:
                    if self.print_out: print('Closest genTop [idx. {}] inside HOTVR'.format(closest_gentop._index))
                    self.gen_top_substructures_check(hotvr_dict, closest_gentop, top_daughters_inside_hotvr, flag_is_top_inside='topIsInside')
                else: 
                    if self.print_out: print('Daugthers {} inside HOTVR of a closed genTop NOT inside HOTVR [deltaR {} > rho/pt {}]'.format(list(map(lambda daughter: daughter.
                    pdgId, top_daughters_inside_hotvr)), deltaR(closest_gentop,hotvr), effective_radius))
                    # print(closest_gentop.pt, hotvr.pt)
                    
                    if any(self.is_inside_hotvr(genp, hotvr) for genp in list(filter( self.is_quark_lepton_gluon, genparticles_not_from_top))):
                        self.gen_top_substructures_check(hotvr_dict, closest_gentop, top_daughters_inside_hotvr, flag_is_top_inside='topIsNotInside_and_has_gluon_or_quark_not_fromTop')
                    else: 
                        self.gen_top_substructures_check(hotvr_dict, closest_gentop, top_daughters_inside_hotvr, flag_is_top_inside='topIsNotInside')   
            # ---
                
            # --- check on W, b, others not from top 
            elif len(genWs_not_from_top) != 0:
                closest_genW = min(genWs_not_from_top, key=lambda genW: deltaR(genW,hotvr))
                if deltaR(closest_genW,hotvr) < effective_radius:
                    if self.print_out: print('Closest genW not from top inside hotvr')
                    W_daughters_inside_hotvr = []
                    for W_daughter in closest_genW.daughters:
                        if deltaR(W_daughter, hotvr) < effective_radius:
                            W_daughters_inside_hotvr.append(W_daughter)
                    if self.print_out: print('Daugthers inside {}'.format(list(map(lambda daughter: daughter.pdgId, W_daughters_inside_hotvr))))
                    self.gen_W_substructures_check(hotvr_dict, W_daughters_inside_hotvr)
                else: hotvr_dict['has_other'] = True 

            elif len(genbs_not_from_top) != 0:
                closest_genb = min(genbs_not_from_top, key=lambda genb: deltaR(genb,hotvr))
                if deltaR(closest_genb,hotvr) < effective_radius:
                    if self.print_out: print('Closest b not from inside hotvr')
                    hotvr_dict['has_b_not_fromTop'] = True
                else: 
                    hotvr_dict['has_other'] = True 
                
            elif any(self.is_inside_hotvr(genp, hotvr) for genp in list(filter( self.is_quark_lepton_gluon, genparticles_not_from_top))):
                if self.print_out: print('Quarks/Gluons not from top inside hotvr')
                hotvr_dict['has_gluon_or_quark_not_fromTop'] = True

            else: 
                hotvr_dict['has_other'] = True 

            hotvr_df = pd.DataFrame(hotvr_dict, index=[hotvr_dict['run_number']])

            self.hotvrjets_df_list.append(hotvr_df)

        return True
