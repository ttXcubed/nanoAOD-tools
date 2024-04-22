from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os, re
from itertools import chain

from utils import getHist, getSFXY, processLabels

ROOT.PyConfig.IgnoreCommandLineOptions = True

def is_relevant_syst_for_shape_corr(flavor_btv, syst, jesSystsForShape=["jes"]):
    """Returns true if a flavor/syst combination is relevant"""
    jesSysts = list(chain(*[("up_" + j, "down_" + j)
                            for j in jesSystsForShape]))

    if flavor_btv == 5:
        return syst in ["central",
                        "up_lf", "down_lf",
                        "up_hfstats1", "down_hfstats1",
                        "up_hfstats2", "down_hfstats2"] + jesSysts
    elif flavor_btv == 4:
        return syst in ["central",
                        "up_cferr1", "down_cferr1",
                        "up_cferr2", "down_cferr2"]
    elif flavor_btv == 0:
        return syst in ["central",
                        "up_hf", "down_hf",
                        "up_lfstats1", "down_lfstats1",
                        "up_lfstats2", "down_lfstats2"] + jesSysts
    else:
        raise ValueError("ERROR: Undefined flavor = %i!!" % flavor_btv)
    return True


class btagSFProducer(Module):
    """Calculate btagging scale factors
    """

    def __init__(
            self, era, algo='deepJet', selectedWPs=['L','M', 'T','shape_corr'],
            sfFileName=None, verbose=0, jesSystsForShape=["jes"],
            nosyst = False, efficiencyMaps = ''
    ):
        self.era = era
        self.algo = algo
        self.selectedWPs = selectedWPs
        self.verbose = verbose
        self.jesSystsForShape = jesSystsForShape
        self.nosyst = nosyst
        self.max_abs_eta = 2.5
        # define measurement type for each flavor
        self.sfFileName = sfFileName
        self.efficiencyMaps = efficiencyMaps
        self.measurement_types = None
        self.supported_wp = None
        supported_btagSF = {
            'deepJet': {
                '2016preVFP': {
                    'inputFileName': "2016preVFP_UL/btagging.json.gz",
                    'measurement_types': {
                        5: "comb",  # b
                        4: "comb",  # c
                        0: "incl"   # light
                    },
                    'supported_wp': ["L", "M", "T", "shape_corr"]
                },
                # update file name once available!
                '2016': {
                    'inputFileName': "2016postVFP_UL/btagging.json.gz",
                    'measurement_types': {
                        5: "comb",  # b
                        4: "comb",  # c
                        0: "incl"   # light
                    },
                    'supported_wp': ["L", "M", "T", "shape_corr"]
                },
                '2017': {
                    'inputFileName': "2017_UL/btagging.json.gz",
                    'measurement_types': {
                        5: "comb",  # b
                        4: "comb",  # c
                        0: "incl"   # light
                    },
                    'supported_wp': ["L", "M", "T", "shape_corr"]
                },
                '2018': {
                    'inputFileName': "2018_UL/btagging.json.gz",
                    'measurement_types': {
                        5: "comb",  # b
                        4: "comb",  # c
                        0: "incl"   # light
                    },
                    'supported_wp': ["L", "M", "T", "shape_corr"]
                },
                '2022': {
                    'inputFileName': "2022/btagging.json.gz",
                    'measurement_types': {
                        5: "comb",  # b
                        4: "comb",  # c
                        0: "light"   # light
                    },
                    'supported_wp': ["L", "M", "T", "shape_corr"]
                },
                '2022EE': {
                    'inputFileName': "2022EE/btagging.json.gz",
                    'measurement_types': {
                        5: "comb",  # b
                        4: "comb",  # c
                        0: "light"   # light
                    },
                    'supported_wp': ["L", "M", "T", "shape_corr"]
                },
            },
        }

        supported_algos = []
        for algo in list(supported_btagSF.keys()):
            if self.era in list(supported_btagSF[algo].keys()):
                supported_algos.append(algo)
        if self.algo in list(supported_btagSF.keys()):
            if self.era in list(supported_btagSF[self.algo].keys()):
                self.sfFileName = self.sfFileName + '/' + supported_btagSF[self.algo][self.era]['inputFileName']
                self.measurement_types = supported_btagSF[self.algo][self.era]['measurement_types']
                self.supported_wp = supported_btagSF[self.algo][self.era]['supported_wp']
            else:
                raise ValueError("ERROR: Algorithm '%s' not supported for era = '%s'! Please choose among { %s }." % (
                    self.algo, self.era, supported_algos))
        else:
            raise ValueError("ERROR: Algorithm '%s' not supported for era = '%s'! Please choose among { %s }." % (
                self.algo, self.era, supported_algos))
        for wp in self.selectedWPs:
            if wp not in self.supported_wp:
                raise ValueError("ERROR: Working point '%s' not supported for algo = '%s' and era = '%s'! Please choose among { %s }." % (
                    wp, self.algo, self.era, self.supported_wp))
        
        # define systematic uncertainties
        self.systs = []
        # --- FIXED WP

        # single data era analysis
        # self.central_and_systs = ["central"]
        # self.systs.append("up")
        # self.systs.append("down")
        # if not self.nosyst:
        #     self.central_and_systs.extend(self.systs_shape_corr)

        # multiple data era analysis
        self.central_and_systs = ["central"]
        if not self.nosyst:
            for syst in ['correlated', 'uncorrelated']: # additional uncertainties for fixed b-tagging WP: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation#Correlations_across_years_2016_2
                self.central_and_systs.append("up_%s" % syst)
                self.central_and_systs.append("down_%s" % syst)
        # ---

        self.systs_shape_corr = []
        for syst in ['lf', 'hf',
                     'hfstats1', 'hfstats2',
                     'lfstats1', 'lfstats2',
                     'cferr1', 'cferr2'] + self.jesSystsForShape:
            self.systs_shape_corr.append("up_%s" % syst)
            self.systs_shape_corr.append("down_%s" % syst)
        self.central_and_systs_shape_corr = ["central"]
        if not self.nosyst:
            self.central_and_systs_shape_corr.extend(self.systs_shape_corr)

        self.branchNames_central_and_systs = {}
        for wp in self.selectedWPs:
            branchNames = {}
            if wp == 'shape_corr':
                central_and_systs = self.central_and_systs_shape_corr
                baseBranchName = 'btagEventWeight_{}_shape'.format(self.algo)
                for central_or_syst in central_and_systs:
                    branchNames[self.getSystForFwk(central_or_syst)] = baseBranchName + '_' + self.getSystForFwk(central_or_syst)
            else:
                central_and_systs = self.central_and_systs
                for central_or_syst in central_and_systs:
                    branchNames[self.getSystForFwk(central_or_syst)] = {}
                    for flav in ['bc', 'light']:
                        baseBranchName = 'btagSF{}_{}_{}'.format(flav, self.algo, wp)
                        branchNames[self.getSystForFwk(central_or_syst)][flav] = baseBranchName + '_' + self.getSystForFwk(central_or_syst)
            self.branchNames_central_and_systs[wp] = branchNames

        self.WP_LABEL = {
            'T': 'tight', 'M': 'medium', 'L': 'loose'
        }
        self.efficiencyMapHist = {}
        for wp in self.selectedWPs:
            if wp != 'shape_corr':
                self.efficiencyMapHist[self.WP_LABEL[wp]] = {}
                for flavor in ['b', 'c', 'usdg']:
                    self.efficiencyMapHist[self.WP_LABEL[wp]][flavor] = ROOT.TH2F()

    def beginJob(self):
        # initialize BTagCorrlibReader
        self.corrlibreader = ROOT.BTagCorrlibReader()
        if not os.path.exists(self.sfFileName): 
            raise ValueError("No SF files for b-tagging found in {}".format(self.sfFileName))
        self.corrlibreader.loadCorrections(os.path.join(self.sfFileName))
        self.readers = {}
        for wp in self.selectedWPs:
            wp_btv = {"l": "L", "m": "M", "t": "T",
                      "shape_corr": "shape_corr"}.get(wp.lower(), None)
            syts = None
            if wp_btv in ["L", "M", "T"]:
                systs = self.systs
            else:
                systs = self.systs_shape_corr

            self.readers[wp_btv] = self.corrlibreader

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # --- this procedure to parse the process depends strongly on the output/input directory! modify it accordingly
        process = ''
        pattern = re.compile(r"/([^/]*?)_TuneCP5")
        match = pattern.search(str(inputFile.GetName()))
        if match: 
            if processLabels(match.group(1)) == None:
                process = match.group(1)
            else: process = processLabels(match.group(1))
            print('Process: {}'.format(process))
        else:
            pattern = re.compile(r"/([^/]*?)_MC")
            match = pattern.search(str(outputFile.GetName()))
            if match:
                process = match.group(1)
            else:
                print('Process name not parsed.')
        # ---

        eff_file_path = os.path.join(self.efficiencyMaps, '{}_efficiencyMap.root'.format(process))
        if not os.path.exists(eff_file_path): 
            eff_file_path = os.path.join(self.efficiencyMaps, 'tt_dilepton_efficiencyMap.root'.format(process))
            print('No b-tagging efficiency map found for the process --> using tt_dilepton efficiency map...')

        for wp in self.selectedWPs:
            if wp != 'shape_corr':
                self.efficiencyMapHist[self.WP_LABEL[wp]] = {}
                for flavor in ['b', 'c', 'udsg']:
                    self.efficiencyMapHist[self.WP_LABEL[wp]][flavor] = getHist(
                        eff_file_path,
                        'efficiency_{}_{}'.format(flavor, self.WP_LABEL[wp])
                    )

        # very low flexibility --> it needs to be changed when both shape and fixed WP are required
        if self.selectedWPs[0] == 'shape_corr' and len(self.selectedWPs) == 1:
            for central_or_syst in list(self.branchNames_central_and_systs.values()):
                for branch in list(central_or_syst.values()):
                    # self.out.branch(branch, "F", lenVar="nJet")
                    self.out.branch(branch, "F")
        else:
            for central_or_syst in list(self.branchNames_central_and_systs.values()):
                for systematics_branches in list(central_or_syst.values()):
                    for branch in list(systematics_branches.values()):
                        self.out.branch(branch, "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getReader(self, wp):
        """
            Get btag scale factor reader.
        """
        wp_btv = {"l": "L", "m": "M", "t": "T",
                  "shape_corr": "shape_corr"}.get(wp.lower(), None)
        if wp_btv is None or wp_btv not in list(self.readers.keys()):
            if self.verbose > 0:
                print(
                    "WARNING: Unknown working point '%s', setting b-tagging SF reader to None!" % wp)
            return None
        return self.readers[wp_btv]

    def getSFs(self, jet_data, syst, reader, measurement_type='auto', wp='', shape_corr=False):
        if reader is None:
            if self.verbose > 0:
                print("WARNING: Reader not available, setting b-tagging SF to -1!")
            for i in range(len(jet_data)):
                yield 1
            raise StopIteration
        for idx, (pt, eta, flavor_btv, discr) in enumerate(jet_data):
            epsilon = 1.e-3
            max_abs_eta = self.max_abs_eta
            if eta <= -max_abs_eta:
                eta = -max_abs_eta + epsilon
            if eta >= +max_abs_eta:
                eta = +max_abs_eta - epsilon
            # evaluate SF
            sf = None
            if shape_corr:
                if is_relevant_syst_for_shape_corr(flavor_btv, syst, self.jesSystsForShape):
                    sf = reader.evaluateBTagShape(self.algo+"_shape",
                        syst, flavor_btv, abs(eta), pt, discr)
                else:
                    sf = reader.evaluateBTagShape(self.algo+"_shape",
                        'central', flavor_btv, abs(eta), pt, discr)
            else:
                if self.measurement_types[flavor_btv] == 'light' and 'correlated' in syst: #no up/down correlated SF for 2022, 2022EE -- 07/03/24
                    sf = 1.
                else:
                    # print(self.algo+"_"+self.measurement_types[flavor_btv],syst, wp, flavor_btv, abs(eta), pt)
                    sf = reader.evaluateBTagWorkingpoint(self.algo+"_"+self.measurement_types[flavor_btv],syst, wp, flavor_btv, abs(eta), pt)
            # check if SF is OK
            if sf < 0.01:
                if self.verbose > 0:
                    print("jet #%i: pT = %1.1f, eta = %1.1f, discr = %1.3f, flavor = %i" % (
                        idx, pt, eta, discr, flavor_btv))
                sf = 1.
            yield sf

    def getEff(self, jet_data, syst, wp=''):
        # retrieve efficiency from etaVSpt map calculated for each jet flavor/wp
        for idx, (pt, eta, flavor_btv, discr) in enumerate(jet_data):
            epsilon = 1.e-3
            max_abs_eta = self.max_abs_eta
            if eta <= -max_abs_eta:
                eta = -max_abs_eta + epsilon
            if eta >= +max_abs_eta:
                eta = +max_abs_eta - epsilon

            if flavor_btv == 5: flavor = 'b'
            elif flavor_btv == 4: flavor = 'c'
            else: flavor = 'udsg'

            efficiency, efficiency_err = getSFXY(self.efficiencyMapHist[wp][flavor], pt, abs(eta))
            if efficiency == 0.:
                if self.verbose > 0:
                    print("jet #%i: pT = %1.1f, eta = %1.1f, discr = %1.3f, flavor = %i" % (
                        idx, pt, eta, discr, flavor_btv))
                    print("No efficiency (pT, eta) available...")
                efficiency, efficiency_err = 1., 1.
            yield efficiency

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        # jets = Collection(event, "Jet")

        discr = None
        if self.algo == "csvv2":
            discr = "btagCSVV2"
        elif self.algo == "deepcsv":
            discr = "btagDeepB"
        elif self.algo == "cmva":
            discr = "btagCMVA"
        elif self.algo == "deepJet":
            discr = "btagDeepFlavB"
        else:
            raise ValueError("ERROR: Invalid algorithm '%s'!" % self.algo)

        for wp in self.selectedWPs:
            reader = self.getReader(wp)
            isShape = (wp == 'shape_corr')
            central_and_systs = (
                self.central_and_systs_shape_corr if isShape else self.central_and_systs)
            for central_or_syst in central_and_systs:
                jet_collection = getattr(event,"selectedJets_nominal")
                if self.isJESvariation(central_or_syst):
                    jet_collection = getattr(event,"selectedJets_"+self.getSystForFwk(central_or_syst))
                
                preloaded_jets = [(jet.pt, jet.eta, jet.hadronFlavour, getattr(jet, discr)) 
                                    for jet in jet_collection]
                # print(central_or_syst, )

                if isShape:
                    ev_weight = 1.
                    scale_factors = list(self.getSFs(
                    preloaded_jets, central_or_syst, reader, 'auto', wp, isShape))
                    for SF in scale_factors:
                        ev_weight *= SF

                else: 
                    ev_weight = {'bc': 1., 'light': 1.}
                    efficiencies = list(self.getEff(
                        preloaded_jets, central_or_syst, self.WP_LABEL[wp]))
                    scale_factors = list(self.getSFs(
                    preloaded_jets, central_or_syst, reader, 'auto', wp))

                    data_weight = {'bc': 1., 'light': 1.}
                    mc_weight = {'bc': 1., 'light': 1.} 
                    # following the definition here: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#b_tagging_efficiency_in_MC_sampl
                    # and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation#Correlations_across_years_2016_2

                    for ijet, jet in enumerate(jet_collection): 
                        flav = 'bc'
                        if jet.hadronFlavour < 4: flav = 'light'

                        if getattr(jet, "b_tagged_{}".format(self.WP_LABEL[wp])):
                            data_weight[flav] *= (scale_factors[ijet] * efficiencies[ijet])
                            mc_weight[flav] *= efficiencies[ijet]
                        else:
                            data_weight[flav] *= (1. - (scale_factors[ijet] * efficiencies[ijet]))
                            mc_weight[flav] *= (1. - efficiencies[ijet])

                    for flav in ['bc', 'light']:
                        if mc_weight[flav] != 0.:
                            ev_weight[flav] = data_weight[flav] / mc_weight[flav]
                        else: ev_weight[flav] = 1.

                        self.out.fillBranch(self.branchNames_central_and_systs[wp][self.getSystForFwk(central_or_syst, wp)][flav], ev_weight[flav])

        return True

    def isJESvariation(self,central_or_syst):
        for syst in self.jesSystsForShape:
            if syst in central_or_syst:
                return True
        return False
        
    def getSystForFwk(self, syst, wp='shape_coor'):
        if syst == 'central':
            if wp != 'shape_corr':
                return Module.globalOptions['year'] 
            else:
                return 'nominal'
        if 'fstats' in syst:
            syst += '_' + Module.globalOptions['year']
        if syst.startswith('up_'):
            syst = syst.replace('up_','') + '_up'
        elif syst.startswith('down_'):
            syst = syst.replace('down_','') + '_down'
        if syst == 'jesUp': syst = 'jesTotalUp'
        elif syst == 'jesDown': syst = 'jesTotalDown'
        if syst == 'uncorrelated_up': syst = Module.globalOptions['year'] + '_up'
        if syst == 'uncorrelated_down': syst = Module.globalOptions['year'] + '_down'
        return syst

