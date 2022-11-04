from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
from itertools import chain
from collections import OrderedDict

ROOT.PyConfig.IgnoreCommandLineOptions = True

class MuonSFProducer(Module):
    """Calculate muon scale factors
    """

    def __init__(
            self, sfFileName=None,
            inputMuonCollection = OrderedDict(),
            nosyst = False,
            verbose=0,
    ):

        self.inputMuonCollection = inputMuonCollection
        self.nosyst = nosyst
        self.sfFileName = sfFileName
        self.verbose = verbose

        self.sys = ['nominal', 'up', 'down']
        self.wp = ['Tight','Medium','Loose']

    def beginJob(self):
        # initialize MuonCorrlibReader
        self.corrlibreader = ROOT.MuonCorrlibReader()
        self.corrlibreader.loadCorrections(self.sfFileName)
 
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for sys in self.sys:
            for wp in ['tight','medium','loose']:
                self.out.branch("tightRelIso_"+wp+"ID_Muons_weight_id_"+sys,"F")
                self.out.branch("tightRelIso_"+wp+"ID_Muons_weight_iso_"+sys,"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    # def getReader(self, wp):
    #     """
    #         Get btag scale factor reader.
    #     """
    #     wp_btv = {"l": "L", "m": "M", "t": "T",
    #               "shape_corr": "shape_corr"}.get(wp.lower(), None)
    #     if wp_btv is None or wp_btv not in list(self.readers.keys()):
    #         if self.verbose > 0:
    #             print(
    #                 "WARNING: Unknown working point '%s', setting b-tagging SF reader to None!" % wp)
    #         return None
    #     return self.readers[wp_btv]

    def getSFs(self, reader, type, year, leptons, syst):
        for idx, lep in enumerate(leptons):
            # evaluate SF
            sf = None
            #print(type, abs(lep.eta), lep.pt, syst)
            sf = reader.evaluateMuonSF(type, year, abs(lep.eta), lep.pt, syst)
            # check if SF is OK
            if sf < 0.01:
                if self.verbose > 0:
                    print("muon #%i: pT = %1.1f, eta = %1.1f" % (
                        idx, lep.pt, lep.eta))
                sf = 1.
            yield sf

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        # jets = Collection(event, "Jet")
        corrlibreader = self.corrlibreader 
        muons_collections = self.inputMuonCollection

        year = self.sfFileName.rsplit('/')[-2]

        for sys, sys_type in zip(self.sys, ['sf','systup','systdown']):
            for wp, wp_type in zip(self.wp,['tight','medium','loose']):
                muons=muons_collections[wp_type](event)
                ev_weight_id, ev_weight_iso = 1., 1.

                scale_factors_id = list(self.getSFs(corrlibreader, 'NUM_'+wp+'ID_DEN_genTracks', year, muons, sys_type))
                if wp_type=='loose': 
                    scale_factors_iso = [1.]
                elif wp_type=='tight':
                    scale_factors_iso = list(self.getSFs(corrlibreader, 'NUM_TightRelIso_DEN_'+wp+'IDandIPCut', year, muons, sys_type))    
                else: 
                    scale_factors_iso = list(self.getSFs(corrlibreader, 'NUM_TightRelIso_DEN_'+wp+'ID', year, muons, sys_type)) 
                
                for sf_id, sf_iso in zip(scale_factors_id, scale_factors_iso):
                    ev_weight_id *= sf_id
                    ev_weight_iso *= sf_iso
                self.out.fillBranch("tightRelIso_"+wp_type+"ID_Muons_weight_id_"+sys, ev_weight_id)
                self.out.fillBranch("tightRelIso_"+wp_type+"ID_Muons_weight_iso_"+sys, ev_weight_iso)

        return True

 