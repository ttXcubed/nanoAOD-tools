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
            inputMuonCollection = None,
            nosyst = False,
            verbose=0,
            outputName = '',
            WP = '',
    ):

        self.inputMuonCollection = inputMuonCollection
        self.nosyst = nosyst
        self.sfFileName = sfFileName
        self.verbose = verbose
        self.outputName = outputName
        self.WP = WP

        self.sys = ['nominal', 'up', 'down']

    def beginJob(self):
        # initialize MuonCorrlibReader
        self.corrlibreader = ROOT.MuonCorrlibReader()
        self.corrlibreader.loadCorrections(self.sfFileName)
 
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for sys in self.sys:
            self.out.branch(self.outputName+"_weight_id_"+sys,"F")
            self.out.branch(self.outputName+"_weight_iso_"+sys,"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

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
        muons = self.inputMuonCollection(event)

        year = self.sfFileName.rsplit('/')[-2]

        wp = ''
        if self.WP=='loose': wp='Loose'
        if self.WP=='medium': wp='Medium'
        if self.WP=='tight': wp='Tight'

        for sys, sys_type in zip(self.sys, ['sf','systup','systdown']):
            ev_weight_id, ev_weight_iso = 1., 1.
            scale_factors_id = list(self.getSFs(corrlibreader, 'NUM_'+wp+'ID_DEN_genTracks', year, muons, sys_type))

            if self.WP=='loose': 
                scale_factors_iso = [1.]
            elif self.WP=='tight':
                scale_factors_iso = list(self.getSFs(corrlibreader, 'NUM_TightRelIso_DEN_'+wp+'IDandIPCut', year, muons, sys_type))    
            else: 
                scale_factors_iso = list(self.getSFs(corrlibreader, 'NUM_TightRelIso_DEN_'+wp+'ID', year, muons, sys_type)) 
                
            for sf_id, sf_iso in zip(scale_factors_id, scale_factors_iso):
                ev_weight_id *= sf_id
                ev_weight_iso *= sf_iso
            self.out.fillBranch(self.outputName+"_weight_id_"+sys, ev_weight_id)
            self.out.fillBranch(self.outputName+"_weight_iso_"+sys, ev_weight_iso)

        return True

 