from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
from itertools import chain
from collections import OrderedDict
from itertools import izip_longest

ROOT.PyConfig.IgnoreCommandLineOptions = True

class ElectronSFProducer(Module):
    """Calculate electron scale factors
    """

    def __init__(
            self, sfFileName=None,
            inputElectronCollection = None,
            nosyst = False,
            verbose=0,
            outputName='',
            WP='',
            ID=''
    ):

        self.inputElectronCollection = inputElectronCollection
        self.nosyst = nosyst
        self.sfFileName = sfFileName
        self.verbose = verbose
        self.outputName = outputName
        self.WP = WP.capitalize()
        self.ID = ID

        self.sys = ['nominal', 'up', 'down']

    def beginJob(self):
        # initialize ElectronCorrlibReader
        self.corrlibreader = ROOT.ElectronCorrlibReader()
        self.corrlibreader.loadCorrections(self.sfFileName)
 
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for sys in self.sys:
            self.out.branch(self.outputName+"_weight_id_"+sys,"F")
            # for reco_type in ['PtAbove20','PtBelow20']:
            #     self.out.branch(self.outputName+"_weight_recoPt"+reco_type+"_"+sys,"F")    
            self.out.branch(self.outputName+"_weight_recoPt_"+sys,"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getSFs(self, reader, type, year, syst, wp, leptons):
        for idx, lep in enumerate(leptons):
            # evaluate SF
            sf = None
            if Module.globalOptions["year"] == '2022':
                year = '2022Re-recoBCD'
            if Module.globalOptions["year"] == '2022EE':
                year = '2022Re-recoE+PromptFG'
            sf = reader.evaluateElectronSF(type, year, syst, wp, lep.eta, lep.pt)
            # check if SF is OK
            
            if sf < 0.01:
                if self.verbose > 0:
                    print("electron #%i: pT = %1.1f, eta = %1.1f" % (
                        idx, lep.pt, lep.eta))
                sf = 1.
            yield sf

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        # jets = Collection(event, "Jet")
        corrlibreader = self.corrlibreader 
        electrons = self.inputElectronCollection(event)

        year = self.sfFileName.rsplit('/')[-2]

        if self.WP == 'Loose' and self.ID == 'MVA': self.WP = 'wploose'
        if self.WP == 'Medium' and self.ID == 'MVA': self.WP = 'wp90iso'
        if self.WP == 'Tight' and self.ID == 'MVA': self.WP = 'wp80iso'

        if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
            electron_label = 'Electron-ID-SF'
        else: electron_label = 'UL-Electron-ID-SF'

        for sys, sys_type in zip(self.sys, ['sf','sfup','sfdown']):
            # --- electron ID SF
            ev_weight_id, ev_weight_reco_pt = 1., 1. 
            if self.WP == 'wploose':  #SF not calculated for MVA WP LOOSE
                scale_factors_id = list(map(lambda electron: 1., electrons))
            else:
                scale_factors_id = list(self.getSFs(corrlibreader, electron_label, year, sys_type, self.WP, electrons))   

            for sf_id in scale_factors_id:
                ev_weight_id *= sf_id
            # ---

            # --- electron pT RECO SF
            electrons_pt_below20 = filter(lambda electron: electron.pt < 20., electrons)
            scale_factors_reco_pt_below20 = list(self.getSFs(corrlibreader, electron_label, year, sys_type, 'RecoBelow20', electrons_pt_below20))

            # differentiation between 2022, 2022EE
            # 2022, 2022EE pT ranges are: 1) below 20 pT; 2) from 20 to 75; 3) above 75
            # UL Run 2 pT ranges are: 1) below 20 pT; 2) above 75
            if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
                electrons_pt_from20to75 = filter(lambda electron: electron.pt >= 20. and electron.pt < 75., electrons) 
                scale_factors_reco_pt_from20to75 = list(self.getSFs(corrlibreader, electron_label, year, sys_type, 'Reco20to75', electrons_pt_from20to75))
                electrons_pt_above75 = filter(lambda electron: electron.pt >= 75., electrons) 
                scale_factors_reco_pt_above75 = list(self.getSFs(corrlibreader, electron_label, year, sys_type, 'RecoAbove75', electrons_pt_above75))

                for (sf_reco_pt_above75, sf_reco_pt_from20to75, sf_reco_pt_below20) in izip_longest(scale_factors_reco_pt_above75, scale_factors_reco_pt_from20to75, scale_factors_reco_pt_below20): 
                    ev_weight_reco_pt *= (sf_reco_pt_below20 or 1) * (sf_reco_pt_from20to75 or 1) * (sf_reco_pt_above75 or 1)

            else:
                electrons_pt_above20 = filter(lambda electron: electron.pt >= 20., electrons) 
                scale_factors_reco_pt_above20 = list(self.getSFs(corrlibreader, electron_label, year, sys_type, 'RecoAbove20', electrons_pt_above20))

                for (sf_reco_pt_above20, sf_reco_pt_below20) in izip_longest(scale_factors_reco_pt_above20, scale_factors_reco_pt_below20): 
                    ev_weight_reco_pt *= (sf_reco_pt_below20 or 1) * (sf_reco_pt_above20 or 1)

            self.out.fillBranch(self.outputName+"_weight_id_"+sys, ev_weight_id)
            self.out.fillBranch(self.outputName+"_weight_recoPt_"+sys, ev_weight_reco_pt)

        return True

