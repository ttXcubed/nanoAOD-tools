from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
from itertools import chain
from collections import OrderedDict

ROOT.PyConfig.IgnoreCommandLineOptions = True

class ElectronSFProducer(Module):
    """Calculate electron scale factors
    """

    def __init__(
            self, sfFileName=None,
            inputElectronCollection = OrderedDict(),
            nosyst = False,
            verbose=0,
    ):

        self.inputElectronCollection = inputElectronCollection
        self.nosyst = nosyst
        self.sfFileName = sfFileName
        self.verbose = verbose

        self.id_type = ['MVA', 'cutBased']
        self.sys = ['nominal', 'up', 'down']
        self.wp = OrderedDict([('MVA', ['wp80iso','wp90iso','wploose']), ('cutBased', ['Tight','Medium','Loose'])]) 

    def beginJob(self):
        # initialize ElectronCorrlibReader
        self.corrlibreader = ROOT.ElectronCorrlibReader()
        self.corrlibreader.loadCorrections(self.sfFileName)
 
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for id_type in self.id_type:
            for sys in self.sys:
                for wp in ['tight', 'medium', 'loose']:
                    self.out.branch(wp+"_"+id_type+"_Electrons_weight_id_"+sys,"F")
                    for reco_type in ['PtAbove20','PtBelow20']:
                        self.out.branch(wp+"_"+id_type+"_Electrons_weight_reco"+reco_type+"_"+sys,"F")    

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getSFs(self, reader, type, year, syst, wp, leptons):
        for idx, lep in enumerate(leptons):
            # evaluate SF
            sf = None
            #print(type, abs(lep.eta), lep.pt, syst)
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
        electrons_collections = self.inputElectronCollection

        year = self.sfFileName.rsplit('/')[-2]

        for id_type in self.id_type:
            for sys, sys_type in zip(self.sys, ['sf','sfup','sfdown']):
               for wp, wp_type in zip(self.wp[id_type], ['tight','medium','loose']): 
                    ev_weight_id, ev_weight_reco_pt_above20, ev_weight_reco_pt_below20 = 1., 1., 1.
                    if wp=='wploose':  #SF not calculated for MVA WP LOOSE
                        scale_factors_id, scale_factors_reco_pt_above20, scale_factors_reco_pt_below20 = [1.], [1.], [1.]
                    else:
                        electrons=electrons_collections[id_type][wp_type](event)
                        scale_factors_id = list(self.getSFs(corrlibreader, 'UL-Electron-ID-SF', year, sys_type, wp, electrons))   
                        electrons_pt_above20 = filter(lambda electron: electron.pt>=20., electrons) 
                        scale_factors_reco_pt_above20 = list(self.getSFs(corrlibreader, 'UL-Electron-ID-SF', year, sys_type, 'RecoAbove20', electrons_pt_above20))
                        electrons_pt_below20 = filter(lambda electron: electron.pt<20., electrons)
                        scale_factors_reco_pt_below20 = list(self.getSFs(corrlibreader, 'UL-Electron-ID-SF', year, sys_type, 'RecoBelow20', electrons_pt_below20))
                    
                        for sf_id, sf_reco_pt_above20, sf_reco_pt_below20 in zip(scale_factors_id, scale_factors_reco_pt_above20, scale_factors_reco_pt_below20): 
                            ev_weight_id *= sf_id
                            ev_weight_reco_pt_above20 *= sf_reco_pt_above20
                            ev_weight_reco_pt_below20 *= sf_reco_pt_below20
                    self.out.fillBranch(wp_type+"_"+id_type+"_Electrons_weight_id_"+sys, ev_weight_id)
                    self.out.fillBranch(wp_type+"_"+id_type+"_Electrons_weight_recoPtAbove20_"+sys, ev_weight_reco_pt_above20)
                    self.out.fillBranch(wp_type+"_"+id_type+"_Electrons_weight_recoPtBelow20_"+sys, ev_weight_reco_pt_below20)

        return True

