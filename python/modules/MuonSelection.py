import math
import sys

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR

class MuonSelection(Module):
    VERYTIGHT = 1
    TIGHT = 1
    MEDIUM = 2
    LOOSE = 3
    NONE = 4
    INV = 5

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Muon"),
        outputName="tightMuons",
        triggerMatch=False,
        muonID=TIGHT,
        muonIso=TIGHT,
        muonMinPt=25.,
        muonMaxEta=2.4,
        storeKinematics=['pt','eta'],
        storeWeights=False,
    ):
        
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.muonMinPt = muonMinPt
        self.muonMaxEta = muonMaxEta
        self.storeKinematics = storeKinematics
        self.storeWeights = storeWeights
        self.triggerMatch = triggerMatch

        if muonID==MuonSelection.MEDIUM or muonIso==MuonSelection.MEDIUM:
            raise Exception("Unsupported ID or ISO")

        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []

        if muonID==MuonSelection.TIGHT:
            self.muonIdFct = lambda muon: muon.tightId==1
            #self.muonIdSF = self.idTightSFHist
        elif muonID==MuonSelection.LOOSE:
            self.muonIdFct = lambda muon: muon.looseId==1
            #self.muonIdSF = self.idLooseSFHist
        elif muonID==MuonSelection.NONE:
            self.muonIdFct = lambda muon: True
            #self.muonIdSF = self.idLooseSFHist        

    def triggerMatched(self, muon, trigger_object):
        if self.triggerMatch:
            trig_deltaR = math.pi
            for trig_obj in trigger_object:
                if abs(trig_obj.id) != 13:
                    continue
                trig_deltaR = min(trig_deltaR, deltaR(trig_obj, muon))
            if trig_deltaR < 0.3:
                return True, trigger_object._idx
            else:
                return False, None
        else:
            return True, None    
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = self.inputCollection(event)

        triggerObjects = self.triggerObjectCollection(event)
        matched_trgObjs_id = []

        selectedMuons = []
        unselectedMuons = []
               
        #https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
        for muon in muons:

            # trigger matching
            hasTriggerObj, matched_trgObjs_id = self.triggerMatched(muon, triggerObjects)
            if len(muons)>1: 
                print(self.triggerMatched(muon, triggerObjects))
                sys.exit(1)
            # matched_trgObjs_id.append(self.triggerMatched(muon, triggerObjects)[1])

            if muon.pt>self.muonMinPt \
            and math.fabs(muon.eta)<self.muonMaxEta \
            and self.muonIdFct(muon) \
            and self.muonIsoFct(muon) \
            and self.triggerMatched(muon, triggerObjects):
            
                selectedMuons.append(muon)
                
            else:
                unselectedMuons.append(muon)


        self.out.fillBranch("n"+self.outputName,len(selectedMuons))
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,map(lambda muon: getattr(muon,variable),selectedMuons))
        
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_id_nominal", weight_id_nominal)
            self.out.fillBranch(self.outputName+"_weight_id_up", weight_id_up)
            self.out.fillBranch(self.outputName+"_weight_id_down", weight_id_down)
            
            self.out.fillBranch(self.outputName+"_weight_iso_nominal", weight_iso_nominal)
            self.out.fillBranch(self.outputName+"_weight_iso_up", weight_iso_up)
            self.out.fillBranch(self.outputName+"_weight_iso_down", weight_iso_down)

        setattr(event,self.outputName,selectedMuons)
        setattr(event,self.outputName+"_unselected",unselectedMuons)

        return True

