import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi
from collections import OrderedDict

class BTagSelection(Module):
    #tight DeepFlav WP (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy)
    LOOSE=0
    MEDIUM=1
    TIGHT=2
    
    def __init__(
        self,
        btaggingWP = {},
        inputCollection=lambda event: Collection(event, "Jet"),
        outputName= 'selectedBJets_', #"btaggedJets",
        jetMinPt=30.,
        jetMaxEta=2.4,
        storeKinematics=['pt', 'eta','phi','mass'],
        storeTruthKeys=[]
    ):
        self.btaggingWP = btaggingWP
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.jetMinPt = jetMinPt
        self.jetMaxEta = jetMaxEta
        self.storeKinematics = storeKinematics
        self.storeTruthKeys = storeTruthKeys
        self.working_points = btaggingWP.keys()

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for wp in self.working_points:
            self.out.branch("n{}_{}".format(self.outputName, wp), "I")
            for variable in self.storeKinematics:
                self.out.branch("{}_{}_{}".format(self.outputName, wp, variable), "F", lenVar="n{}_{}".format(self.outputName, wp))
            if not Module.globalOptions["isData"]:
                for variable in self.storeTruthKeys:
                    self.out.branch("{}_{}_{}".format(self.outputName, wp, variable), "F", lenVar="n{}_{}".format(self.outputName, wp))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        # only info of b-tagged jets will be stored in the ntuplizer
        bJets = OrderedDict([('tight',[]), ('medium',[]), ('loose',[])])
        lJets = []

        for jet in jets:
            # probably this check is not needed since the input collection
            # is the selectedJets which have already this pT, eta requirements
            if jet.pt<self.jetMinPt:
                lJets.append(jet)
                continue
            if math.fabs(jet.eta) > self.jetMaxEta:
                lJets.append(jet)
                continue

            if jet.btagDeepFlavB > self.btaggingWP['tight']:
                bJets['tight'].append(jet), bJets['medium'].append(jet), bJets['loose'].append(jet)
                setattr(jet, "b_tagged_tight", True)
                setattr(jet, "b_tagged_medium", True)  # A tight jet is also medium
                setattr(jet, "b_tagged_loose", True)   # And also loose
            elif jet.btagDeepFlavB > self.btaggingWP['medium']:
                bJets['medium'].append(jet), bJets['loose'].append(jet)
                setattr(jet, "b_tagged_medium", True)
                setattr(jet, "b_tagged_loose", True)   # A medium jet is also loose
                setattr(jet, "b_tagged_tight", False)
            elif jet.btagDeepFlavB > self.btaggingWP['loose']:
                bJets['loose'].append(jet)
                setattr(jet, "b_tagged_loose", True)
                setattr(jet, "b_tagged_medium", False)
                setattr(jet, "b_tagged_tight", False)
            else:
                setattr(jet, "b_tagged_loose", False)
                setattr(jet, "b_tagged_medium", False)
                setattr(jet, "b_tagged_tight", False)


        for wp in self.working_points:
            self.out.fillBranch("n{}_{}".format(self.outputName, wp), len(bJets[wp]))
            for variable in self.storeKinematics:
                self.out.fillBranch("{}_{}_{}".format(self.outputName, wp, variable), map(lambda jet: getattr(jet, variable), bJets[wp]))

            if not Module.globalOptions["isData"]:
                for variable in self.storeTruthKeys:
                    self.out.fillBranch("{}_{}_{}".format(self.outputName, wp, variable), map(lambda jet: getattr(jet, variable), bJets[wp]))
            setattr(event, self.outputName + '_' + wp, bJets[wp])

        return True
