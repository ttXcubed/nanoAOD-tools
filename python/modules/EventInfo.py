import os
import sys
import math
import json
import ROOT
import random
import yaml
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Event, InputTree
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class EventInfo(Module):
    def __init__(
        self,
        storeVariables = [],
    ):
        self.storeVariables = storeVariables

        self.genEventSumw = 0
        self.genEventCount = 0
        self.genEventSumw2 = 0
        self.LHEScaleSumw = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        self.genEventInfoInFile = False

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for variable in self.storeVariables:
            variable[0](self.out)

        runsTree = InputTree(inputFile.Get("Runs"))
        runsEvent = Event(runsTree, 0)
        if hasattr(runsEvent, "genEventSumw"):
            self.genEventInfoInFile = True

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        if not Module.globalOptions["isData"]:            
            if self.genEventInfoInFile == True:
                runsTree = InputTree(inputFile.Get("Runs"))
                runsEvent = Event(runsTree, 0)
                self.genEventSumw = runsEvent.genEventSumw
                self.genEventCount = runsEvent.genEventCount
                self.genEventSumw2 = runsEvent.genEventSumw2
                self.LHEScaleSumw = [runsEvent.LHEScaleSumw[i] for i in range(9)]

            nGenEventCount_parameter = ROOT.TParameter(float)("genEventCount", self.genEventCount)
            nGenWeight_parameter = ROOT.TParameter(float)("sumGenWeights", self.genEventSumw)
            genEventSumw2_parameter = ROOT.TParameter(float)("sumGenWeights2", self.genEventSumw2)
            LHEScaleSumw_parameter = []
            # -- in case the parameters are not already stored in "Runs" TTree of the input file,
            # -- a manual calculation is required --> https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2017UL/doc_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1.html

            if self.genEventInfoInFile == True:
                for i in range(9): LHEScaleSumw_parameter.append(ROOT.TParameter(float)("LHEScaleSumw_{}".format(i), self.LHEScaleSumw[i]))
            else:
                for i in range(9): LHEScaleSumw_parameter.append(ROOT.TParameter(float)("LHEScaleSumw_{}".format(i), self.LHEScaleSumw[i]/self.genEventSumw))

            outputFile.cd()
            nGenEventCount_parameter.Write()
            nGenWeight_parameter.Write() 
            genEventSumw2_parameter.Write()
            for i in range(9): LHEScaleSumw_parameter[i].Write()
        else:
            pass
        
    def analyze(self, event):

        if not Module.globalOptions["isData"]:
            # -- in case the parameters are not already stored in "Runs" TTree of the input file,
            # -- a manual calculation is required
            if self.genEventInfoInFile == False:
                self.genEventSumw += event.Generator_weight
                self.genEventCount += 1
                self.genEventSumw2 += event.Generator_weight ** 2
                for i in range(9): self.LHEScaleSumw[i] += event.LHEScaleWeight[i] * event.Generator_weight

        else:
            with open(os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/13TeV_UL_Era_runNumber.yaml") as yaml_f:
                try:
                    era_yaml = yaml.safe_load(yaml_f)
                except yaml.YAMLError as exc:
                    print(exc)

            for era in era_yaml[Module.globalOptions['year']]:
                if event.run >= era_yaml[Module.globalOptions['year']][era][0] and event.run <= era_yaml[Module.globalOptions['year']][era][1]:
                    Module.globalOptions['era'] = era

        for variable in self.storeVariables:
            variable[1](self.out, event)
        return True
