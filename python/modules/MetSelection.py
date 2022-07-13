from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class MetSelection(Module):
    def __init__(
         self,
         outputName="MET",
         storeKinematics=['pt', 'phi']
     ):
        self.outputName = outputName
        self.storeKinematics = storeKinematics

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable, "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName + "_" + variable, getattr(event, "MET_"+variable))

        return True