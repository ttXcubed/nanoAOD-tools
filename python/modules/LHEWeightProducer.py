from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class LHEWeightProducer(Module):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for yt in range(0, 601, 50):
            self.out.branch("weights_yt_{}".format(yt), "F")
        self.out.branch("weights_yt_sm", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #electrons = Collection(event, "Electron")
        #muons = Collection(event, "Muon")
        #jets = Collection(event, "Jet")
        #eventSum = ROOT.TLorentzVector()
        #for lep in muons:
        #    eventSum += lep.p4()
        #for lep in electrons:
        #    eventSum += lep.p4()
        #for j in filter(self.jetSel, jets):
        #    eventSum += j.p4()

        for yt in range(0, 601, 50):
            if hasattr(event, "LHEWeights_yt_{}".format(yt)):
                ytweight = getattr(event, "LHEWeights_yt_{}".format(yt))
            else: ytweight = 1
            self.out.fillBranch("weights_yt_{}".format(yt), ytweight)
        if hasattr(event, "LHEWeights_yt_sm"):
            ytweight = getattr(event, "LHEWeights_yt_sm")
        else: ytweight = 1
        self.out.fillBranch("weights_yt_sm", ytweight)
        return True