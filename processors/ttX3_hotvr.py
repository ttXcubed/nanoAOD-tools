import argparse
import json
import math
import os
import random
import sys
from collections import OrderedDict

import numpy as np
import ROOT
from PhysicsTools.NanoAODTools.modules import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import (
    Collection, Object)
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import \
    PostProcessor

parser = argparse.ArgumentParser()


parser.add_argument('--isData', dest='isData',
                    action='store_true', default=False)
parser.add_argument('--isSignal', dest='isSignal',
                    action='store_true', default=False)
parser.add_argument('--nosys', dest='nosys',
                    action='store_true', default=False)
parser.add_argument('--invid', dest='invid',
                    action='store_true', default=False)
parser.add_argument('--year', dest='year',
                    action='store', type=str, default='2017', choices=[
                        '2016','2016preVFP','2017','2018', '2022', '2022EE'])
parser.add_argument('-i','--input', dest='inputFiles', action='append', default=[])
parser.add_argument('--maxEvents', dest='maxEvents', type=int, default=None)
parser.add_argument('--trigger', dest='trigger', type=str, default=None, choices=['mumu','emu','ee'])
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "isData:",args.isData
print "isSignal:",args.isSignal
print "evaluate systematics:",not args.nosys
print "invert lepton id/iso:",args.invid
print "inputs:",len(args.inputFiles)
print "year:", args.year
print "output directory:", args.output[0]
if args.maxEvents:
    print 'max number of events', args.maxEvents

globalOptions = {
    "isData": args.isData,
    "isSignal": args.isSignal,
    "year": args.year,
    "trigger": args.trigger, 
    "era": None
}

Module.globalOptions = globalOptions

isMC = not args.isData
isPowheg = 'powheg' in args.inputFiles[0].lower()
isPowhegTTbar = 'TTTo' in args.inputFiles[0] and isPowheg

#recommended pT threshold for the subleading lepton --> https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2020/085 (4 top in dilepton final state)
minMuonPt =     {'2016': 15., '2016preVFP': 15., '2017': 15., '2018': 15.,
                 '2022': 15., '2022EE': 15.}
minElectronPt = {'2016': 15., '2016preVFP': 15., '2017': 15., '2018': 15.,
                 '2022': 15., '2022EE': 15.}

if args.isData:
    if args.year == '2016preVFP':
        with open(os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/GoldenJSON/13TeV_UL2016_GoldenJSON.txt", 'r') as f:
            data_json = json.load(f)
    # goldenJSON ReReco 2022 --> https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json
    elif args.year == '2022' or args.year == '2022EE':
        with open(os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/GoldenJSON/13TeV_ReReco2022_GoldenJSON.txt", 'r') as f:
            data_json = json.load(f)
    else:
        with open(os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/GoldenJSON/13TeV_UL"+args.year+"_GoldenJSON.txt", 'r') as f:
            data_json = json.load(f)

    filtered_data = dict()
    for run in data_json.keys():
        filtered_data[int(run)] = list()
        for lumi in data_json[run]:
            filtered_data[int(run)].extend(range(lumi[0], lumi[1]+1))


#b-tagging working point
b_tagging_wpValues = {
    '2016preVFP': {'loose': 0.0614, 'medium': 0.3093, 'tight': 0.7221}, #https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016preVFP/
    '2016': {'loose': 0.0480, 'medium': 0.2489, 'tight': 0.6377}, #https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016postVFP/
    '2017': {'loose': 0.0532, 'medium': 0.3040, 'tight': 0.7476}, #https://btv-wiki.docs.cern.ch/ScaleFactors/UL2017/
    '2018': {'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100}, #https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    '2022': {'loose': 0.0583, 'medium': 0.3086, 'tight': 0.7183}, #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
    '2022EE': {'loose': 0.0614, 'medium': 0.3196, 'tight': 0.73}, #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
}

b_tagging_efficiencyMap = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btag_efficiency_maps/2016",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btag_efficiency_maps/2016preVFP",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btag_efficiency_maps/2017",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btag_efficiency_maps/2018",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btag_efficiency_maps/2022",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btag_efficiency_maps/2018",
}

b_tagging_SF =    os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/btagSF/"

# --- JERC webpage: https://cms-jerc.web.cern.ch/
jesUncertaintyFilesRegrouped = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_V2_MC_UncertaintySources_AK4PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_V2_MC_UncertaintySources_AK4PFPuppi.txt",
}
jerResolutionFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi.txt"
}
jerSFUncertaintyFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK4PFchs.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_JRV1_MC_SF_AK4PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_JRV1_MC_SF_AK4PFPuppi.txt",
}

jesAK8UncertaintyFilesRegrouped = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16_V7_MC_UncertaintySources_AK8PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16APV_V7_MC_UncertaintySources_AK8PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_V5_MC_UncertaintySources_AK8PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_V2_MC_UncertaintySources_AK8PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_V2_MC_UncertaintySources_AK8PFPuppi.txt",
}
jerAK8ResolutionFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_JRV1_MC_PtResolution_AK8PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK8PFPuppi.txt",
}
jerAK8SFUncertaintyFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK8PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK8PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_JRV1_MC_SF_AK8PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_JRV1_MC_SF_AK8PFPuppi.txt"
}

jesHOTVRUncertaintyFilesRegrouped = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16_V7_MC_UncertaintySources_AK4PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16APV_V7_MC_UncertaintySources_AK4PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_V5_MC_UncertaintySources_AK4PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_V5_MC_UncertaintySources_AK4PFPuppi.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_V2_MC_UncertaintySources_AK4PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_V2_MC_UncertaintySources_AK4PFPuppi.txt",
}
jerHOTVRResolutionFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK4PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK4PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK4PFPuppi.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi.txt"
}
jerHOTVRSFUncertaintyFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK4PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK4PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK4PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK4PFPuppi.txt",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22_22Sep2023_JRV1_MC_SF_AK4PFPuppi.txt",
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer22EE_22Sep2023_JRV1_MC_SF_AK4PFPuppi.txt",
}
# ---

#dilepton triggers SF provided by TOP PAG --> https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#Dilepton_triggers
dileptonTriggerSFFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2016postVFP_UL/TriggerSF_2016postVFP_ULv2.root",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2016preVFP_UL/TriggerSF_2016preVFP_ULv2.root",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2017_UL/TriggerSF_2017_ULv2.root",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2018_UL/TriggerSF_2018_ULv2.root",
    '2022':       "", #CHANGE!!!
    '2022EE':       "", #CHANGE!!!
}

#muon and electron SF files from the central POG repository--> https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG
muonSFFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2016postVFP_UL/muon_Z.json.gz",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2016preVFP_UL/muon_Z.json.gz",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2017_UL/muon_Z.json.gz", 
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2018_UL/muon_Z.json.gz", 
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2022/muon_Z.json.gz", 
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2022EE/muon_Z.json.gz", 
}

# EGM POG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/POG/EGM/
electronSFFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2016postVFP/electron.json.gz",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2016preVFP/electron.json.gz",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2017/electron.json.gz", 
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2018/electron.json.gz",
    '2022':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2022/electron.json.gz", 
    '2022EE':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2022EE/electron.json.gz", 
}

#xgb models
xgb_models = {
    '2022EE': os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/bdt/2022EE/hadronicTopVSQCD_bdt_nTrees1000_maxDepth5_learningRate0.05_minChildWeight0.5.bin",
    '2022': os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/bdt/2022/hadronicTopVSQCD_bdt_nTrees500_maxDepth5_learningRate0.08_minChildWeight0.5.bin",
    '2018': os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/bdt/2018/hadronicTopVSQCD_bdt_nTrees1000_maxDepth3_learningRate0.05_minChildWeight0.5.bin",
    '2017': os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/bdt/2017/hadronicTopVSQCD_bdt_nTrees500_maxDepth5_learningRate0.08_minChildWeight0.5.bin", 
    '2016': os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/bdt/2016/hadronicTopVSQCD_bdt_nTrees1000_maxDepth5_learningRate0.03_minChildWeight0.5.bin",
    '2016preVFP': os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/bdt/2016preVFP/hadronicTopVSQCD_bdt_nTrees500_maxDepth5_learningRate0.05_minChildWeight0.5.bin",
}

##### LEPTON MODULES
muon_dict, electron_dict = {}, {}
def leptonSequence():
    seq = [
        MuonSelection(
            inputCollection=lambda event: Collection(event, "Muon"),
            outputName_list=["tightRelIso_tightID_Muons", "tightRelIso_mediumID_Muons", "tightRelIso_looseID_Muons"],
            triggerMatch=True,
            storeKinematics=['pt', 'eta', 'charge', 'phi', 'mass', 'pfRelIso04_all', 'pfRelIso03_all'],
            muonMinPt=minMuonPt[args.year],
            muonMaxEta=2.4,
            storeTruthKeys = ['genPartIdx'],
        ),
               
        #MuonVeto(
        #    inputCollection=lambda event: event.unselectedMuons,
        #    outputName = "vetoMuons",
        #    muonMinPt = 10.,
        #    muonMaxEta = 2.4,
        #    storeKinematics=['pt','eta','charge','phi','mass'],
        #),

        ElectronSelection(
            inputCollection = lambda event: Collection(event, "Electron"),
            id_type = ['MVA', 'cutBased'],
            #outputName_list = ["tight_MVA_Electrons","medium_MVA_Electrons","loose_MVA_Electrons"],
            triggerMatch=True,
            electronMinPt = minElectronPt[args.year],
            electronMaxEta = 2.4,
            storeKinematics=['pt', 'eta', 'charge', 'phi', 'mass', 'dxy', 'dz', 'pfRelIso03_all', 'miniPFRelIso_all', 'jetRelIso'],
            storeTruthKeys = ['genPartIdx'], 
        ),

        #ElectronVeto(
        #    inputCollection=lambda event: event.unselectedElectrons,
        #    outputName = "vetoElectrons",
        #    electronMinPt = 10.,
        #    electronMaxEta = 2.4,
        #    storeKinematics=['pt','eta','charge','phi','mass'],
        #),

        # EventSkim(selection=lambda event: (event.nMuon + event.nElectron) > 1 ),
        # EventSkim(selection=lambda event: (event.nMuon + event.nloose_MVA_Electrons) > 1 ),
    ]

    if not Module.globalOptions["isData"]:
        for wp, muon_collection in muon_dict.items():
            seq.extend([
                MuonSFProducer(
                    muonSFFiles[args.year],
                    inputMuonCollection = muon_collection,
                    nosyst = False,
                    outputName = 'tightRelIso_'+wp+'ID_Muons',
                    WP = wp,
                    ),
            ])
        
        for wp in electron_dict.keys():
            for id_type, electron_collection in electron_dict[wp].items():
                seq.extend([
                    ElectronSFProducer(
                        electronSFFiles[args.year],
                        inputElectronCollection = electron_collection,
                        nosyst = False,
                        outputName = wp+'_'+id_type+'_Electrons',
                        WP = wp,
                        ID = id_type
                    )
                ])

    return seq

muon_dict = {
    'tight': lambda event: event.tightRelIso_tightID_Muons,
    'medium': lambda event: event.tightRelIso_mediumID_Muons,
    'loose': lambda event: event.tightRelIso_looseID_Muons
}

electron_dict = {
    'tight': 
        {"MVA": lambda event: event.tight_MVA_Electrons, "cutBased": lambda event: event.tight_cutBased_Electrons},
    'medium':
        {"MVA": lambda event: event.medium_MVA_Electrons, "cutBased": lambda event: event.medium_cutBased_Electrons},
    'loose': 
        {"MVA": lambda event: event.loose_MVA_Electrons, "cutBased": lambda event: event.loose_cutBased_Electrons}
}

muon_collection_for_selection_and_cleaning = 'tightRelIso_looseID_Muons'
electron_collection_for_selection_and_cleaning = 'medium_MVA_Electrons'
if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
    electron_collection_for_selection_and_cleaning = 'medium_MVA_Electrons'
    muon_collection_for_selection_and_cleaning = 'tightRelIso_mediumID_Muons'

#####

##### TRIGGER MODULES - DILEPTON TRIGGER
def trigger():
    store_weights_trigger = False
    if not Module.globalOptions["isData"]:
        store_weights_trigger = True
    seq = [
        TriggerSelection(
            dileptonTriggerSFFiles[args.year],
            inputMuonCollection = lambda event: Collection(event, "Muon"),
            inputElectronCollection = lambda event: Collection(event, "Electron"),
            outputName="trigger",
            storeWeights=store_weights_trigger,
            thresholdPt=15. 
        ),
        # EventSkim(selection=lambda event: (event.trigger_general_flag)),
    ]
    return seq
#####

##### JET MODULES   
selectedFatJets_dict, selectedJets_dict, selectedBJets_dict = {}, {}, {}
def jetSelection(jetDict):
    seq = []
    
    for systName, (jetCollection, fatjetCollection, hotvrjetCollection, subhotvrjetCollection) in jetDict.items():
        seq.extend([
            JetSelection(
                inputCollection= jetCollection, 
                leptonCollectionDRCleaning=lambda event: getattr(event, muon_collection_for_selection_and_cleaning) + getattr(event, electron_collection_for_selection_and_cleaning), 
                jetMinPt=30.,
                jetMaxEta=2.4,
                dRCleaning=0.4,
                jetId=JetSelection.LOOSE,
                storeKinematics=['pt', 'eta','phi','mass','btagDeepFlavB', 'area',
                                 'minDPhiClean', 'minDRClean', 'jetId'],
                outputName_list=["selectedJets_"+systName, "unselectedJets_"+systName],
                metInput = lambda event: Object(event, "MET"),
                storeTruthKeys = ['hadronFlavour','partonFlavour'],
            ),
            #TODO: every ak8 will also be ak4 -> some cross cleaning required
            JetSelection(
                inputCollection= fatjetCollection, 
                leptonCollectionDRCleaning=lambda event,sys=systName: getattr(event, muon_collection_for_selection_and_cleaning) + getattr(event, electron_collection_for_selection_and_cleaning),
                jetMinPt=400., 
                jetMaxEta=2.4,
                dRCleaning=0.8,
                jetId=JetSelection.LOOSE,
                storeKinematics=['pt', 'eta','phi','mass', 'genJetAK8Idx', 'deepTag_TvsQCD', 'deepTag_WvsQCD',
                                  'particleNet_TvsQCD', 'particleNet_WvsQCD', 'particleNet_QCD',
                                    'particleNet_mass', 'btagDeepB', 'tau2', 'tau3', 'tau1', 'msoftdrop', 'area',
                                    'minDPhiClean', 'minDRClean', 'jetId'],
                outputName_list=["selectedFatJets_"+systName,"unselectedFatJets_"+systName],
                metInput = lambda event: Object(event, "MET"),
                storeTruthKeys = ['hadronFlavour','nBHadrons', 'nCHadrons', 'genJetAK8Idx'],
            ),
            JetSelection(
                inputCollection= hotvrjetCollection, 
                leptonCollectionDRCleaning=lambda event: getattr(event, muon_collection_for_selection_and_cleaning) + getattr(event, electron_collection_for_selection_and_cleaning),
                jetMinPt=200., 
                jetMaxEta= 2.4,
                dRCleaning=None,
                jetId=JetSelection.NONE,
                storeKinematics=['pt','eta','phi','mass', 
                                 'tau2', 'tau3', 'tau1', 'area', 
                                 'nConstituents', 'subJetIdx1', 'subJetIdx2', 'subJetIdx3',
                                 'minDPhiClean', 'minDRClean'],
                outputName_list=["selectedHOTVRJets_"+systName, "unselectedHOTVRJets_"+systName],
                metInput = lambda event: Object(event, "MET"),
                # storeTruthKeys = ['hadronFlavour','partonFlavour'],
                )
        ])
        # if systName == 'nominal':
        seq.extend([
            JetSelection(
                inputCollection= subhotvrjetCollection,
                leptonCollectionDRCleaning=lambda event: getattr(event, muon_collection_for_selection_and_cleaning) + getattr(event, electron_collection_for_selection_and_cleaning),
                jetMinPt=30., 
                jetMaxEta= 2.4,
                dRCleaning=None,
                jetId=JetSelection.NONE,
                storeKinematics=['pt', 'eta', 'phi', 'mass', '_index', 'area', 'minDPhiClean', 'minDRClean'],
                outputName_list=["selectedHOTVRSubJets_"+systName, "unselectedHOTVRSubJets_"+systName],
                metInput = lambda event: Object(event, "MET"),
                # storeTruthKeys = ['hadronFlavour','partonFlavour'],
        )])


        seq.append(
            BTagSelection(
                btaggingWP=b_tagging_wpValues[args.year],
                inputCollection=lambda event,sys=systName: getattr(event,"selectedJets_"+sys),
                outputName='selectedBJets_'+systName,
                jetMinPt=30.,
                jetMaxEta=2.4,
                storeKinematics=['pt', 'eta','phi','mass', 'area', 'minDPhiClean', 'minDRClean'],
                storeTruthKeys = ['hadronFlavour', 'partonFlavour', 'genJetIdx', 'jetId'],
            )
        )


    seq.extend([
        JetSelection(
            inputCollection= lambda event: Collection(event,"HOTVRJet"), 
            leptonCollectionDRCleaning=lambda event: getattr(event, muon_collection_for_selection_and_cleaning) + getattr(event, electron_collection_for_selection_and_cleaning),
            jetMinPt=200., 
            jetMaxEta= 2.4,
            dRCleaning=None,
            jetId=JetSelection.NONE,
            storeKinematics=['pt','eta','phi','mass', 'btagDeepB', 'tau2', 'tau3', 'tau1', 'area', 
                             'nConstituents', 'subJetIdx1', 'subJetIdx2', 'subJetIdx3',
                             'minDPhiClean', 'minDRClean'],
            outputName_list=["preselectedHOTVRJets","preunselectedHOTVRJets"],
            metInput = lambda event: Object(event, "MET"),
            # storeTruthKeys = ['hadronFlavour','partonFlavour'],
            )])

    systNames = jetDict.keys()
   
    # at least 1 HOTVR
    # seq.append(
    #    EventSkim(selection=lambda event: 
    #        getattr(event, "npreselectedHOTVRJets") > 0
    #    ),
    # )
   
    #at least 2 AK8 jets
    #seq.append(
    #    EventSkim(selection=lambda event, systNames=systNames: 
    #        any([getattr(event, "nselectedFatJets_"+systName) >= 2 for systName in systNames])
    #    )
    #)
    
    # fixed working point event reweighting and shape implemented
    # https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#b_tagging_efficiency_in_MC_sampl
    if isMC:
        jesUncertForBtag = ['jes'+syst.replace('Total','') for syst in jesUncertaintyNames]
        seq.append(
            btagSFProducer(
                era=args.year,
                # algo = 'deepJet',
                sfFileName = b_tagging_SF,
                efficiencyMaps = b_tagging_efficiencyMap[str(args.year)],
                jesSystsForShape = jesUncertForBtag,
                nosyst = args.nosys,
                selectedWPs = ['L', 'M', 'T']
            )
        )

    return seq

# selectedFatJets_dict = {systName: lambda event,sys=systName: getattr(event,"selectedFatJets_"+sys)}
# selectedJets_dict = {systName: lambda event,sys=systName: getattr(event,"selectedJets_"+sys)}
# selectedBJets_dict = {systName: [lambda event,sys=systName: getattr(event,"selectedBJets_"+sys+"_tight"),lambda event,sys=systName: getattr(event,"selectedBJets_"+sys+"_medium"),lambda event,sys=systName: getattr(event,"selectedBJets_"+sys+"_loose")]}
#####

##### EVENT INFO MODULE
# storeVariables = []
storeVariables = [[lambda tree: tree.branch("run", "L"),
        lambda tree, event: tree.fillBranch("run", event.run)],
        [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
        event: tree.fillBranch("PV_npvs", event.PV_npvs)],
        [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
        event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
        [lambda tree: tree.branch("luminosityBlock", "I"), lambda tree,
        event: tree.fillBranch("luminosityBlock", event.luminosityBlock)],
        [lambda tree: tree.branch("eventID", "L"), lambda tree,
        event: tree.fillBranch("eventID", event.event)]
    ]
if not Module.globalOptions["isData"]:
    storeVariables.extend([
        [lambda tree: tree.branch("genweight", "F"), 
         lambda tree, event: tree.fillBranch("genweight", event.Generator_weight)],
        [lambda tree: tree.branch("LHE_HT", "F"), 
         lambda tree, event: tree.fillBranch("LHE_HT", event.LHE_HT)]
        ])
    if Module.globalOptions["year"] in ['2016', '2016preVFP', '2017']:
        for l1prefiringw in ['Nom', 'Up', 'Dn']:
            storeVariables.append([
                lambda tree, l1prefiringw=l1prefiringw: tree.branch("L1PreFiringWeight_{}".format(l1prefiringw), "F"),
                lambda tree, event, l1prefiringw=l1prefiringw: tree.fillBranch("L1PreFiringWeight_{}".format(l1prefiringw), getattr(event, "L1PreFiringWeight_{}".format(l1prefiringw)))])

if not Module.globalOptions["year"] == '2022' and not Module.globalOptions["year"] == '2022EE':
    storeVariables.append([
        lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"),
          lambda tree, event: tree.fillBranch("fixedGridRhoFastjetAll", event.fixedGridRhoFastjetAll)])
else:
    storeVariables.append([
        lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"),
          lambda tree, event: tree.fillBranch("fixedGridRhoFastjetAll", event.Rho_fixedGridRhoFastjetAll)])

#####ANALYZER CHAIN : DATA ARE FILTERED FROM GOLDENJSON RUNS
if not Module.globalOptions["isData"]:  
    analyzerChain = [
        EventInfo(storeVariables=storeVariables),
        # EventSkim(selection=lambda event: event.nTrigObj > 0),
        MetFilter(
            globalOptions=globalOptions,
            outputName="MET_filter"
            ),
        ]
else:
    analyzerChain = [
        EventInfo(storeVariables=storeVariables),
        EventSkim(selection=lambda event: event.run in filtered_data.keys()),
        EventSkim(selection=lambda event: event.luminosityBlock in filtered_data[event.run]) ,
        # EventSkim(selection=lambda event: event.nTrigObj > 0),
            MetFilter(
            globalOptions=globalOptions,
            outputName="MET_filter"
        ),
]
#####

analyzerChain.extend(trigger())
analyzerChain.extend(leptonSequence())

#####JETMET UNCERTAINTIES MODULE
if args.isData:
    analyzerChain.extend(
        jetSelection({
            "nominal": (lambda event: Collection(event,"Jet"),
                        lambda event: Collection(event,"FatJet"),
                        lambda event: Collection(event,"HOTVRJet"),
                         lambda event: Collection(event,"HOTVRSubJet"))
        })
    )

    jetDict = {
        "nominal": (lambda event: event.selectedJets_nominal,
                    lambda event: event.selectedFatJets_nominal,
                    lambda event: event.selectedHOTVRJets_nominal,
                    lambda event: event.selectedHOTVRSubJets_nominal)
    }

else:
    # pu weight producer are set to AUTO mode for 2022, 2022EE --> https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#True_and_observed_pileup_in_data
    # 1) no MC distributions available (?) - 07/03/2024
    # 2) no plus/minus variation
    analyzerChain.append(PUWeightProducer_dict[args.year]())

    if args.nosys:
        jesUncertaintyNames = []
    else:
        
        jesUncertaintyNames = ["Total","Absolute","EC2","BBEC1", "HF","RelativeBal","FlavorQCD" ]
        for jesUncertaintyExtra in ["RelativeSample","HF","Absolute","EC2","BBEC1"]:
            jesUncertaintyNames.append(jesUncertaintyExtra+"_"+args.year.replace("preVFP",""))
        
        jesUncertaintyNames = ["Total"]
            
        print "JECs: ",jesUncertaintyNames
        
    #TODO: apply type2 corrections? -> improves met modelling; in particular for 2018
    if Module.globalOptions["year"] == '2022' or Module.globalOptions["year"] == '2022EE':
        rhoInput = lambda event: event.Rho_fixedGridRhoFastjetAll
    else: 
        rhoInput = lambda event: event.fixedGridRhoFastjetAll

    analyzerChain.extend([
        JetMetUncertainties(
            jesUncertaintyFilesRegrouped[args.year],
            jerResolutionFiles[args.year],
            jerSFUncertaintyFiles[args.year],
            jesUncertaintyNames = jesUncertaintyNames, 
            metInput = lambda event: Object(event, "MET"),
            rhoInput = rhoInput,
            jetCollection = lambda event: Collection(event,"Jet"),
            lowPtJetCollection = lambda event: Collection(event,"CorrT1METJet"),
            genJetCollection = lambda event: Collection(event,"GenJet"),
            muonCollection = lambda event: Collection(event,"Muon"),
            electronCollection = lambda event: Collection(event,"Electron"),
            propagateJER = False, #not recommended
            outputJetPrefix = 'jets_',
            outputMetPrefix = 'met_',
            jetKeys=['jetId', 'nConstituents','btagDeepFlavB','hadronFlavour','partonFlavour','genJetIdx', 'area'],
        ),
        JetMetUncertainties(
            jesAK8UncertaintyFilesRegrouped[args.year],
            jerAK8ResolutionFiles[args.year],
            jerAK8SFUncertaintyFiles[args.year],
            jesUncertaintyNames = jesUncertaintyNames, 
            metInput = None,
            rhoInput = rhoInput,
            jetCollection = lambda event: Collection(event,"FatJet"),
            lowPtJetCollection = lambda event: [],
            genJetCollection = lambda event: Collection(event,"GenJetAK8"),
            muonCollection = lambda event: Collection(event,"Muon"),
            electronCollection = lambda event: Collection(event,"Electron"),
            propagateToMet = False, #no met variations
            propagateJER = False, #not recommended
            outputJetPrefix = 'fatjets_',
            outputMetPrefix = None,
            jetKeys=['jetId', 'btagDeepB','deepTag_TvsQCD','deepTag_WvsQCD',
                     'particleNet_TvsQCD','particleNet_WvsQCD','particleNet_mass', 
                     'particleNet_QCD','hadronFlavour', 'genJetAK8Idx', 'nBHadrons', 
                     'tau2', 'tau3', 'tau1', 'msoftdrop', 'nCHadrons', 'nBHadrons', 'area'],  #'nConstituents'
        ), 
        JetHOTVRUncertainties(
            jesHOTVRUncertaintyFilesRegrouped[args.year],
            jerHOTVRResolutionFiles[args.year],
            jerHOTVRSFUncertaintyFiles[args.year],
            jesUncertaintyNames = jesUncertaintyNames, 
            rhoInput = rhoInput,
            jetCollection = lambda event: Collection(event,"HOTVRJet"),
            subjetCollection = lambda event: Collection(event,"HOTVRSubJet"),
            lowPtJetCollection = lambda event: [],
            genJetCollection = lambda event: Collection(event,"GenHOTVRJet"),
            genSubJetCollection = lambda event: Collection(event,"GenHOTVRSubJet"),
            outputJetPrefix = 'hotvrjets_',
            outputSubJetPrefix = 'hotvrsubjets_',
            jetKeys=['nConstituents', 'area', 'tau2', 'tau3', 'tau1', 'subJetIdx1', 'subJetIdx2', 'subJetIdx3'],
        )
    ])

    jetDict = {
        "nominal": (lambda event: event.jets_nominal,
                    lambda event: event.fatjets_nominal,
                    lambda event: event.hotvrjets_nominal,
                    lambda event: event.hotvrsubjets_nominal)
    }
    
    if not args.nosys:
        jetDict["jerUp"] = (
            lambda event: event.jets_jerUp,
            lambda event: event.fatjets_jerUp,
            lambda event: event.hotvrjets_jerUp,
             lambda event: event.hotvrsubjets_jerUp,)
        jetDict["jerDown"] = (
            lambda event: event.jets_jerDown,
            lambda event: event.fatjets_jerDown,
            lambda event: event.hotvrjets_jerDown,
            lambda event: event.hotvrsubjets_jerDown)
        
        for jesUncertaintyName in jesUncertaintyNames:
            jetDict['jes'+jesUncertaintyName+"Up"] = (
                lambda event,sys=jesUncertaintyName: getattr(event,"jets_jes"+sys+"Up"),
                lambda event,sys=jesUncertaintyName: getattr(event,"fatjets_jes"+sys+"Up"),
                lambda event,sys=jesUncertaintyName: getattr(event,"hotvrjets_jes"+sys+"Up"),
                lambda event,sys=jesUncertaintyName: getattr(event,"hotvrsubjets_jes"+sys+"Up"),)
            jetDict['jes'+jesUncertaintyName+"Down"] = (
                lambda event,sys=jesUncertaintyName: getattr(event,"jets_jes"+sys+"Down"),
                lambda event,sys=jesUncertaintyName: getattr(event,"fatjets_jes"+sys+"Down"),
                lambda event,sys=jesUncertaintyName: getattr(event,"hotvrjets_jes"+sys+"Down"),
                lambda event,sys=jesUncertaintyName: getattr(event,"hotvrsubjets_jes"+sys+"Down"),)
    
    analyzerChain.extend(
        jetSelection(jetDict)
    )
#####


##### GENERATION MODULE
if isMC:
    # for systName,(jetCollection,fatjetCollection) in jetDict.items():
        analyzerChain.extend( [
            GenParticleModule(
                inputGenCollection=lambda event: Collection(event, "GenPart"),
                inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
                inputGenJetCollection=lambda event: Collection(event, "GenJet"),
                inputMuonCollection=lambda event: getattr(event, muon_collection_for_selection_and_cleaning),
                inputElectronCollection=lambda event: getattr(event, electron_collection_for_selection_and_cleaning),
                outputName="genPart",
                storeKinematics= ['pt','eta','phi','mass'],
            ),
        ])
####


##### EVENT RECONSTRUCTION MODULE
triggers = {'ee': lambda event: event.trigger_ee_flag, 'emu': lambda event: event.trigger_emu_flag, 'mumu': lambda event: event.trigger_mumu_flag}

for systName in jetDict.keys():
    event_reco_inputs = {
        'inputTriggersCollection': triggers,
        'inputMuonCollection': lambda event: getattr(event, muon_collection_for_selection_and_cleaning),
        'inputElectronCollection': lambda event: getattr(event, electron_collection_for_selection_and_cleaning),
        'inputJetCollection': lambda event: getattr(event, "selectedJets_{}".format(systName)),
        'inputBJetCollection': lambda event: event.selectedBJets_nominal_loose,
        'inputFatJetCollection': lambda event: getattr(event, "selectedFatJets_{}".format(systName)),
        'inputMETCollection': [],  # to be included!
        'inputHOTVRJetCollection': lambda event: getattr(event, "selectedHOTVRJets_{}".format(systName)),
        'inputHOTVRSubJetCollection': lambda event: getattr(event, "selectedHOTVRSubJets_{}".format(systName)),
        'systName': systName
    }
    if not Module.globalOptions["isData"]: 
        event_reco_inputs['inputGenTopCollection'] = lambda event: event.genTops

    analyzerChain.extend([
        EventReconstruction(**event_reco_inputs),
        # EventSkim(selection=lambda event: (event.event_selection_OS_dilepton_cut))
    ])

    # #### XGB EVALUATION MODULE
    analyzerChain.append(
        XGBEvaluationProducer(
            modelPath=xgb_models[args.year],
            inputHOTVRJetCollection=lambda event: getattr(event, "selectedHOTVRJets_{}".format(systName)),
            outputName="scoreBDT",
            outputJetPrefix='selectedHOTVRJets_'+systName
        )
    )

##### HOTVR/AK8 JET COMPOSITION MODULE 
if not Module.globalOptions["isData"]:
    # analyzerChain.append(
    #     HOTVR_MVA(
    #         inputHOTVRJetCollection = lambda event: getattr(event,"selectedHOTVRJets_nominal"),
    #         inputGenParticleCollections = {
    #             'gentops': lambda event: event.genTops, 
    #             'genWs_not_from_top': lambda event: event.gen_w_bosons_not_from_top, 
    #             'genbs_not_from_top': lambda event: event.gen_b_quarks_not_from_top, 
    #             'genparticles_not_from_top': lambda event: event.gen_particles_not_from_top
    #             },
    #         inputSubHOTVRJetCollection = lambda event: getattr(event,"selectedHOTVRSubJets_nominal"),
    #     )
    # )
    analyzerChain.append(
        HOTVRJetComposition(
            inputHOTVRJetCollection = lambda event: getattr(event,"selectedHOTVRJets_nominal"),
            inputGenParticleCollections = {
                'gentops': lambda event: event.genTops, 
                'genWs_not_from_top': lambda event: event.gen_w_bosons_not_from_top, 
                'genbs_not_from_top': lambda event: event.gen_b_quarks_not_from_top, 
                'genparticles_not_from_top': lambda event: event.gen_particles_not_from_top
                },
            inputSubHOTVRJetCollection = lambda event: getattr(event,"selectedHOTVRSubJets_nominal"),
            outputJetPrefix='selectedHOTVRJets_nominal'
        )
    )
    analyzerChain.append(
        FatJetComposition(
            inputFatJetCollection = lambda event: getattr(event,"selectedFatJets_nominal"),
            inputGenParticleCollections = {
                'gentops': lambda event: event.genTops, 
                'genWs_not_from_top': lambda event: event.gen_w_bosons_not_from_top, 
                'genbs_not_from_top': lambda event: event.gen_b_quarks_not_from_top, 
                'genparticles_not_from_top': lambda event: event.gen_particles_not_from_top
                },
            outputJetPrefix='selectedFatJets_nominal'
        )
    )


if not args.isData:
    #analyzerChain.append(GenWeightProducer())
    if isPowhegTTbar:
        analyzerChain.append(
            TopPtWeightProducer(
                mode=TopPtWeightProducer.DATA_NLO
            )
        )

p = PostProcessor(
    args.output[0],
    args.inputFiles,
    cut="",#"(nJet>1)&&((nElectron+nMuon)>1)", #at least 2 jets + 2 leptons
    modules=analyzerChain,
    friend=True,
    maxEntries = args.maxEvents,
    # noOut = True
)

p.run()

