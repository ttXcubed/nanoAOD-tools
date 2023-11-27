import os
import sys
import math
import argparse
import random
import ROOT
import numpy as np
import json


from collections import OrderedDict
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor \
    import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *

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
                    action='store', type=str, default='2017', choices=['2016','2016preVFP','2017','2018'])
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
    "era": None
}

Module.globalOptions = globalOptions

isMC = not args.isData
isPowheg = 'powheg' in args.inputFiles[0].lower()
isPowhegTTbar = 'TTTo' in args.inputFiles[0] and isPowheg

#recommended pT threshold for the subleading lepton --> https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2020/085 (4 top in dilepton final state)
minMuonPt =     {'2016': 15., '2016preVFP': 15., '2017': 15., '2018': 15.}
minElectronPt = {'2016': 15., '2016preVFP': 15., '2017': 15., '2018': 15.}

if args.isData:
    with open(os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/GoldenJSON/13TeV_UL"+args.year+"_GoldenJSON.txt", 'r') as f:
        data_json = json.load(f)

    filtered_data=dict()
    for run in data_json.keys():
        filtered_data[int(run)] = list()
        for lumi in data_json[run]:
            filtered_data[int(run)].extend(range(lumi[0],lumi[1]+1))



#b-tagging working point
b_tagging_wpValues = {
    '2016preVFP': [0.0614, 0.3093, 0.7221], #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    '2016': [0.0480, 0.2489, 0.6377], #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16postVFP
    '2017': [0.0532, 0.3040, 0.7476], #https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    '2018': [0.0490, 0.2783, 0.7100] #https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
}

jesUncertaintyFilesRegrouped = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt"
}
jerResolutionFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt"
}
jerSFUncertaintyFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK4PFchs.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt"
}

jesAK8UncertaintyFilesRegrouped = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16_V7_MC_UncertaintySources_AK8PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16APV_V7_MC_UncertaintySources_AK8PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_V5_MC_UncertaintySources_AK8PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.txt"
}
jerAK8ResolutionFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt"
}
jerAK8SFUncertaintyFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK8PFPuppi.txt",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.txt",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK8PFPuppi.txt",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt"
}

#dilepton triggers SF provided by TOP PAG --> https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#Dilepton_triggers
dileptonTriggerSFFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2016postVFP_UL/TriggerSF_2016postVFP_ULv2.root",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2016preVFP_UL/TriggerSF_2016preVFP_ULv2.root",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2017_UL/TriggerSF_2017_ULv2.root",
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/dilepton_triggerSF/2018_UL/TriggerSF_2018_ULv2.root"    
}

#muon and electron SF files from the central POG repository--> https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG
muonSFFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2016postVFP_UL/muon_Z.json.gz",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2016preVFP_UL/muon_Z.json.gz",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2017_UL/muon_Z.json.gz", 
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/muon/2018_UL/muon_Z.json.gz"      
}

electronSFFiles = {
    '2016':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2016postVFP/electron.json.gz",
    '2016preVFP': os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2016preVFP/electron.json.gz",
    '2017':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2017/electron.json.gz", 
    '2018':       os.environ['CMSSW_BASE']+"/src/PhysicsTools/NanoAODTools/data/electron/2018/electron.json.gz"      
}

##### LEPTON MODULES
def leptonSequence():
    seq = [
        MuonSelection(
            inputCollection=lambda event: Collection(event, "Muon"),
            outputName_list=["tightRelIso_tightID_Muons","tightRelIso_mediumID_Muons","tightRelIso_looseID_Muons"],
            triggerMatch=True,
            storeKinematics=['pt','eta','charge','phi','mass'],
            muonMinPt=minMuonPt[args.year],
            muonMaxEta=2.4,
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
            storeKinematics=['pt','eta','charge','phi','mass'],#, 
        ),

        #ElectronVeto(
        #    inputCollection=lambda event: event.unselectedElectrons,
        #    outputName = "vetoElectrons",
        #    electronMinPt = 10.,
        #    electronMaxEta = 2.4,
        #    storeKinematics=['pt','eta','charge','phi','mass'],
        #),

        EventSkim(selection=lambda event: (len(event.tightRelIso_looseID_Muons) + ( len(event.loose_cutBased_Electrons) + len(event.loose_MVA_Electrons) ) > 1 )),
    ]

    if not Module.globalOptions["isData"]:
        seq.extend([
            MuonSFProducer(
            muonSFFiles[args.year],
            inputMuonCollection = OrderedDict([('tight',lambda event: event.tightRelIso_tightID_Muons), ('medium',lambda event: event.tightRelIso_mediumID_Muons), ('loose',lambda event: event.tightRelIso_looseID_Muons)]),
            nosyst = False
            ),
            ElectronSFProducer(
            electronSFFiles[args.year],
            inputElectronCollection = OrderedDict([ ("MVA", OrderedDict([('tight',lambda event: event.tight_MVA_Electrons), ('medium',lambda event: event.medium_MVA_Electrons), ('loose',lambda event: event.loose_MVA_Electrons)]) ), \
               ("cutBased", OrderedDict([('tight',lambda event: event.tight_cutBased_Electrons), ('medium',lambda event: event.medium_cutBased_Electrons), ('loose',lambda event: event.loose_cutBased_Electrons)])) ]),
            nosyst = False
            ),
        ])

    return seq
#####

##### TRIGGER MODULES
def trigger():
	
    if args.isData:
        if args.trigger=='emu':
            seq = [
                ElectronMuonTriggerSelection(
                    outputName="trigger",
                ),
                EventSkim(selection=lambda event: (event.trigger_emu_flag)),
            ]
            return seq
        elif args.trigger=='mumu':
            seq = [
                DoubleMuonTriggerSelection(
                    outputName="trigger",
                ),	
                EventSkim(selection=lambda event: (event.trigger_mumu_flag)),
            ]
            return seq
        elif args.trigger=='ee':
            seq = [
                DoubleElectronTriggerSelection(
                    outputName="trigger",
                ),
                EventSkim(selection=lambda event: (event.trigger_ee_flag)),
            ] 
            return seq
    else: 
        seq = [
            DoubleLeptonTriggerSelection(
                dileptonTriggerSFFiles[args.year],
                inputMuonCollection = lambda event: event.tightRelIso_looseID_Muons,
                inputElectronCollection = lambda event: event.loose_MVA_Electrons,
                outputName="trigger",
                storeWeights=True,
                thresholdPt=15. #minDiMuonPt[args.year], PRELIMINARY VALUE
            ),	
            EventSkim(selection=lambda event: (event.trigger_general_flag)),
        ]
        return seq
#####

##### JET MODULES   
def jetSelection(jetDict):
    seq = []
    
    for systName,(jetCollection,fatjetCollection) in jetDict.items():
        seq.extend([
            JetSelection(
                inputCollection= jetCollection, 
                leptonCollectionDRCleaning=lambda event: event.tightRelIso_looseID_Muons+event.loose_MVA_Electrons,#event.tightRelIso_tightID_Muons+event.tight_MVA_Electrons,
                jetMinPt=30.,
                jetMaxEta=2.4,
                dRCleaning=0.4,
                jetId=JetSelection.LOOSE,
                storeKinematics=['pt', 'eta','phi','mass'],
                outputName_list=["selectedJets_"+systName,"unselectedJets_"+systName],
                fatFlag=False,
                metInput = lambda event: Object(event, "MET"),
            ),
            #TODO: every ak8 will also be ak4 -> some cross cleaning required
            JetSelection(
                inputCollection= fatjetCollection, 
                leptonCollectionDRCleaning=lambda event,sys=systName: event.tightRelIso_looseID_Muons+event.loose_MVA_Electrons,
                jetMinPt=400., 
                jetMaxEta=2.4,
                dRCleaning=0.8,
                jetId=JetSelection.LOOSE,
                storeKinematics=['pt', 'eta','phi','mass'],
                outputName_list=["selectedFatJets_"+systName,"unselectedFatJets_"+systName],
        fatFlag=True,
        metInput = lambda event: Object(event, "MET"),
            )
        ])
        
        
        seq.append(
            BTagSelection(
                b_tagging_wpValues[args.year],
                inputCollection=lambda event,sys=systName: getattr(event,"selectedJets_"+sys),
                outputName_list=["selectedBJets_"+systName+"_tight", "selectedBJets_"+systName+"_medium", "selectedBJets_"+systName+"_loose"],
                jetMinPt=30.,
                jetMaxEta=2.4,
                workingpoint = [],
                storeKinematics=['pt', 'eta','phi','mass'],
                storeTruthKeys = ['hadronFlavour','partonFlavour'],
            )
        )
        
    systNames = jetDict.keys()
   
    #at least 4 AK4 jets and 2 AK8 jets
    #seq.append(
    #    EventSkim(selection=lambda event, systNames=systNames: 
    #        any([getattr(event, "nselectedJets_"+systName) >= 2 for systName in systNames])
    #    ),
    #)
   
    
    #seq.append(
    #    EventSkim(selection=lambda event, systNames=systNames: 
    #        any([len(filter(lambda jet: jet.isBTagged,getattr(event,"selectedJets_"+systName))) >= 2 for systName in systNames])
    #   )
    #)

    #seq.append(
    #    EventSkim(selection=lambda event, systNames=systNames: 
    #        any([getattr(event, "nselectedBJets_"+systName+"_loose") >= 2 for systName in systNames])
    #    ),
    #)
	   
    #at least 2 AK8 jets
    #seq.append(
    #    EventSkim(selection=lambda event, systNames=systNames: 
    #        any([getattr(event, "nselectedFatJets_"+systName) >= 2 for systName in systNames])
    #    )
    #)
    
    if isMC:
        jesUncertForBtag = ['jes'+syst.replace('Total','') for syst in jesUncertaintyNames]
        seq.append(
            btagSFProducer(
            era=args.year,
            jesSystsForShape = jesUncertForBtag,
            nosyst = args.nosys
            )
        )
        
    return seq
#####

##### EVENT INFO MODULE
if not Module.globalOptions["isData"]:
    storeVariables = [[lambda tree: tree.branch("genweight", "F"),
                        lambda tree,
                        event: tree.fillBranch("genweight",
                        event.Generator_weight)],
        [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
        event: tree.fillBranch("PV_npvs", event.PV_npvs)],
        [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
        event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
        [lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,
        event: tree.fillBranch("fixedGridRhoFastjetAll",
        event.fixedGridRhoFastjetAll)],
    ]
	
else: 
    storeVariables = [[lambda tree: tree.branch("run", "I"),
                            lambda tree,
                            event: tree.fillBranch("run",
                            event.run)],
        [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
        event: tree.fillBranch("PV_npvs", event.PV_npvs)],
        [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
        event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
        [lambda tree: tree.branch("luminosityBlock", "I"), lambda tree,
        event: tree.fillBranch("luminosityBlock", event.luminosityBlock)],
        [lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,
        event: tree.fillBranch("fixedGridRhoFastjetAll",
        event.fixedGridRhoFastjetAll)],
    ]
	
if not Module.globalOptions["isData"]:	
    analyzerChain = [EventInfo(storeVariables=storeVariables),
    EventSkim(selection=lambda event: event.nTrigObj > 0),
    MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
        ),
    ]
else:                   
    analyzerChain = [EventSkim(selection=lambda event: event.run in filtered_data.keys()),
    EventSkim(selection=lambda event: event.luminosityBlock in filtered_data[event.run]) ,
    EventInfo(storeVariables=storeVariables),
    EventSkim(selection=lambda event: event.nTrigObj > 0),
        MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
        ),
]
#####

analyzerChain.extend(leptonSequence())
analyzerChain.extend(trigger())

#####JETMET UNCERTAINTIES MODULE
if args.isData:
    analyzerChain.extend(
        jetSelection({
            "nominal": (lambda event: Collection(event,"Jet"),lambda event: Collection(event,"FatJet"))
        })
    )

else:
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
    analyzerChain.extend([
        JetMetUncertainties(
            jesUncertaintyFilesRegrouped[args.year],
            jerResolutionFiles[args.year],
            jerSFUncertaintyFiles[args.year],
            jesUncertaintyNames = jesUncertaintyNames, 
            metInput = lambda event: Object(event, "MET"),
            rhoInput = lambda event: event.fixedGridRhoFastjetAll,
            jetCollection = lambda event: Collection(event,"Jet"),
            lowPtJetCollection = lambda event: Collection(event,"CorrT1METJet"),
            genJetCollection = lambda event: Collection(event,"GenJet"),
            muonCollection = lambda event: Collection(event,"Muon"),
            electronCollection = lambda event: Collection(event,"Electron"),
            propagateJER = False, #not recommended
            outputJetPrefix = 'jets_',
            outputMetPrefix = 'met_',
            jetKeys=['jetId', 'nConstituents','btagDeepFlavB','hadronFlavour','partonFlavour','genJetIdx'],
        ),
        JetMetUncertainties(
            jesAK8UncertaintyFilesRegrouped[args.year],
            jerAK8ResolutionFiles[args.year],
            jerAK8SFUncertaintyFiles[args.year],
            jesUncertaintyNames = jesUncertaintyNames, 
            metInput = None,
            rhoInput = lambda event: event.fixedGridRhoFastjetAll,
            jetCollection = lambda event: Collection(event,"FatJet"),
            lowPtJetCollection = lambda event: [],
            genJetCollection = lambda event: Collection(event,"GenJetAK8"),
            muonCollection = lambda event: Collection(event,"Muon"),
            electronCollection = lambda event: Collection(event,"Electron"),
            propagateToMet = False, #no met variations
            propagateJER = False, #not recommended
            outputJetPrefix = 'fatjets_',
            outputMetPrefix = None,
            jetKeys=['jetId', 'btagDeepB','deepTag_TvsQCD','deepTag_WvsQCD','particleNet_TvsQCD','particleNet_WvsQCD','particleNet_mass', 'particleNet_QCD','hadronFlavour','genJetAK8Idx', 'nBHadrons', 'tau2', 'tau3'],  #'nConstituents'
        )
    ])

    jetDict = {
        "nominal": (lambda event: event.jets_nominal,lambda event: event.fatjets_nominal)
    }
    
    if not args.nosys:
        jetDict["jerUp"] = (lambda event: event.jets_jerUp,lambda event: event.fatjets_jerUp)
        jetDict["jerDown"] = (lambda event: event.jets_jerDown,lambda event: event.fatjets_jerDown)
        
        for jesUncertaintyName in jesUncertaintyNames:
            jetDict['jes'+jesUncertaintyName+"Up"] = (lambda event,sys=jesUncertaintyName: getattr(event,"jets_jes"+sys+"Up"),lambda event,sys=jesUncertaintyName: getattr(event,"fatjets_jes"+sys+"Up"))
            jetDict['jes'+jesUncertaintyName+"Down"] = (lambda event,sys=jesUncertaintyName: getattr(event,"jets_jes"+sys+"Down"),lambda event,sys=jesUncertaintyName: getattr(event,"fatjets_jes"+sys+"Down"))
    
    analyzerChain.extend(
        jetSelection(jetDict)
    )
#####

##### EVENT RECONSTRUCTION MODULE  - TO DO 
# analyzerChain.append(
#     EventReconstructionModule(
#         inputMuonCollection = OrderedDict([('tight',lambda event: event.tightRelIso_tightID_Muons), ('medium',lambda event: event.tightRelIso_mediumID_Muons), ('loose',lambda event: event.tightRelIso_looseID_Muons)]),
#         inputElectronCollection = OrderedDict([ ("MVA", OrderedDict([('tight',lambda event: event.tight_MVA_Electrons), ('medium',lambda event: event.medium_MVA_Electrons), ('loose',lambda event: event.loose_MVA_Electrons)]) ), \
#                ("cutBased", OrderedDict([('tight',lambda event: event.tight_cutBased_Electrons), ('medium',lambda event: event.medium_cutBased_Electrons), ('loose',lambda event: event.loose_cutBased_Electrons)])) ]),
#         inputJetCollection=lambda event,sys=systName: getattr(event,"selectedFatJets_"+sys),
#         inputFatJetCollection=lambda event,sys=systName: getattr(event,"selectedJets_"+sys),
#     )
# )
#####


##### GENERATION MODULES
if args.isSignal:
    analyzerChain.extend( [
        GenParticleModule_Signal(
            inputCollection=lambda event: Collection(event, "GenPart"),
            inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
            inputGenJetCollection=lambda event: Collection(event, "GenJet"),
            inputFatJetCollection= lambda event: event.selectedFatJets_nominal, #lambda event: event.diHadronic_selectedAK8s,
            inputJetCollection= lambda event: event.selectedJets_nominal, #lambda event: event.diHadronic_cleanedAK4s,
            inputMuonCollection=lambda event: event.tightRelIso_looseID_Muons,
            inputElectronCollection=lambda event: event.loose_MVA_Electrons,
            outputName="genPart",
            storeKinematics= ['pt','eta','phi','mass'],
        ),
        EventSkim(selection=lambda event: event.ngenPart == 21),
    ])

if isMC:
    analyzerChain.extend( [
        GenParticleModule(
            inputCollection=lambda event: Collection(event, "GenPart"),
            inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
            inputGenJetCollection=lambda event: Collection(event, "GenJet"),
            inputFatJetCollection= lambda event: event.selectedFatJets_nominal, #lambda event: event.diHadronic_selectedAK8s,
            inputJetCollection= lambda event: event.selectedJets_nominal, #lambda event: event.diHadronic_cleanedAK4s,
            inputMuonCollection=lambda event: event.tightRelIso_looseID_Muons,
            inputElectronCollection=lambda event: event.loose_MVA_Electrons,
            outputName="genPart",
            storeKinematics= ['pt','eta','phi','mass'],
        ),
        #EventSkim(selection=lambda event: event.ngenPart == 21),
    ])
#####

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
    maxEntries = args.maxEvents
)

p.run()

