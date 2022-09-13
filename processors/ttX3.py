import os
import sys
import math
import argparse
import random
import ROOT
import numpy as np

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
    "year": args.year
}

Module.globalOptions = globalOptions

isMC = not args.isData
isPowheg = 'powheg' in args.inputFiles[0].lower()
isPowhegTTbar = 'TTTo' in args.inputFiles[0] and isPowheg

minMuonPt =     {'2016': 15., '2016preVFP': 15., '2017': 15., '2018': 15.}
minElectronPt = {'2016': 15., '2016preVFP': 15., '2017': 15., '2018': 15.}


jesUncertaintyFilesRegrouped = {
    '2016':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt",
    '2016preVFP': "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt",
    '2017':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt",
    '2018':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt"
}
jerResolutionFiles = {
    '2016':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2016preVFP': "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2017':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK4PFchs.txt",
    '2018':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt"
}
jerSFUncertaintyFiles = {
    '2016':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt",
    '2016preVFP': "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt",
    '2017':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK4PFchs.txt",
    '2018':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt"
}



jesAK8UncertaintyFilesRegrouped = {
    '2016':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16_V7_MC_UncertaintySources_AK8PFPuppi.txt",
    '2016preVFP': "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL16APV_V7_MC_UncertaintySources_AK8PFPuppi.txt",
    '2017':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_V5_MC_UncertaintySources_AK8PFPuppi.txt",
    '2018':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.txt"
}
jerAK8ResolutionFiles = {
    '2016':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2016preVFP': "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2017':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.txt",
    '2018':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt"
}
jerAK8SFUncertaintyFiles = {
    '2016':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16_JRV3_MC_SF_AK8PFPuppi.txt",
    '2016preVFP': "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.txt",
    '2017':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL17_JRV3_MC_SF_AK8PFPuppi.txt",
    '2018':       "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt"
}

def leptonSequence():
    seq = [
        MuonSelection(
            inputCollection=lambda event: Collection(event, "Muon"),
            outputName_list=["tightMuons","mediumMuons","looseMuons"],
            storeKinematics=['pt','eta','charge','phi','mass'],
            storeWeights=True,
            muonMinPt=minMuonPt[args.year],
            muonMaxEta=2.4,
            triggerMatch=True,
            #muonID= MuonSelection.TIGHT,
            #muonIso= MuonSelection.INV if args.invid else MuonSelection.TIGHT,
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
            outputName_list = ["tightElectrons","mediumElectrons","looseElectrons"],
            #electronID = ElectronSelection.INV if args.invid else ElectronSelection.WP90,
            electronMinPt = minElectronPt[args.year],
            electronMaxEta = 2.4,
            storeKinematics=['pt','eta','charge','phi','mass'],
            storeWeights=True,
        ),

        #ElectronVeto(
        #    inputCollection=lambda event: event.unselectedElectrons,
        #    outputName = "vetoElectrons",
        #    electronMinPt = 10.,
        #    electronMaxEta = 2.4,
        #    storeKinematics=['pt','eta','charge','phi','mass'],
        #),
	
	#EventSkim(selection=lambda event: (len(event.tightMuons) > 0 or len(event.tightElectrons)) > 0 ),
	#EventSkim(selection=lambda event: (len(event.looseMuons) > 0 or len(event.looseElectrons)) > 0 ),
	EventSkim(selection=lambda event: (len(event.looseMuons) + len(event.looseElectrons)) == 2 ),

        DoubleLeptonTriggerSelection(
            outputName="Trigger",
            storeWeights=False,
            thresholdPt=15. #minDiMuonPt[args.year], PRELIMINARY VALUE
         ),
 
 	EventSkim(selection=lambda event: (event.Trigger_general_flag)),
 
 	#TriggerMatching(
 	#	inputCollectionMuon = lambda event: event.tightMuons,
        #	inputCollectionElectron = lambda event: event.tightElectrons,
        #	storeWeights=False,
       # 	outputName = "TriggerObjectMatching", 
        #	triggerMatch=True,
        #	thresholdPt=15.
 	#),       
        
        #EventSkim(selection=lambda event: (event.Trigger_general_flag and event.TriggerObjectMatching_flag )),
        #EventSkim(selection=lambda event: (len(event.looseMuons) + len(event.looseElectrons)) == 0),
        
    ]
    return seq


    
def jetSelection(jetDict):
    seq = []
    
    for systName,(jetCollection,fatjetCollection) in jetDict.items():
        seq.extend([
            JetSelection(
                inputCollection= jetCollection, #lambda event: Collection(event, "Jet"),
                leptonCollectionDRCleaning=lambda event: event.looseMuons+event.looseElectrons,#event.tightMuons+event.tightElectrons,
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
                inputCollection= fatjetCollection, #lambda event: Collection(event, "FatJet"),
                leptonCollectionDRCleaning=lambda event,sys=systName: event.looseMuons+event.looseElectrons,
                jetMinPt=400., #ak8 only stored > 175 GeV
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
                inputCollection=lambda event,sys=systName: getattr(event,"selectedJets_"+sys),
                flagName="isBTagged",
                outputName_list=["selectedBJets_"+systName+"_tight", "selectedBJets_"+systName+"_medium", "selectedBJets_"+systName+"_loose"],
                jetMinPt=30.,
                jetMaxEta=2.4,
                workingpoint = [],#BTagSelection.TIGHT,
                storeKinematics=['pt', 'eta','phi','mass'],
                storeTruthKeys = ['hadronFlavour','partonFlavour'],
            )
        )
        
    systNames = jetDict.keys()
   
    #at least 4 AK4 jets and 2 AK8 jets
    seq.append(
        EventSkim(selection=lambda event, systNames=systNames: 
            any([getattr(event, "nselectedJets_"+systName) >= 2 for systName in systNames])
        ),
    )
   
    
    #seq.append(
    #    EventSkim(selection=lambda event, systNames=systNames: 
    #        any([len(filter(lambda jet: jet.isBTagged,getattr(event,"selectedJets_"+systName))) >= 2 for systName in systNames])
    #   )
    #)

    seq.append(
        EventSkim(selection=lambda event, systNames=systNames: 
            any([getattr(event, "nselectedBJets_"+systName+"_loose") >= 2 for systName in systNames])
        ),
    )
	    
    #at least 2 AK8 jets
    seq.append(
        EventSkim(selection=lambda event, systNames=systNames: 
            any([getattr(event, "nselectedFatJets_"+systName) >= 2 for systName in systNames])
        )
    )
    
    #TODO: btagging SF producer might have a bug
    '''
    if isMC:
        jesUncertForBtag = ['jes'+syst.replace('Total','') for syst in jesUncertaintyNames]
        # to remove once breakdown available
        if args.year != '2016preVFP':
            jesUncertForBtag = ['jes']
        seq.append(
            btagSFProducer(
                era=args.year,
                jesSystsForShape = jesUncertForBtag,
                nosyst = args.nosys
            )
        )
    '''

    
            
    return seq

storeVariables = [[lambda tree: tree.branch("genweight", "F"),
                           lambda tree,
                           event: tree.fillBranch("genweight",
                           event.Generator_weight)],
    #[lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
     #event: tree.fillBranch("PV_npvs", event.PV_npvs)],
    #[lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
     #event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
    #[lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,
     #event: tree.fillBranch("fixedGridRhoFastjetAll",
                            #event.fixedGridRhoFastjetAll)],
]
		                   
analyzerChain = [EventInfo(storeVariables=storeVariables),
EventSkim(selection=lambda event: event.nTrigObj > 0),
    MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
    ),
]

analyzerChain.extend(leptonSequence())


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
            jetKeys=['jetId', 'btagDeepB','deepTag_TvsQCD','deepTag_WvsQCD','particleNet_TvsQCD','particleNet_WvsQCD','particleNet_mass', 'particleNet_QCD','hadronFlavour','genJetAK8Idx', 'nBHadrons'],  #'nConstituents'
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


if not args.isData:
    #analyzerChain.append(GenWeightProducer())
    if isPowhegTTbar:
        analyzerChain.append(
            TopPtWeightProducer(
                mode=TopPtWeightProducer.DATA_NLO
            )
        )


storeVariables = [
    [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
     event: tree.fillBranch("PV_npvs", event.PV_npvs)],
    [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
     event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
    [lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,
     event: tree.fillBranch("fixedGridRhoFastjetAll",
                            event.fixedGridRhoFastjetAll)],
]


if not globalOptions["isData"]:
    storeVariables.append([lambda tree: tree.branch("genweight", "F"),
                           lambda tree,
                           event: tree.fillBranch("genweight",
                           event.Generator_weight)])

    '''
    L1prefirWeights =  ['Dn', 'Nom', 'Up', 'ECAL_Dn', 'ECAL_Nom', 'ECAL_Up',
                        'Muon_Nom', 'Muon_StatDn', 'Muon_StatUp', 'Muon_SystDn', 'Muon_SystUp']

    for L1prefirWeight in L1prefirWeights:
        storeVariables.append([
            lambda tree, L1prefirWeight=L1prefirWeight: tree.branch('L1PreFiringWeight_{}'.format(L1prefirWeight.replace('Dn','Down').replace('Nom','Nominal')), "F"),
            lambda tree, event, L1prefirWeight=L1prefirWeight: tree.fillBranch('L1PreFiringWeight_{}'.format(L1prefirWeight.replace('Dn','Down').replace('Nom','Nominal')),
                                                                               getattr(event,'L1PreFiringWeight_{}'.format(L1prefirWeight)))
        ])
    '''
    #analyzerChain.append(EventInfo(storeVariables=storeVariables))

##############################SISTEMA PER LE SYSTEMATICHE

#analyzerChain.append( 
#	EventReconstruction(
#		inputMuonCollection=lambda event: [event.tightMuons, event.mediumMuons, event.looseMuons],
#        	inputElectronCollection=lambda event: [event.tightElectrons, event.mediumElectrons, event.looseElectrons],
#        	inputJetCollection=lambda event: event.selectedJets_nominal,
#        	inputFatJetCollection=lambda event: event.selectedFatJets_nominal,
#        	inputBJetCollection=lambda event: [event.selectedBJets_nominal_tight, event.selectedBJets_nominal_medium, event.selectedBJets_nominal_loose],  #ALL THE BJETS ARE AK4 JETS
#        	inputMETCollection=lambda event: event.met_nominal,
#        	outputName="diHadronic",
#        	storeKinematics_jets= ['pt','eta','phi','mass'],
#        	storeKinematics_leptons= ['pt','eta','phi','mass','charge','leptonFlavour'],
#	)		
#)


if args.isSignal:
	analyzerChain.extend( [
		GenParticleModule_Signal(
			inputCollection=lambda event: Collection(event, "GenPart"),
			inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
			inputGenJetCollection=lambda event: Collection(event, "GenJet"),
			inputFatJetCollection= lambda event: event.selectedFatJets_nominal, #lambda event: event.diHadronic_selectedAK8s,
        		inputJetCollection= lambda event: event.selectedJets_nominal, #lambda event: event.diHadronic_cleanedAK4s,
			inputMuonCollection=lambda event: event.looseMuons,
			inputElectronCollection=lambda event: event.looseElectrons,
			#eventReco_flags = lambda event: event.diHadronic_event_selection_flags,
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
			inputMuonCollection=lambda event: event.looseMuons,
			inputElectronCollection=lambda event: event.looseElectrons,
			#inputMuonCollection
			#eventReco_flags = lambda event: event.diHadronic_event_selection_flags,
			outputName="genPart",
			storeKinematics= ['pt','eta','phi','mass'],
		),
		#EventSkim(selection=lambda event: event.ngenPart == 21),
	])

'''
analyzerChain.extend([
	EventSkim(selection=lambda event: event.MET_energy > 200),
	EventSkim(selection=lambda event: event.selectedJets_nominal_HT > 1000)
])
'''

p = PostProcessor(
    args.output[0],
    args.inputFiles,
    cut="",#"(nJet>1)&&((nElectron+nMuon)>1)", #at least 2 jets + 2 leptons
    modules=analyzerChain,
    friend=True,
    maxEntries = args.maxEvents
)

p.run()

