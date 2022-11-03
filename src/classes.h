#include "PhysicsTools/NanoAODTools/interface/PyJetResolutionWrapper.h"
#include "PhysicsTools/NanoAODTools/interface/PyJetResolutionScaleFactorWrapper.h"
#include "PhysicsTools/NanoAODTools/interface/PyJetParametersWrapper.h"
#include "PhysicsTools/NanoAODTools/interface/WeightCalculatorFromHistogram.h"
#include "PhysicsTools/NanoAODTools/interface/ReduceMantissa.h"

#include "PhysicsTools/NanoAODTools/src/BTagCorrlibReader.h"
#include "PhysicsTools/NanoAODTools/src/MuonCorrlibReader.h"
#include "PhysicsTools/NanoAODTools/src/ElectronCorrlibReader.h"

PyJetResolutionWrapper jetRes;
PyJetResolutionScaleFactorWrapper jetResScaleFactor;
PyJetParametersWrapper jetParams;
WeightCalculatorFromHistogram wcalc;
ReduceMantissaToNbitsRounding red(12);

BTagCorrlibReader btagSFReader;
MuonCorrlibReader muonSFReader;
ElectronCorrlibReader electronSFReader;
