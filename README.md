# ttX^3 specific documentation

[![ttX^3 workflow](https://github.com/ttXcubed/nanoAOD-tools/actions/workflows/main.yml/badge.svg)](https://github.com/ttXcubed/nanoAOD-tools/actions/workflows/main.yml)

Tests are defined in [test/runCMSSWTest.sh](test/runCMSSWTest.sh):


## First steps

### Setup a fork

Go to [https://github.com/ttXcubed/nanoAOD-tools](https://github.com/ttXcubed/nanoAOD-tools) and fork the repo under your name.
You can now work and develop new things in your fork w/o interfering with others. Once finished a pull request (PR) can be 
created to merge your developments with the main repo again making them available to others. 
Please note to select the right base repository (`ttXcubed/nanoAOD-tools`) and branch (currently `ttx3`).

If you have already a fork of `nanoAOD-tools` just create a new branch for developments and merge the updates from the main repo if any.
The following adds a short cut (`ttx3-main`) to the main repo, creates a branch `ttx3` in your fork, and merges/pulls any updates from the main repo in branch `ttx3`.

```
git remote add ttx3-main https://github.com/ttXcubed/nanoAOD-tools.git
git fetch ttx3-main
git checkout -b ttx3
git pull ttx3-main ttx3
```


### Checkout instructions: CMSSW_11_1_7

```
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_11_1_7
cd CMSSW_11_1_7/src
git clone -b ttx3 git@github.com:<yourname>/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools
cmsenv
scram b
```

## Run an example

Example files:
* TTbar dilepton (ca 60k): `root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/TTTo2L2Nu_HT500Njet7_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/230000/2F43988F-FE4C-934B-B3FC-8343D149B68F.root`
* BSM signal (1.6k): `root://cms-xrd-global.cern.ch//store/user/mullerd/tt_heavyHiggs_M-500_TuneCP5_PSweights_13TeV_pythia8/RunIISummer20UL17_NanoAODv9/220206_135734/0000/NanoAODv9_1.root`
* SM 4-tops (38k): `root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v2/260000/92382329-7C4F-D241-808B-A6C678805119.root`

Full sample files can be found under [doc/SignalSamples](doc/SignalSamples.md).

The main post-processing script for producing final ntuples for plotting and statistical interpretation is:
```
python PhysicsTools/NanoAODTools/processors/ttX3.py \
    -i <root input file> \
    --year <year>  \
    [--isSignal] [--isData] [--nosys] [--invid] [--maxEvents=<N>] \
    <output directory>
```

The script accepts the following arguments:
* `-i` the input root file (nanoAODv9 or higher; needs to created from miniAODv2)
* `--year` needs to be one of the following: '2016','2016preVFP','2017','2018' (default: '2017')
* `--isSignal` optional flag to store additional information for the signal (e.g. parton/particle level observables, LHE weights) (default: false)
* `--isData` optional flag to remove gen-level information when running on data (default: false)
* `--nosys` optional flag to remove all systematics (default: false)
* `--invid` optional flag to invert the lepton ID/isolation which can be used for QCD estimation (default: false)
* `--maxEvents=<N>` optional flag to reduce the number of events for testing

Feel free to create other post-processing scripts for studies on the side!

## Analysis modules

The analysis-specific modules can be found under [python/modules](https://github.com/WbWbX/nanoAOD-tools/tree/wbwbxUL/python/modules).
In general, the input and output collections of modules should be configurable so that they can be easily reused, i.e. for evaluating systematics. A few important ones are
* [MuonSelection](python/modules/MuonSelection.py)/[MuonVeto](python/modules/MuonVeto.py): modules to select tight isolated muon candidates and veto additional loose muons (similar modules exists for electrons). For MC, the module also produces weights to evaluate the uncertainties on the lepton selection and reconstruction efficiencies.
* [JetMetUncertainties](python/modules/JetMetUncertainties.py): module to evalutate JEC/JER/MET uncertainties. Separate collections of objects are created for each variation.
* [JetSelection](python/modules/JetSelection.py): module to select jets within acceptance
* [MetFilter](python/modules/MetFilter.py): rejects events with spurious MET



# Official nanoAOD-tools documentation
Tools for working with NanoAOD (requiring only python + root, not CMSSW)

## Checkout instructions: standalone

You need to setup python 2.7 and a recent ROOT version first.

    git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git NanoAODTools
    cd NanoAODTools
    bash standalone/env_standalone.sh build
    source standalone/env_standalone.sh

Repeat only the last command at the beginning of every session.

Please never commit neither the build directory, nor the empty init.py files created by the script.

## Checkout instructions: CMSSW

    cd $CMSSW_BASE/src
    git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
    cd PhysicsTools/NanoAODTools
    cmsenv
    scram b

## General instructions to run the post-processing step

The script to run the post-processing step is `scripts/nano_postproc.py`.

The basic syntax of the command is the following:

    python scripts/nano_postproc.py /path/to/output_directory /path/to/input_tree.root

Here is a summary of its features:
* the `-s`,`--postfix` option is used to specify the suffix that will be appended to the input file name to obtain the output file name. It defaults to *_Friend* in friend mode, *_Skim* in full mode.
* the `-c`,`--cut` option is used to pass a string expression (using the same syntax as in TTree::Draw) that will be used to select events. It cannot be used in friend mode.
* the `-J`,`--json` option is used to pass the name of a JSON file that will be used to select events. It cannot be used in friend mode.
* if run with the `--full` option (default), the output will be a full nanoAOD file. If run with the `--friend` option, instead, the output will be a friend tree that can be attached to the input tree. In the latter case, it is not possible to apply any kind of event selection, as the number of entries in the parent and friend tree must be the same.
* the `-b`,`--branch-selection` option is used to pass the name of a file containing directives to keep or drop branches from the output tree. The file should contain one directive among `keep`/`drop` (wildcards allowed as in TTree::SetBranchStatus) or `keepmatch`/`dropmatch` (python regexp matching the branch name) per line, as shown in the [this](python/postprocessing/examples/keep_and_drop.txt) example file.
  * `--bi` and `--bo` allows to specify the keep/drop file separately for input and output trees.  
* the `--justcount` option will cause the script to printout the number of selected events, without actually writing the output file.

Please run with `--help` for a complete list of options.

## How to write and run modules

It is possible to import modules that will be run on each entry passing the event selection, and can be used to calculate new variables that will be included in the output tree (both in friend and full mode) or to apply event filter decisions.

We will use `python/postprocessing/examples/exampleModule.py` as an example. The module definition [file](python/postprocessing/examples/exampleModule.py), containing a simple constructor
```
   exampleModuleConstr = lambda : exampleProducer(jetSelection= lambda j : j.pt > 30)
```
should be imported using the following syntax:

```
python scripts/nano_postproc.py outDir /eos/cms/store/user/andrey/f.root -I PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule exampleModuleConstr
```

Let us now examine the structure of the `exampleProducer` module class. All modules must inherit from `PhysicsTools.NanoAODTools.postprocessing.framework.eventloop.Module`.
* the `__init__` constructor function should be used to set the module options.
* the `beginFile` function should create the branches that you want to add to the output file, calling the `branch(branchname, typecode, lenVar)` method of `wrappedOutputTree`. `typecode` should be the ROOT TBranch type ("F" for float, "I" for int etc.). `lenVar` should be the name of the variable holding the length of array branches (for instance, `branch("Electron_myNewVar","F","nElectron")`). If the `lenVar` branch does not exist already - it can happen if you create a new collection, see an example [here](python/postprocessing/examples/collectionMerger.py)) - it will be automatically created.
* the `analyze` function is called on each event. It should return `True` if the event is to be retained, `False` if it should be dropped.

### Keep/drop branches
See the effect of keep/drop instructions by running:
```
python scripts/nano_postproc.py outDir /eos/cms/store/user/andrey/f.root -I PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule exampleModuleConstr -s _exaModu_keepdrop --bi scripts/keep_and_drop_input.txt --bo scripts/keep_and_drop_output.txt
```
comparing to the previous command (without `--bi` and `--bo`).
The output branch created by _exampleModuleConstr_ produces the same result in both cases. But this one drops all other branches when creating output tree. It also runs faster.

The event interface, defined in `PhysicsTools.NanoAODTools.postprocessing.framework.datamodule`, allows to dynamically construct views of objects organized in collections, based on the branch names, for instance:

    electrons = Collection(event, "Electron")
    if len(electrons)>1: print electrons[0].someVar+electrons[1].someVar
    electrons_highpt = filter(lambda x: x.pt>50, electrons)

and this will access the elements of the `Electron_someVar`, `Electron_pt` branch arrays. Event variables can be accessed simply by `event.someVar`, for instance `event.rho`.

The output branches should be filled calling the `fillBranch(branchname, value)` method of `wrappedOutputTree`. `value` should be the desired value for single-value branches, an iterable with the correct length for array branches. It is not necessary to fill the `lenVar` branch explicitly, as this is done automatically using the length of the passed iterable.


### mht producer
Now, let's have a look at another example, `python/postprocessing/examples/mhtjuProducerCpp.py`, [file](python/postprocessing/examples/mhtjuProducerCpp.py). Similarly, it should be imported using the following syntax:

```
python scripts/nano_postproc.py outDir /eos/cms/store/user/andrey/f.root -I PhysicsTools.NanoAODTools.postprocessing.examples.mhtjuProducerCpp mhtju
```
This module has the same structure of its producer as `exampleProducer`, but in addition it utilizes a C++ code to calculate the mht variable, `src/mhtjuProducerCppWorker.cc`. This code is loaded in the `__init__` method of the producer.


