#!/bin/bash
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
source $VO_CMS_SW_DIR/cmsset_default.sh

# cd /afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/
# eval `scramv1 runtime -sh`

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
eval `scramv1 runtime -sh`
cd "${SCRIPT_DIR}/.." || exit 1

echo "${SCRIPT_DIR}/.."
echo "ARGS:" "$@"

hostname
date
pwd
python3 trainBDT.py "$@"
date
