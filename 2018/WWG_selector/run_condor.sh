#!/bin/bash
ulimit -s unlimited
set -e
cd /afs/cern.ch/work/y/yian/work/GZ/CMSSW_10_6_20/src/PhysicsTools/NanoAODTools/WWG/2018/WWG_selector
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
line=`sed -n -e "$1p" file `
echo $line
python WWG_postproc.py $line
