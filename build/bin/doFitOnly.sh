#!/bin/bash
export AMPGENROOT="/afs/cern.ch/work/p/phdargen/lhcb-analysis-b2psikpipi/AmpGen/"
export HOME="/afs/cern.ch/user/p/phdargen/"
#source "/cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_96python3 x86_64-centos7-gcc8-opt"
source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
#$AMPGENROOT/bin/build/myFitter "$@"
#echo "Executing $AMPGENROOT/build/bin/myFitter $@"
#echo "$@"
#$AMPGENROOT/build/bin/myFitter "$@"
touch result.root
touch log.txt
touch model.txt
touch Fit_weights.root
./myFitter.exe "$@"
#./myPlotter.exe "$@"
