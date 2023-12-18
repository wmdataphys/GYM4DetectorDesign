#!/usr/bin/sh

source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh -n new
cd /work/
root -q 'Fun4All_G4_EICDetector.C("PARAMS", "SETTINGS", NEVENTS, ETAMIN, ETAMAX, MOMMIN, MOMMAX, "OUTPUTFILE")'
if [ NEVENTS -gt -1 ]
then
	root -q 'analysis_resolution.C("ExecuteSingularity_g4tracking_eval.root", 0, 1.4)'
fi
exit
