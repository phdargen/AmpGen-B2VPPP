#!/bin/bash

#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_XsS.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_XsA.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_Xs2A.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_XsV.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0

./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_psi.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_XS.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_X2S.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_XA.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0

#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_ZA.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_Z2A.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_Zs.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_thresholds/mod_Zs2.txt --outDir=mod_thresholds --addAmpList="" --saveWeights=0
