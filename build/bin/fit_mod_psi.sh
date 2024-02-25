#!/bin/bash                                                                                                                                                                                                                    
./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_omnes.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0 
./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_omnes_1.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_1.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_omnes_2.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_2.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_omnes_3.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
./myFitter.exe main_psi.txt bestNoExotic.txt mod_psi/mod_psi_3.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

##
#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_omnes.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_omnes_1.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_1.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_omnes_2.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_2.txt --outDir=mod_psi --addAmpList="" --saveWeights=0

#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_omnes_3.txt omnes/omnes.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
#./myFitter.exe main_psi.txt bestModel5.txt mod_psi/mod_psi_3.txt --outDir=mod_psi --addAmpList="" --saveWeights=0
