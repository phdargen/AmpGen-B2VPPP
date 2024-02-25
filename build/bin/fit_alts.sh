#!/bin/bash                                                                                                                                                                                       

#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_addK.txt myTest3/modelBest5_fixed.txt freeResoParams.txt --outDir=altModels/1  --addAmpList=""
#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_addExotic.txt myTest3/modelBest5_fixed.txt freeResoParams.txt --outDir=altModels/2  --addAmpList=""

#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extK.txt myTest3/modelBest5_fixed.txt freeResoParams.txt --outDir=altModels/3  --addAmpList=""
./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extExotic.txt myTest3/modelBest5_fixed.txt freeResoParams.txt --outDir=altModels/4  --addAmpList=""
./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extNR.txt myTest3/modelBest5_fixed.txt freeResoParams.txt --outDir=altModels/5  --addAmpList=""

#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extX.txt --outDir=altModels/extX  --addAmpList=""
#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extXs.txt --outDir=altModels/extXs  --addAmpList=""
#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extZ.txt --outDir=altModels/extZ  --addAmpList=""
#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extZ2.txt --outDir=altModels/extZ2  --addAmpList=""
#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extZs.txt --outDir=altModels/extZs  --addAmpList="" 

#./myFitter.exe main_psi.txt bestModel5.txt altModels/altModel_extExotic.txt --outDir=altModels/extExotic  --addAmpList=""

#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_addK.txt --outDir=altModels/addK  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_addExotic.txt --outDir=altModels/addExotic  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1           
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extK.txt --outDir=altModels/extK  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1      
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extK2.txt --outDir=altModels/extK2  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1                       
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extNR.txt --outDir=altModels/extNR  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1                         
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extX.txt --outDir=altModels/extX  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1                         
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extXs.txt --outDir=altModels/extXs  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1                         
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extZ.txt --outDir=altModels/extZ  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extZ2.txt --outDir=altModels/extZ2  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1                  
#./myPlotter.exe main_psi.txt bestModel5.txt altModels/altModel_extZs.txt --outDir=altModels/extZs  --addAmpList="" --nBins=70 --markerSize=0.75 --lineWidthFit=4 --lineWidthAmps=2 --smoothFitProj=2 --addChi2PerPlot=1
