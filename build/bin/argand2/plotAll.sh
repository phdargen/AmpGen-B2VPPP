#!/bin/bash

./LineshapePlotter.exe argand2/all/XS_spline.txt --plotSpline=2 --outDir=argand2/all/XS  --Spline::NormBinBW=3
./LineshapePlotter.exe argand2/all/XS2_spline.txt --plotSpline=2 --outDir=argand2/all/XS2  --Spline::NormBinBW=3
./LineshapePlotter.exe argand2/all/X2S_spline.txt --plotSpline=2 --outDir=argand2/all/X2S  --Spline::NormBinBW=3
./LineshapePlotter.exe argand2/all/XA_spline.txt --plotSpline=2 --outDir=argand2/all/XA  --Spline::NormBinBW=5

./LineshapePlotter.exe argand2/all/XsS_spline.txt --plotSpline=2 --outDir=argand2/all/XsS  --Spline::NormBinBW=1 --ArgandLine=L
./LineshapePlotter.exe argand2/all/XsA2_spline.txt --plotSpline=2 --outDir=argand2/all/XsA2  --Spline::NormBinBW=4
./LineshapePlotter.exe argand2/all/Xs2A_spline.txt --plotSpline=2 --outDir=argand2/all/Xs2A  --Spline::NormBinBW=4 --ArgandLine=L --ArgandMinX=-1
./LineshapePlotter.exe argand2/all/XsA_spline.txt --plotSpline=2 --outDir=argand2/all/XsA  --Spline::NormBinBW=3 --AmpMax=2 --ArgandMaxY=2
./LineshapePlotter.exe argand2/all/XsV_spline.txt --plotSpline=2 --outDir=argand2/all/XsV  --Spline::NormBinBW=4

./LineshapePlotter.exe argand2/all/Z_spline.txt --plotSpline=2 --outDir=argand2/all/Z  --Spline::NormBinBW=3
