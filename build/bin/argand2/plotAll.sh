 #!/bin/bash

./LineshapePlotter.exe argand2/all/XS_spline.txt --plotSpline=2 --outDir=argand2/all/XS  --Spline::NormBinBW=3 --Lineshape="GSpline" --ArgandMinY=-0.5 --ArgandMaxY=1.2 --ArgandMinX=-0.5 --ArgandMaxX=1.2
#./LineshapePlotter.exe argand2/all/XS2_spline.txt --plotSpline=2 --outDir=argand2/all/XS2  --Spline::NormBinBW=3 --Lineshape="GSpline"
./LineshapePlotter.exe argand2/all/X2S_spline.txt --plotSpline=2 --outDir=argand2/all/X2S  --Spline::NormBinBW=3 --Lineshape="GSpline" --ArgandMinY=-0.7 --ArgandMaxY=1.7 --ArgandMinX=-1.2 --ArgandMaxX=1.2 --AmpMax=1.5
./LineshapePlotter.exe argand2/all/XA_spline.txt --plotSpline=2 --outDir=argand2/all/XA  --Spline::NormBinBW=5 --Lineshape="GSpline"  --ArgandMinY=-0.7 --ArgandMaxY=1.7 --ArgandMinX=-1.2 --ArgandMaxX=1.2 --AmpMax=1.5
./LineshapePlotter.exe argand2/all/XV_spline.txt --plotSpline=2 --outDir=argand2/all/XV  --Spline::NormBinBW=5 --Lineshape="GSpline" --ArgandMinY=-2 --ArgandMaxY=2 --ArgandMinX=-2 --ArgandMaxX=2 --AmpMax=2

./LineshapePlotter.exe argand2/all/XsS_spline.txt --plotSpline=2 --outDir=argand2/all/XsS  --Spline::NormBinBW=1 --ArgandLine=L --Lineshape="GSpline" --fitBW=0 --ArgandMinY=-0.6 --ArgandMaxY=1.3 --ArgandMinX=-0.6 --ArgandMaxX=1.3
#./LineshapePlotter.exe argand2/all/XsA2_spline.txt --plotSpline=2 --outDir=argand2/all/XsA2  --Spline::NormBinBW=4
./LineshapePlotter.exe argand2/all/Xs2A_spline.txt --plotSpline=2 --outDir=argand2/all/Xs2A  --Spline::NormBinBW=4 --ArgandLine=L --Lineshape="GSpline" --ArgandMinY=-0.5 --ArgandMaxY=1.2 --ArgandMinX=-0.5 --ArgandMaxX=1.2 
./LineshapePlotter.exe argand2/all/XsA_spline.txt --plotSpline=2 --outDir=argand2/all/XsA  --Spline::NormBinBW=3 --AmpMax=2 --Lineshape="GSpline" --ArgandMinY=-1. --ArgandMaxY=1.6 --ArgandMinX=-1.3 --ArgandMaxX=1.3 --ArgandLine=L --fitBW=0
./LineshapePlotter.exe argand2/all/XsV_spline.txt --plotSpline=2 --outDir=argand2/all/XsV  --Spline::NormBinBW=4 --Lineshape="GSpline" --AmpMax=2 --ArgandMinY=-0.3 --ArgandMaxY=1.5 --ArgandMinX=-0.5 --ArgandMaxX=1.3  --ArgandLine=L

./LineshapePlotter.exe argand2/all/Z_spline.txt --plotSpline=2 --outDir=argand2/all/Z  --Spline::NormBinBW=3 --ArgandMinY=-0.4 --ArgandMaxY=1.3 --ArgandMinX=-0.7 --ArgandMaxX=1.2
./LineshapePlotter.exe argand2/all/Z2_spline.txt --plotSpline=2 --outDir=argand2/all/Z2  --Spline::NormBinBW=3 --fitBW=0  --AmpMax=2 --ArgandMinY=-0.3 --ArgandMaxY=1.5 --ArgandMinX=-0.6 --ArgandMaxX=1.2
./LineshapePlotter.exe argand2/all/Zs_spline.txt --plotSpline=2 --outDir=argand2/all/Zs  --Spline::NormBinBW=3 --AmpMax=1 --ArgandMinY=-0.2 --ArgandMaxY=.8 --ArgandMinX=-0.5 --ArgandMaxX=0.5
