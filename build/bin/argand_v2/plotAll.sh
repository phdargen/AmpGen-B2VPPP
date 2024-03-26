 #!/bin/bash

#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/XS_spline_2.txt --plotSpline=2 --outDir=argand_v2/all/XS  --Spline::NormBinBW=3 --Lineshape="GSpline" --ArgandMinY=-0.3 --ArgandMaxY=1.4 --ArgandMinX=-0.5 --ArgandMaxX=1.2 --ArgandLine=L --fig_name="Fig11"
#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/X2S_spline.txt --plotSpline=2 --outDir=argand_v2/all/X2S  --Spline::NormBinBW=3 --Lineshape="GSpline" --ArgandMinY=-0.5 --ArgandMaxY=1.5 --ArgandMinX=-1 --ArgandMaxX=1. --AmpMax=1.5 --ArgandLine=L  --fig_name="Fig13"
#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/XA_spline.txt --plotSpline=2 --outDir=argand_v2/all/XA  --Spline::NormBinBW=5 --Lineshape="GSpline" --ArgandMinY=-0.6 --ArgandMaxY=1.7 --ArgandMinX=-1.1 --ArgandMaxX=1.2 --AmpMax=1.5 --ArgandLine=L --fig_name="Fig12"
#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/XV_spline.txt --plotSpline=2 --outDir=argand_v2/all/XV  --Spline::NormBinBW=5 --Lineshape="GSpline" --ArgandMinY=-0.5 --ArgandMaxY=1.5 --ArgandMinX=-0.5 --ArgandMaxX=1.5 --AmpMax=2 --ArgandLine=L --fig_name="Fig14"

#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/Xs2A_spline.txt --plotSpline=2 --outDir=argand_v2/all/Xs2A  --Spline::NormBinBW=4 --ArgandLine=L --Lineshape="GSpline" --ArgandMinY=-0.5 --ArgandMaxY=1.2 --ArgandMinX=-0.5 --ArgandMaxX=1.2 --AmpMax=1.4 --ArgandLine=L  --fig_name="Fig16"
#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/XsA_spline.txt --plotSpline=2 --outDir=argand_v2/all/XsA  --Spline::NormBinBW=3 --AmpMax=2 --Lineshape="GSpline" --ArgandMinY=-.8 --ArgandMaxY=1.6 --ArgandMinX=-1.1 --ArgandMaxX=1.3 --ArgandLine=L --fig_name="Fig15"
#./LineshapePlotter.exe rw_X_phsp.txt argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/XsV_spline.txt --plotSpline=2 --outDir=argand_v2/all/XsV  --Spline::NormBinBW=4 --Lineshape="GSpline" --ArgandMinY=-0.2 --ArgandMaxY=1. --ArgandMinX=-0.2 --ArgandMaxX=1.  --ArgandLine=L --fig_name="Fig17"

#./LineshapePlotter.exe argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/Z_spline.txt --plotSpline=2 --outDir=argand_v2/all/Z --Spline::NormBinBW=3 --ArgandMinY=-0.3 --ArgandMaxY=1.4 --ArgandMinX=-1 --ArgandMaxX=0.9 --ArgandLine=L --fig_name="Fig18"
./LineshapePlotter.exe argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/Z2_spline.txt --plotSpline=2 --outDir=argand_v2/all/Z2 --Spline::NormBinBW=3 --fitBW=0  --AmpMax=2 --ArgandMinY=-0.3 --ArgandMaxY=1.5 --ArgandMinX=-0.6 --ArgandMaxX=1.2 --fig_name="Fig19"
#./LineshapePlotter.exe argand_v2/model_baselineModel_v2_fixed.txt argand_v2/all/Zs_spline.txt --plotSpline=2 --outDir=argand_v2/all/Zs --Spline::NormBinBW=3 --AmpMax=1 --ArgandMinY=-0.1 --ArgandMaxY=.7 --ArgandMinX=-0.5 --ArgandMaxX=0.3 --ArgandLine=L --fig_name="Fig20"

