#include "AmpGen/HyperPlot/HyperBinningMakerSmartLikelihood.h"


HyperBinningMakerSmartLikelihood::HyperBinningMakerSmartLikelihood(const HyperCuboid& binningRange,const HyperPointSet& data) :
  HyperBinningMaker(binningRange, data)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerSmartLikelihood() Constructor"<<std::endl;  
}

void HyperBinningMakerSmartLikelihood::makeBinning(){
  
  int dimension = _binningDimensions.size();  
  
  int nBins = 0;
  int unchanged = 0;

  while (smartLikelihoodSplitAll() != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= 2.0*dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins" << std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "likelihood binning algorithm complete " << std::endl;

}

HyperBinningMakerSmartLikelihood::~HyperBinningMakerSmartLikelihood(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerSmartLikelihood() Constructor"<<std::endl; 
}
