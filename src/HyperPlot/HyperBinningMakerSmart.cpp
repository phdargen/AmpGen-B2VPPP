#include "AmpGen/HyperPlot/HyperBinningMakerSmart.h"


HyperBinningMakerSmart::HyperBinningMakerSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerSmart() Constructor";  
}

void HyperBinningMakerSmart::makeBinning(){
   
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim; 
  if (splitDim >= dimension) splitDim = 0; 
  
  int nBins = 0;
  int unchanged = 0;  

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim) << std::endl;

  while (smartSplitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0; 
    if (unchanged > dimension) break;
    nBins = getNumBins(); 
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<< std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<< std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Smart binning algorithm complete"<< std::endl;

}

HyperBinningMakerSmart::~HyperBinningMakerSmart(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerSmart() Constructor";  
}


