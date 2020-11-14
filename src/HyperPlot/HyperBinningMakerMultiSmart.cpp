#include "AmpGen/HyperPlot/HyperBinningMakerMultiSmart.h"



HyperBinningMakerMultiSmart::HyperBinningMakerMultiSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerMultiSmart() Constructor";  
  //makeBinning(startingDim);
}

void HyperBinningMakerMultiSmart::makeBinning(){
   
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim; 
  if (splitDim >= dimension) splitDim = 0; 
  
  int nBins = 0;
  int unchanged = 0;  

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim) << std::endl;

  while (smartMultiSplitAll(_binningDimensions.at(splitDim)) != -1){
    finishedIteration();
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0; 
    if (unchanged > dimension) break;
    nBins = getNumBins(); 
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<< std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<< std::endl;
  }
  
  INFO_LOG << "Gone as far as possible with MultiSplit - now trying SmartSplit " << std::endl;
  unchanged = 0;  

  while (smartSplitAll(_binningDimensions.at(splitDim), 0.5) != -1){
    finishedIteration();
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

HyperBinningMakerMultiSmart::~HyperBinningMakerMultiSmart(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerMultiSmart() Constructor";  
}
