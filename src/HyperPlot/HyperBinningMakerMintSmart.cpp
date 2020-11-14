#include "AmpGen/HyperPlot/HyperBinningMakerMintSmart.h"



HyperBinningMakerMintSmart::HyperBinningMakerMintSmart(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim) :
  HyperBinningMaker(binningRange, data),
  _startingDim(startingDim)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerMintSmart() Constructor"<<std::endl; 
  //makeBinning(startingDim);
}

void HyperBinningMakerMintSmart::makeBinning(){
  
  int dimension = _binningDimensions.size();
  
  int splitDim = _startingDim;
  if (splitDim >= dimension) splitDim = 0;

  int nBins = 0;
  int unchanged = 0;

  if (s_printBinning == true) INFO_LOG << "Splitting all bins in dimension " << _binningDimensions.at(splitDim)<<std::endl;

  while ( splitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<<std::endl;
  }

  nBins = 0;
  unchanged = 0;

  while (smartSplitAll(_binningDimensions.at(splitDim), 0.5) != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged > dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
    splitDim++;
    if( splitDim == dimension ) splitDim = 0;
    if (s_printBinning == true) INFO_LOG << "Trying to split all bins in dimension " << _binningDimensions.at(splitDim)<<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Mint binning algorithm complete "<<std::endl;

}

HyperBinningMakerMintSmart::~HyperBinningMakerMintSmart(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerMintSmart() Constructor"<<std::endl;  
}
