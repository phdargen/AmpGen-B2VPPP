#include "AmpGen/HyperPlot/HyperBinningMakerMintRandomise.h"



HyperBinningMakerMintRandomise::HyperBinningMakerMintRandomise(const HyperCuboid& binningRange,const HyperPointSet& data) :
  HyperBinningMaker(binningRange, data)
{
  WELCOME_LOG << "Good day from the HyperBinningMakerRandomise() Constructor"<<std::endl;  
}

void HyperBinningMakerMintRandomise::makeBinning(){
  
  if (s_printBinning == true) INFO_LOG << "Splitting all bins in random dimensions"<<std::endl;
  int dimension = _binningDimensions.size();  

  int nBins = 0;
  int unchanged = 0;

  while (splitAllRandomise() != 0){
    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= 2.0*dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins"<<std::endl;
  }

  if (s_printBinning == true) INFO_LOG << "Random binning algorithm complete "<<std::endl;


}

HyperBinningMakerMintRandomise::~HyperBinningMakerMintRandomise(){
  GOODBYE_LOG << "Goodbye from the HyperBinningMakerRandomise() Constructor"<<std::endl;  
}
