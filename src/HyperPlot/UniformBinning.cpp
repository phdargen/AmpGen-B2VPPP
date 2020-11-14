#include "AmpGen/HyperPlot/UniformBinning.h"


///The only constructor
UniformBinning::UniformBinning(HyperCuboid limits, int nLocalBins) :
  _limits    (limits),
  _nLocalBins(_limits.getDimension(), nLocalBins)
{
  setBinningType("UniformBinning");
  setDimension  (_limits.getDimension());
  
  WELCOME_LOG << "Hello from the UniformBinning() Constructor";
}


UniformBinning::UniformBinning(HyperCuboid limits, std::vector<int> nLocalBins) :
  _limits    (limits),
  _nLocalBins(nLocalBins)
{
  setBinningType("UniformBinning");
  setDimension  (_limits.getDimension());

  WELCOME_LOG << "Hello from the UniformBinning() Constructor";
}



int UniformBinning::getNumLocalBins(int dimension) const{
  return _nLocalBins.at(dimension); 
}

int UniformBinning::getNumBins() const{
  int nBins = 1;
  for (int i = 0; i < getDimension(); i++){
    nBins *= getNumLocalBins(i);
  }
  return nBins;
}


int UniformBinning::getGlobalBinNumber( std::vector<int> localBinNumbers ) const{

  int dimension = getDimension();

  int multiplier = 1;
  int binNumber  = 0;

  for (int i = 0; i < dimension; i++){
    binNumber += localBinNumbers.at(i)*multiplier;
    multiplier *= getNumLocalBins(i);
  }
  
  return binNumber;

}

std::vector<int> UniformBinning::getLocalBinNumbers( int globalBinNumber ) const{

  int nBins = getNumBins();

  if ( globalBinNumber >= nBins || globalBinNumber < 0 ){
    ERROR_LOG << "UniformBinning::getLocalBinNumbers - The global bin number you have given is out of range = " << globalBinNumber << std::endl;
  }
 
  int dimension = getDimension();
  
  std::vector<int> localBinNums(dimension, -1.0);

  int multiplier = nBins;

  for (int i = dimension - 1; i >= 0; i--){

    multiplier = multiplier/getNumLocalBins(i); 

    int localBinNum = floor(double(globalBinNumber)/double(multiplier));
    
    globalBinNumber -= multiplier*localBinNum;

    localBinNums.at(i) = localBinNum;
  }

  return localBinNums;

}

double UniformBinning::getLowBinEdgeLocal(int dim, int localBinNum) const{
  
  double low    = _limits.getLowCorner ().at(dim);
  double high   = _limits.getHighCorner().at(dim);
  double nbins  = _nLocalBins            .at(dim);
  double width = (high - low)/nbins;

  return low + width*localBinNum;

}

double UniformBinning::getHighBinEdgeLocal(int dim, int localBinNum) const{
  
  double low    = _limits.getLowCorner ().at(dim);
  double high   = _limits.getHighCorner().at(dim);
  double nbins  = _nLocalBins            .at(dim);
  double width = (high - low)/nbins;

  return low + width*(localBinNum + 1.0);

}

HyperPoint UniformBinning::getLowCorner(int globalBinNum) const{
  HyperPoint corner(getDimension());
  std::vector<int> localBinNum = getLocalBinNumbers(globalBinNum);
  for (int i = 0; i < getDimension(); i++){
    corner.at( i ) = getLowBinEdgeLocal( i, localBinNum.at(i) );
  }
  return corner;
}

HyperPoint UniformBinning::getHighCorner(int globalBinNum) const{
  HyperPoint corner(getDimension());
  std::vector<int> localBinNum = getLocalBinNumbers(globalBinNum);
  for (int i = 0; i < getDimension(); i++){
    corner.at( i ) = getHighBinEdgeLocal( i, localBinNum.at(i) );
  }
  return corner;
}

void UniformBinning::load(TString filename, TString option){
  //not implemented this yet so added these lines to stop
  //compiler warnings
  filename = filename;
  option = option;
}

BinningBase* UniformBinning::clone() const{
  return new UniformBinning(*this);
}

void UniformBinning::save(TString filename) const{
  //not implemented this yet so added these lines to stop
  //compiler warnings
  filename = filename;
}

void UniformBinning::save() const{

}

void UniformBinning::mergeBinnings( const BinningBase& other ){
  //not implemented this yet so added these lines to stop
  //compiler warnings
  other.getDimension();
}

int UniformBinning::getLocalBinNumber(int dim, double val) const{
  
  double low    = _limits.getLowCorner ().at(dim);
  double high   = _limits.getHighCorner().at(dim);
  double nbins  = _nLocalBins            .at(dim);
  double width = (high - low)/nbins; 
  
  return floor( val - low )/width; 

}

std::vector<int> UniformBinning::getLocalBinNumbers(const HyperPoint& coords) const{
  
  std::vector<int> localBinNums(getDimension(), -1);

  for (int i = 0; i < getDimension(); i++){
    localBinNums.at(i) = getLocalBinNumber(i, coords.at(i));
  }
  return localBinNums;

}


int UniformBinning::getBinNum(const HyperPoint& coords) const{
  
  return getGlobalBinNumber( getLocalBinNumbers(coords) );

}

HyperVolume UniformBinning::getBinHyperVolume(int binNumber) const{
  HyperCuboid cube( getLowCorner(binNumber), getHighCorner(binNumber) );
  return HyperVolume(cube);
}



HyperPoint  UniformBinning::getAverageBinWidth() const{

  HyperPoint binwidth(_limits.getHighCorner() - _limits.getLowCorner());
  for (int i = 0; i < getDimension(); i++){
    binwidth.at(i)/= double(getNumLocalBins(i));
  }
  return binwidth;

}

HyperCuboid UniformBinning::getLimits()          const{
  return _limits;
}

UniformBinning::~UniformBinning(){


}


