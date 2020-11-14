#include "AmpGen/HyperPlot/CyclicPhaseBins.h"

CyclicPhaseBins::CyclicPhaseBins(){

}

CyclicPhaseBins::CyclicPhaseBins(std::vector<double> binEdges){
  
  double min = binEdges.at(0);
  for (unsigned i = 1; i < binEdges.size(); i++){

    double val = binEdges.at(i);
    if ( val < min ){
      ERROR_LOG << "CyclicPhaseBins::CyclicPhaseBins The bin edges MUST increase" << std::endl;
    }
    if ( val > min + TMath::Pi()*2.0 ){
      ERROR_LOG << "CyclicPhaseBins::CyclicPhaseBins This bin edge is more than 2pi bigger than the smallest bin edge!" << std::endl;
    }
  
  }

  _binEdges = binEdges;  

}
  
int CyclicPhaseBins::getBinNumber(double phase) const{
  
  double min = _binEdges.at(0);  
  
  bool tooSmall   = (phase - min) < 0.0;
  bool tooBig     = (phase - min) > TMath::Pi()*2.0;

  while (tooSmall || tooBig){
    if (tooSmall) {
      phase += 2.0*TMath::Pi();
      tooSmall = (phase - min) < 0.0;
    }
    if (tooBig  ) {
      phase -= 2.0*TMath::Pi();
      tooBig = (phase - min) > TMath::Pi()*2.0;
    }
  }
  
  for (unsigned i = 1; i < _binEdges.size(); i++){
    if (phase < _binEdges.at(i)) return i - 1; 
  }

  return _binEdges.size() - 1;
 
}

void CyclicPhaseBins::setUniformCiSiBins(int nBinPairs               ){
  
  std::vector<double> binEdges;

  for (int i = 0; i < 2*nBinPairs; i++){
    binEdges.push_back( - TMath::Pi() + (TMath::Pi()/double(nBinPairs))*i );
  }

  _binEdges = binEdges;
}

void CyclicPhaseBins::setBinEdges       (std::vector<double> binEdges){
  _binEdges = binEdges; 
}

double CyclicPhaseBins::getLowBinBoundary (int bin) const{

  return _binEdges.at(bin);

}

double CyclicPhaseBins::getHighBinBoundary (int bin) const{

  double highEdge = 0.0;

  if (bin == ((int)_binEdges.size() - 1) ){
    highEdge = _binEdges.at(0) + TMath::Pi()*2.0;
  }
  else{
    highEdge = _binEdges.at(bin+1);
  }

  return highEdge;
}

double CyclicPhaseBins::getLowBinBoundary (double phase) const{

  int binNum = getBinNumber(phase);
  
  double lowEdge = getLowBinBoundary(binNum);
  
  while (phase < lowEdge){
    lowEdge -= TMath::Pi()*2.0;
  }

  while (phase > lowEdge + TMath::Pi()*2.0){
    lowEdge += TMath::Pi()*2.0;
  }
  
  return lowEdge;
}

double CyclicPhaseBins::getHighBinBoundary(double phase) const{

  int binNum = getBinNumber(phase);

  double highEdge = getHighBinBoundary(binNum);

  while (phase > highEdge){
    highEdge += TMath::Pi()*2.0;
  }

  while (phase < highEdge - TMath::Pi()*2.0){
    highEdge -= TMath::Pi()*2.0;
  }
  
  return highEdge;
}

int CyclicPhaseBins::getNumBins() const{

  return (int)_binEdges.size();

}


CyclicPhaseBins::~CyclicPhaseBins(){

}
