#include "AmpGen/HyperPlot/HyperBinningPainter1D.h"

/** Construct a 1D HyperBinningPainter for a given HyperBinningHistogram */
HyperBinningPainter1D::HyperBinningPainter1D(HyperHistogram* histogram) :
  HyperBinningPainter(histogram)
{
  WELCOME_LOG << "Good day from the HyperBinningPainter1D() Constructor";  

}

/** Turn the HyperBinningHistogram into a TH1D */
TH1D* HyperBinningPainter1D::getHistogram(TString histname){

  int nBins = _histogram->getBinning().getNumBins();

  std::vector<double> binEdges;
  for (int i = 0; i < nBins; i++){
    double min = _histogram->getBinning().getBinHyperVolume(i).getMin(0);
    binEdges.push_back(min);
  }
  binEdges.push_back( _histogram->getBinning().getMax(0) );

  std::sort(binEdges.begin(), binEdges.end());

  double* aryBinEdges = new double [binEdges.size()];
  for (unsigned i = 0; i < binEdges.size(); i++) {
    VERBOSE_LOG << "Bin edge " << binEdges.at(i);
    aryBinEdges[i] = binEdges.at(i);
  }

  TH1D* tempHist = new TH1D(histname, histname, nBins, aryBinEdges);
  tempHist->GetXaxis()->SetTitle("Val");
  tempHist->GetYaxis()->SetTitle("Frequency");

  for (int i = 0; i < nBins; i++){
    int bin = _histogram->getBinning().getBinNum( HyperPoint( tempHist->GetXaxis()->GetBinCenter(i+1) ) );
  
    if (_density == true) {
      double volume = _histogram->getBinning().getBinHyperVolume(bin).volume();
      tempHist->SetBinContent(i+1, _histogram->getBinContent(bin)/volume);
      tempHist->SetBinError  (i+1, _histogram->getBinError(bin  )/volume);
    }
    else{
      double val = _histogram->getBinContent(bin);
      double err = _histogram->getBinError  (bin);
      tempHist->SetBinContent(i+1, val);
      tempHist->SetBinError  (i+1, err);
      std::cout << "Bin " << i << ":   Content = " << val << "    Error = " << err << std::endl;
    }
  }
  
  return tempHist;

}

/** Draw the HyperBinningHistogram  */
void HyperBinningPainter1D::draw(TString path, TString option){
  
  option = option; //remove conpilation warning

  TH1D* hist = getHistogram("tempToDraw");

  RootPlotter1D plotter(hist);
  plotter.plot(path, "E");

  delete hist;

}

/** Destructor  */
HyperBinningPainter1D::~HyperBinningPainter1D(){

}



