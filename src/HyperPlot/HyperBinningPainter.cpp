#include "AmpGen/HyperPlot/HyperBinningPainter.h"

/** Construct a HyperBinningPainter for a given HyperBinningHistogram */
HyperBinningPainter::HyperBinningPainter(HyperHistogram* histogram) :
  _histogram(histogram),
  _density(false)
{
  WELCOME_LOG << "Good day from the HyperBinningPainter() Constructor";  
}

/** Draw the HyperBinningHistogram */
void HyperBinningPainter::draw(TString path, TString option){
  
  option = option; //remove conpilation warning


  int nBins = _histogram->getNBins();
  
  //double volumeSum = 0.0;
  //for (int i = 0; i < nBins; i++){
  //  double volume = _histogram->getBinning().getBinHyperVolume(i).volume();
  //  volumeSum += volume;
  //}
//
  //double* edgeArray = new double [nBins + 1];
  //edgeArray[0]=0.0;
//
  //double binSum = 0.0;
//
  //for (int i = 0; i < nBins; i++){
  //  double volume = _histogram->getBinning().getBinHyperVolume(i).volume();
  //  binSum += 1.0;//volume/volumeSum;
  //  edgeArray[i+1] = binSum;
  //}
//
  TH1D tempHist("tempHist", "tempHist", nBins, 0.0, nBins);
  tempHist.GetXaxis()->SetTitle("Bin Number");
  //tempHist.GetYaxis()->SetTitle("Frequency");

  //int dim = _histogram->getBinning().getDimension();
  
  //TCanvas canvas("canvas", "canvas", 400.0, dim*80.0 + 2.0*30.0 + 300.0);



  //TPad padPlot   ("padPlot", "padPlot"      , 0.0, (dim*80.0 + 2.0*30.0)/(dim*80.0 + 2.0*30.0 + 300.0) , 1.0, 1.0);
  //canvas.cd();
  //padPlot.Draw();
//
  //TPad padBinning("padBinning", "padBinning", 0.0, 0.0 , 1.0, (dim*80.0 + 2.0*30.0)/(dim*80.0 + 2.0*30.0 + 300.0) );
  //canvas.cd();
  //padBinning.Draw();

  for (int i = 0; i < nBins; i++){
    if (_density == true) {
      double volume = _histogram->getBinning().getBinHyperVolume(i).volume();
      tempHist.SetBinContent(i+1, _histogram->getBinContent(i)/volume);
      tempHist.SetBinError  (i+1, _histogram->getBinError(i)/volume);
    }
    else{
      tempHist.SetBinContent(i+1, _histogram->getBinContent(i));
      tempHist.SetBinError  (i+1, _histogram->getBinError(i));
    }
  }

  RootPlotter1D plotter(&tempHist);
  plotter.plot( path );

 // plotter.plot("", "", &padPlot);
 // _histogram->getBinning().drawBinning("", &padBinning);
 // 
//
 // canvas.Draw();
//
 // gStyle->SetPaperSize(400.0, dim*80.0 + 2.0*30.0 + 300.0);
//
 // padBinning.Print(path + "_bin" +  Plotter::s_imageformat);
 // padPlot.Print(path + "_data" +  Plotter::s_imageformat);
//
 // canvas.Print(path + Plotter::s_imageformat);

  //INFO_LOG << "I can't plot N dimensional histograms just yet";

}
/** Destructor */
HyperBinningPainter::~HyperBinningPainter(){

}


