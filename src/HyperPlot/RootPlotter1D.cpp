#include "AmpGen/HyperPlot/RootPlotter1D.h"
#include "AmpGen/HyperPlot/GlobalFunctions.h"

double RootPlotter1D::s_ratioMax = -999.999;
double RootPlotter1D::s_ratioMin = -999.999;

/** Construct a RootPlotter with one histogram on a canvas with specified width and height */
RootPlotter1D::RootPlotter1D(TH1* histogram, double width, double height) :
  RootPlotter(histogram, width, height),
  _ratioMax(s_ratioMax),
  _ratioMin(s_ratioMin)
{

  _lMargin = 0.14;
  _rMargin = 0.10;
  _tMargin = 0.10;
  _bMargin = 0.14;
  WELCOME_LOG << "Hello from the Plotter() constructor!";
}


/** Set the style of a TH1D */
void RootPlotter1D::setHistogramStyle(TH1* histogram, bool setMinMax){
     
  VERBOSE_LOG << histogram << "  " << histogram->GetName();

  histogram->SetTitle("");

  histogram->GetXaxis()->SetTitleSize(_xAxisTitleSize);
  histogram->GetYaxis()->SetTitleSize(_yAxisTitleSize);
  
  if ( histogram->GetXaxis()->GetLabelSize() != 0.0 ) 
    histogram->GetXaxis()->SetLabelSize(_xAxisLabelSize);
  if ( histogram->GetYaxis()->GetLabelSize() != 0.0 ) 
    histogram->GetYaxis()->SetLabelSize(_yAxisLabelSize);
  

  histogram->GetXaxis()->SetTickLength(_xAxisTickLength);
  histogram->GetXaxis()->SetTickLength(_yAxisTickLength);
  
  if ( histogram->GetXaxis()->GetLabelOffset() != 999 ) 
    histogram->GetXaxis()->SetLabelOffset(_xAxisLabelOffset);
  if ( histogram->GetYaxis()->GetLabelOffset() != 999 ) 
    histogram->GetYaxis()->SetLabelOffset(_yAxisLabelOffset);
  

  histogram->GetXaxis()->SetTitleOffset(_xAxisTitleOffset);
  histogram->GetYaxis()->SetTitleOffset(_yAxisTitleOffset);

  histogram->GetXaxis()->SetTitle(_xAxisName);
  histogram->GetYaxis()->SetTitle(_yAxisName);

  double nBins = (double)histogram->GetNbinsX();
  double markerSize = 15.0/nBins;
  if (markerSize > 0.3) markerSize = 0.3;
  histogram->SetMarkerSize(markerSize);

  if (setMinMax) histogram->SetMaximum(getGlobalMax());
  if (setMinMax) histogram->SetMinimum(getGlobalMin());

  if (_bMargin < 0.01) histogram->GetXaxis()->SetNoExponent();

  VERBOSE_LOG << "Global Min = " << getGlobalMin();
  VERBOSE_LOG << "Global Max = " << getGlobalMax();

}


/** Add a line that spans the canvas in the vertical direction, and is at xmin 
in the horizontal direction */
void RootPlotter1D::addVerticalLine(double xpos, int style, int colour){
  double ymin = getGlobalMin();
  double ymax = getGlobalMax();
  VERBOSE_LOG << "Adding TLine from (" << xpos << ", " << ymin << ") to (" << xpos << ", " << ymax << ")";
  TLine* line = new TLine(xpos ,ymin ,xpos, ymax );
  line->SetLineColor(colour);
  line->SetLineStyle(style);
  _objToPlot.push_back(line);
}

/** Add a line that spans the canvas in the horizontal direction, and is at xmin 
in the vertical direction */
void RootPlotter1D::addHorizontalLine(double ypos, int style, int colour){
  double xlow  = getHistogram(0)->GetXaxis()->GetBinLowEdge(1);
  double xhigh = getHistogram(0)->GetXaxis()->GetBinUpEdge (getHistogram(0)->GetXaxis()->GetNbins());
  TLine* line = new TLine(xlow ,ypos ,xhigh, ypos );
  line->SetLineColor(colour);
  line->SetLineStyle(style);
  _objToPlot.push_back(line);
}

/** Add a box that spans the canvas in the horizontal direction, and goes from xmin
to xmax in the vertical direction */
void RootPlotter1D::addHorizontalBox(double ypos, double width, int fillColour){

  double xlow  = getHistogram(0)->GetXaxis()->GetBinLowEdge(1);
  double xhigh = getHistogram(0)->GetXaxis()->GetBinUpEdge (getHistogram(0)->GetXaxis()->GetNbins());
  TBox* box = new TBox(xlow ,ypos - (width*0.5) ,xhigh, ypos + (width*0.5) );
  box->SetLineColor(0);
  box->SetFillColor(fillColour);
  _objToPlot.push_back(box);
}

/** Add a box that spans the canvas in the vertical direction, and goes from xmin
to xmax in the horizontal direction */
void RootPlotter1D::addVerticalBox(double xmin, double xmax, int fillColour, int fillstyle){
  double ymin = getGlobalMin();
  double ymax = getGlobalMax();

  TBox* box = new TBox(xmin ,ymin ,xmax, ymax );
  box->SetLineColor(0);
  box->SetFillColor(fillColour);
  box->SetFillStyle(fillstyle);
  _objToPlot.push_back(box);
}

/** find the minimum value accross all histogrmas */
double RootPlotter1D::getGlobalMin(){

  if ( _forcedMin != -999.999) return  _forcedMin;

  double globalMin = getHistogram(0)->GetBinContent(getHistogram(0)->GetMinimumBin());

  for (unsigned int h=1; h<_histograms.size(); h++)
  {
    double min = getHistogram(h)->GetBinContent(getHistogram(h)->GetMinimumBin());
    if (min < globalMin) globalMin = min;
  }

  if (globalMin > 0) globalMin = globalMin*0.95;
  else               globalMin = globalMin*1.05;

  return globalMin;
}

/** find the maximum value accross all histogrmas */
double RootPlotter1D::getGlobalMax(){

  if ( _forcedMax != -999.999) return  _forcedMax;

  double globalMax = getHistogram(0)->GetBinContent(getHistogram(0)->GetMaximumBin());
  for (unsigned int h=1; h<_histograms.size(); h++)
  {
    double max = getHistogram(h)->GetBinContent(getHistogram(h)->GetMaximumBin());
    if (max > globalMax) globalMax = max;
  }

  if (globalMax > 0) globalMax = globalMax*1.05;
  else               globalMax = globalMax*0.95;



  return globalMax;

}



RootPlotter1D::~RootPlotter1D(){

}







void RootPlotter1D::plotWithRatio(TString plotDirectory, TString plotOptions, TPad* pad){

  double splitHeight = 0.3;
  double scaleFactor = (1.0 - splitHeight) / splitHeight;

  if (pad == 0) pad = _canvas;  

  TString upperPadName = (TString)_canvas->GetName() + "_upper";
  TString lowerPadName = (TString)_canvas->GetName() + "_lower";
  TPad* upperPad = new TPad(upperPadName, upperPadName, 0.0, splitHeight, 1.0, 1.0);
  TPad* lowerPad = new TPad(lowerPadName, lowerPadName, 0.0, 0.0, 1.0, splitHeight);

  plotRatio("", "",lowerPad, scaleFactor);
  plot     ("", plotOptions,upperPad );
  
  pad->cd(); 
  upperPad->Draw();
  lowerPad->Draw();
  if (plotDirectory != "") pad->Print(plotDirectory + s_imageformat);  
}  


void RootPlotter1D::plotRatio(TString plotDirectory, TString plotOptions, TPad* pad, double scaleFactor, double* returnMin, double* returnMax){

  if (pad == 0) pad = _canvas;
  pad->cd();
  TH1D* divHist0 = new TH1D( gDivideTH1D((TH1D*)getHistogram(1), (TH1D*)getHistogram(0), (TString)"division_1" + getHistogram(0)->GetName() ));
  divHist0->GetXaxis()->SetTitle( getHistogram(0)->GetXaxis()->GetTitle() );
  divHist0->GetYaxis()->SetTitle( getHistogram(0)->GetYaxis()->GetTitle() );

  //if there is a scale factor it means that be are adding this underneath another plot.
  //We therefore adjust the histogram ticks to be the same in both

  if (scaleFactor != 1.0){   
    divHist0->GetXaxis()->SetTickLength(divHist0->GetXaxis()->GetTickLength()*scaleFactor);
    divHist0->GetYaxis()->SetNdivisions(502);
  }

  RootPlotter1D divPlotter(divHist0);
  divPlotter.setHistogramOwnership(true);

  for (unsigned int i = 2; i < _histograms.size(); i++) {
    TH1D* divHist = new TH1D( gDivideTH1D((TH1D*)getHistogram(i), (TH1D*)getHistogram(0), "division_" + makeString<int>(i) + getHistogram(i)->GetName()) );
    divPlotter.add(divHist);
  }

  //if there is a scale factor it means that be are adding this underneath another plot.
  //We therefore adjust all the margins and scale the labels to be the same in both sets
  //of plots

  if (scaleFactor != 1.0){ 
    this->     _bMargin = 0.015;
    divPlotter._tMargin = 0.07;
    divPlotter._bMargin       = divPlotter._bMargin       * scaleFactor;
    divPlotter._xAxisLabelSize = divPlotter._xAxisLabelSize * scaleFactor;
    divPlotter._yAxisLabelSize = divPlotter._yAxisLabelSize * scaleFactor;
    divPlotter._xAxisTitleSize = divPlotter._xAxisTitleSize * scaleFactor;
    divPlotter._yAxisTitleSize = divPlotter._yAxisTitleSize * scaleFactor;
  }

  double min = divPlotter.getRatioMin(&divPlotter);
  double max = divPlotter.getRatioMax(&divPlotter);
  
  if ( fabs(min-1) > fabs(max-1) ) max = 1 + fabs(min-1);

  divPlotter.setMin(2.0-max);
  divPlotter.setMax(max);

  if(returnMin != 0) *returnMin = 2.0-max;
  if(returnMax != 0) *returnMax = max;

  divPlotter.addHorizontalLine(1.0, 2);
  divPlotter.plot("", plotOptions, pad);
  

  if (plotDirectory != "") pad->Print(plotDirectory + s_imageformat);
}

void RootPlotter1D::plotPulls(TString plotDirectory, TString plotOptions , TPad* pad, double scaleFactor){

  if (pad == 0) pad = _canvas;
  pad->cd();
  TH1D* divHist0 = new TH1D( gPullTH1D((TH1D*)getHistogram(1), (TH1D*)getHistogram(0), (TString)"division_1" + getHistogram(0)->GetName() ));
  divHist0->GetXaxis()->SetTitle( getHistogram(0)->GetXaxis()->GetTitle() );
  divHist0->GetYaxis()->SetTitle( getHistogram(0)->GetYaxis()->GetTitle() );

  //if there is a scale factor it means that be are adding this underneath another plot.
  //We therefore adjust the histogram ticks to be the same in both

  if (scaleFactor != 1.0){   
    divHist0->GetXaxis()->SetTickLength(divHist0->GetXaxis()->GetTickLength()*scaleFactor);
    divHist0->GetYaxis()->SetNdivisions(502);
  }

  RootPlotter1D divPlotter(divHist0);
  for (unsigned int i = 2; i < _histograms.size(); i++) {
    TH1D* divHist = new TH1D( gPullTH1D((TH1D*)getHistogram(i), (TH1D*)getHistogram(0), "division_" + makeString<int>(i) + getHistogram(i)->GetName()) );
    divPlotter.add(divHist);
  }

  //if there is a scale factor it means that be are adding this underneath another plot.
  //We therefore adjust all the margins and scale the labels to be the same in both sets
  //of plots

  if (scaleFactor != 1.0){ 
    this->     _bMargin = 0.0001;
    divPlotter._tMargin = 0.0001;
    divPlotter._bMargin       = divPlotter._bMargin       * scaleFactor;
    divPlotter._xAxisLabelSize = divPlotter._xAxisLabelSize * scaleFactor;
    divPlotter._yAxisLabelSize = divPlotter._yAxisLabelSize * scaleFactor;    
    divPlotter._xAxisTitleSize = divPlotter._xAxisTitleSize * scaleFactor;
    divPlotter._yAxisTitleSize = divPlotter._yAxisTitleSize * scaleFactor;    
    divPlotter.setMax(5.0);
    divPlotter.setMin(-5.0);
  }



  divPlotter.addHorizontalBox(0.0, 2.0, 17);  
  divPlotter.addHorizontalLine(0.0, 2);
  divPlotter.plot("", plotOptions, pad);
  

  if (plotDirectory != "") pad->Print(plotDirectory + s_imageformat);

}

void RootPlotter1D::plotWithPulls(TString plotDirectory, TString plotOptions , TPad* pad){
  double splitHeight = 0.3;
  double scaleFactor = (1.0 - splitHeight) / splitHeight;

  if (pad == 0) pad = _canvas;  

  TString upperPadName = (TString)_canvas->GetName() + "_upper";
  TString lowerPadName = (TString)_canvas->GetName() + "_lower";
  TPad* upperPad = new TPad(upperPadName, upperPadName, 0.0, splitHeight, 1.0, 1.0);
  TPad* lowerPad = new TPad(lowerPadName, lowerPadName, 0.0, 0.0, 1.0, splitHeight);

  plotPulls("", plotOptions,lowerPad, scaleFactor);
  plot     ("", plotOptions,upperPad);
  
  pad->cd(); 
  upperPad->Draw();
  lowerPad->Draw();
  if (plotDirectory != "") pad->Print(plotDirectory + s_imageformat);  
}

double RootPlotter1D::getRatioMin(RootPlotter1D* ratioPlotter){

  if (  _ratioMin != -999.999) return  _ratioMin;
  if ( s_ratioMin != -999.999) return s_ratioMin;

  ratioPlotter->setMin(-999.999);

  double min = ratioPlotter->getGlobalMin() - 1;
  double absMin = fabs(min);

  double base = pow(10, floor(log10(absMin)));
  int unit = ceil(double(absMin)/base);
  double roundedMin = double(unit)*base;

  if (min < 0.0) roundedMin = -roundedMin;

  return 1.0 + roundedMin;

}

double RootPlotter1D::getRatioMax(RootPlotter1D* ratioPlotter){

  if (  _ratioMax != -999.999) return  _ratioMax;
  if ( s_ratioMax != -999.999) return s_ratioMax;
   
  ratioPlotter->setMax(-999.999);

  double max = ratioPlotter->getGlobalMax() - 1;
  double absMax = fabs(max);

  double base = pow(10, floor(log10(absMax)));
  int unit = ceil(double(absMax)/base);
  double roundedMax = double(unit)*base;

  if (max < 0.0) roundedMax = -roundedMax;

  return 1.0 + roundedMax;

}
