#include "AmpGen/HyperPlot/RootPlotter2D.h"
#include "AmpGen/HyperPlot/GlobalFunctions.h"

/** Construct a RootPlotter with one histogram on a canvas with specified width and height */
RootPlotter2D::RootPlotter2D(TH1* histogram, double width, double height) :
  RootPlotter(histogram, width, height)
{
  WELCOME_LOG << "Hello from the Plotter() constructor!";
  //_yAxisLabelOffset = 0.0;
  _yAxisTitleOffset = 1.3;
  //_xAxisLabelOffset = 0.9;
  _lMargin = 0.20;
  _rMargin = 0.18;
  _tMargin = 0.10;
  _bMargin = 0.14;
}

/** Set the style of a TH2D */
void RootPlotter2D::setHistogramStyle(TH1* histogram, bool setMinMax){

  histogram->SetTitle("");
  
  histogram->GetXaxis()->SetTitleSize(_xAxisTitleSize);
  histogram->GetYaxis()->SetTitleSize(_yAxisTitleSize);

  histogram->GetXaxis()->SetLabelSize(_xAxisLabelSize);
  histogram->GetYaxis()->SetLabelSize(_yAxisLabelSize);

  histogram->GetXaxis()->SetTickLength(_xAxisTickLength);
  histogram->GetYaxis()->SetTickLength(_yAxisTickLength);

  histogram->GetXaxis()->SetLabelOffset(_xAxisLabelOffset);
  histogram->GetYaxis()->SetLabelOffset(_yAxisLabelOffset);

  histogram->GetXaxis()->SetTitleOffset(_xAxisTitleOffset);
  histogram->GetYaxis()->SetTitleOffset(_yAxisTitleOffset);

  histogram->GetXaxis()->SetTitle(_xAxisName);
  histogram->GetYaxis()->SetTitle(_yAxisName);


  if (setMinMax) histogram->SetMaximum(getGlobalMax());
  if (setMinMax) histogram->SetMinimum(getGlobalMin()); 
  if (setMinMax) histogram->GetZaxis()->SetRangeUser(getGlobalMin(), getGlobalMax());

}

/** Add a box that spans the canvas in the vertical direction, and goes from xmin
to xmax in the horizontal direction */
void RootPlotter2D::addVerticalBox(double xmin, double xmax, int fillColour, int fillstyle){
  double ymin = getHistogram(0)->GetYaxis()->GetBinLowEdge(0);
  double ymax = getHistogram(0)->GetYaxis()->GetBinUpEdge( getHistogram(0)->GetYaxis()->GetNbins() );

  TBox* box = new TBox(xmin ,ymin ,xmax, ymax );
  box->SetLineColor(0);
  box->SetFillColor(fillColour);
  box->SetFillStyle(fillstyle);
  _objToPlot.push_back(box);
}

/** find the minimum value accross all histogrmas */
double RootPlotter2D::getGlobalMin(){

  if (_forcedMin != -999.999) return _forcedMin;
  double globalMin = getHistogram(0)->GetBinContent(getHistogram(0)->GetMinimumBin());

  for (unsigned int h=1; h<_histograms.size(); h++)
  {
    double min = getHistogram(h)->GetBinContent(getHistogram(h)->GetMinimumBin());
    if (min < globalMin) globalMin = min;
  }

  if (globalMin > 0) globalMin = globalMin*0.95;
  else globalMin = globalMin*1.05;

  return globalMin;
}

/** find the maximum value accross all histogrmas */
double RootPlotter2D::getGlobalMax(){

  if (_forcedMax != -999.999) return _forcedMax;
  double globalMax = getHistogram(0)->GetBinContent(getHistogram(0)->GetMaximumBin());
  for (unsigned int h=1; h<_histograms.size(); h++)
  {
    double max = getHistogram(h)->GetBinContent(getHistogram(h)->GetMaximumBin());
    if (max > globalMax) globalMax = max;
  }

  globalMax = globalMax*1.05;

  return globalMax;

}

/** Add a line that spans the canvas in the vertical direction, and is at xmin 
in the horizontal direction */
void RootPlotter2D::addVerticalLine(double xpos, int style){
  double ymin = getHistogram(0)->GetYaxis()->GetBinLowEdge(1);
  double ymax = getHistogram(0)->GetYaxis()->GetBinUpEdge (getHistogram(0)->GetXaxis()->GetNbins());
  INFO_LOG << "Adding TLine from (" << xpos << ", " << ymin << ") to (" << xpos << ", " << ymax << ")";
  TLine* line = new TLine(xpos ,ymin ,xpos, ymax );
  line->SetLineColor(1);
  line->SetLineStyle(style);
  _objToPlot.push_back(line);
}

/** Add a line that spans the canvas in the Horizontal direction, and is at ymin 
in the horizontal direction */
void RootPlotter2D::addHorizontalLine(double ypos, int style){
  double xlow  = getHistogram(0)->GetXaxis()->GetBinLowEdge(1);
  double xhigh = getHistogram(0)->GetXaxis()->GetBinUpEdge (getHistogram(0)->GetXaxis()->GetNbins());
  TLine* line = new TLine(xlow ,ypos ,xhigh, ypos );
  line->SetLineColor(1);
  line->SetLineStyle(style);
  _objToPlot.push_back(line);
}

/** destructor */
RootPlotter2D::~RootPlotter2D(){

}
