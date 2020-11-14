#include "AmpGen/HyperPlot/Plotter.h"

TString Plotter::s_imageformat     = ".pdf";
TString Plotter::s_imageformat2    = "";
TString Plotter::s_legend_position = "RightTop";
int     Plotter::s_plotterCount    = 0;
double  Plotter::s_forcedMax       = -999.999;
double  Plotter::s_forcedMin       = -999.999;

/** 
Set the image format for all plots to be made
*/
void Plotter::setImageFormat(TString format){
  s_imageformat = format;
}

/** 
Contruct a plotter with some name, width and height
*/
Plotter::Plotter(TString canvasName, double width, double height) :
  _canvas(0),
  _legend(0),
  _forcedMax(s_forcedMax),
  _forcedMin(s_forcedMin),
  _lMargin(0.11),
  _rMargin(0.15),
  _tMargin(0.11),
  _bMargin(0.15),
  _xAxisTitleOffset(1.0),
  _yAxisTitleOffset(1.0),
  _xAxisLabelOffset(0.010),
  _yAxisLabelOffset(0.005),
  _xAxisTickLength(0.035),
  _yAxisTickLength(0.015),
  _xAxisLabelSize(0.06),
  _yAxisLabelSize(0.06),
  _xAxisTitleSize(0.06),
  _yAxisTitleSize(0.06),
  _histogramOwnership(0),
  _objectOwnership(1),
  _usePresetColours(1),
  _allImageFormats (0)
{
  WELCOME_LOG << "Hello from the Plotter() constructor!";
  TString canvasNameUnique = canvasName; canvasNameUnique += s_plotterCount;
  _canvas = new TCanvas(canvasNameUnique, canvasName, width, height);
  _colours.push_back(1);
  _colours.push_back(kRed);
  _colours.push_back(kMagenta);
  _colours.push_back(kGreen);
  _colours.push_back(kCyan);
  _colours.push_back(kYellow);
  _colours.push_back(kBlue);
  s_plotterCount++;
}

/**
Copy another plotter - this will hold the same pointers, 
so make sure that the new plotter doesn't own any of them
*/
Plotter::Plotter(const Plotter& other) :
  _canvas( new TCanvas( other._canvas )),
  _legend(other._legend),
  _forcedMax(other._forcedMax),
  _forcedMin(other._forcedMin),
  _objToPlot(other._objToPlot),
  _histograms(other._histograms),
  _colours(other._colours),
  _lMargin(other._lMargin),
  _rMargin(other._rMargin),
  _tMargin(other._tMargin),
  _bMargin(other._bMargin),
  _xAxisLabelOffset(other._xAxisLabelOffset),
  _yAxisLabelOffset(other._yAxisLabelOffset),
  _histogramOwnership(0),
  _objectOwnership(0),
  _usePresetColours(1),
  _allImageFormats (0)
{

}

/** Set defaults on the TCanvas */
void Plotter::setCanvasDefaults(TPad* pad){

  pad->SetLeftMargin(_lMargin);
  pad->SetBottomMargin(_bMargin);
  pad->SetRightMargin(_rMargin);
  pad->SetTopMargin(_tMargin);
  
}

void Plotter::setColor(int i, int color){
  
  while (i >= (int)_colours.size()) {
    _colours.push_back( getColor(i) );
  }
  _colours.at(i) = color;

}

/** how many objects have been added to the plotter */
int Plotter::getNumObjects(){

  return _objToPlot.size();

}

int Plotter::getColor(int i){

  VERBOSE_LOG << "Colour " << i;

  if (i >= (int)_colours.size()) {
    int j = i%_colours.size();
    int k = (int)(i/_colours.size());
    return _colours.at(j) + 2*k;
  }
  return _colours.at(i);
}

/** add a histogram to the plotter */
void Plotter::add(TObject* histogram){
  _histograms.push_back(histogram);
}

/** add a dot to some point in the canvas

The TString 'shape' can be

  - square (draws a TBox     at that (x,y) coord )
  - circle (draws a TEllipse at that (x,y) coord )

 */
void Plotter::addDot(double xpos, double ypos, double size, int colour, TString shape, double sizeY){
  
  if (sizeY == 0.0) sizeY = size; 

  if (shape == "square") {
    TBox* box = new TBox(xpos-size*0.5,ypos-sizeY*0.5,xpos+size*0.5,ypos+sizeY*0.5);
    box->SetFillColor(colour);
    _objToPlot.push_back(box);
  }
  if (shape == "circle") {
    TEllipse* box = new TEllipse(xpos,ypos,size*0.5,sizeY*0.5);
    box->SetFillColor(colour);
    _objToPlot.push_back(box);
  }
  
}

/** use a log scale for the x-axis */
void Plotter::logX(bool log){
  _canvas->SetLogx(log);
}

/** use a log scale for the y-axis */
void Plotter::logY(bool log){
  _canvas->SetLogy(log);
}

/** use a log scale for the z-axis */
void Plotter::logZ(bool log){
  _canvas->SetLogz(log);   
}

/** Add an object to the Plotter (e.g. TMarker, TBox...).
By default, the plotter takes ownership of the object */
void Plotter::addObject(TObject* obj){
  _objToPlot.push_back(obj);
}

void Plotter::setPropertiesFromTH1(TH1* hist){
  
  _xAxisTitleOffset = hist->GetXaxis()->GetTitleOffset(); //(1.0),
  _yAxisTitleOffset = hist->GetYaxis()->GetTitleOffset(); //(1.0),
  _xAxisLabelOffset = hist->GetXaxis()->GetLabelOffset(); //(0.010),
  _yAxisLabelOffset = hist->GetYaxis()->GetLabelOffset(); //(0.005),
  _xAxisTickLength  = hist->GetXaxis()->GetTickLength(); //(0.035),
  _yAxisTickLength  = hist->GetYaxis()->GetTickLength(); //(0.015),
  _xAxisLabelSize   = hist->GetXaxis()->GetLabelSize(); //(0.06),
  _yAxisLabelSize   = hist->GetYaxis()->GetLabelSize(); //(0.06),
  _xAxisTitleSize   = hist->GetXaxis()->GetTitleSize(); //(0.06),
  _yAxisTitleSize   = hist->GetYaxis()->GetTitleSize(); //(0.06),

}

/**<
 Destructor - delete all the histograms and objects that I own.
*/
Plotter::~Plotter(){
  if (_objectOwnership == 1){
    for(unsigned int i = 0; i < _objToPlot.size(); i++) delete _objToPlot.at(i);
  }
  if (_histogramOwnership == 1){
    for(unsigned int i = 0; i < _histograms.size(); i++) delete _histograms.at(i);
  }
  delete _canvas;
  if (_legend != 0) delete _legend;
}

