#include "AmpGen/HyperPlot/RootPlotter.h"

/** Construct a RootPlotter with one histogram on a canvas with specified width and height */
RootPlotter::RootPlotter(TH1* histogram, double width, double height) :
  Plotter("canvas_" + (TString)histogram->GetName(), width, height),
  _xAxisName(histogram->GetXaxis()->GetTitle()),
  _yAxisName(histogram->GetYaxis()->GetTitle())
{
  _histograms.push_back(histogram);

  WELCOME_LOG << "Hello from the Plotter() constructor!";
}


/** Get one of the histograms that has been added to the plotter */
TH1* RootPlotter::getHistogram(int i){
  return (TH1*) _histograms.at(i);
}

/** Get a reference to the draw options of histogram i -  used for getting and setting */
TString& RootPlotter::drawOptions(int i){
  while ((int)_drawOptions.size() <= i) _drawOptions.push_back("");
  return _drawOptions.at(i);
}

/** Get a reference to the draw options of object i -  used for getting and setting */
TString& RootPlotter::objDrawOptions(int i){
  while ((int)_objDrawOptions.size() <= i) _objDrawOptions.push_back("");
  return _objDrawOptions.at(i);
}

/** Add a text string at a point (x,y) on the histogram */
void RootPlotter::addText(TString text, double x, double y, int alignh, int alignv, double size, int ndc, int color){
  
  //1left 2center 3right
  
  TLatex* textL  = new TLatex(x, y, text);

  if (ndc == true) textL->SetNDC();
  textL->SetTextAlign(alignh*10 + alignv);
  textL->SetTextSize(size);
  textL->SetTextColor(color);
  
  this->addObject(textL);

}

/** Add a legend to the Canvas - not really tested */
void RootPlotter::drawLegend(){


  if (_legend != 0) {delete _legend; _legend = 0;}
  
  if      (s_legend_position == "TopRight"   ) _legend = new TLegend(0.55,0.6,0.85,0.78);
  else if (s_legend_position == "BottomRight") _legend = new TLegend(0.55,0.2,0.85,0.38);
  else if (s_legend_position == "TopLeft"    ) _legend = new TLegend(0.15,0.6,0.5,0.78);
  else if (s_legend_position == "BottomLeft" ) _legend = new TLegend(0.6,0.6,0.85,0.8);
  
  VERBOSE_LOG << "Legend position is " << s_legend_position;

  if (_legend != 0){
    _legend->SetFillColor(0);
  
    for (unsigned int i = 1; i < _histograms.size(); i++){
        _legend->AddEntry(getHistogram(i),getHistogram(i)->GetTitle(),"l");
    }
  
    _legend->Draw();
  }

}

/** plot stacked  */
void RootPlotter::plotStacked(TPad* pad, double scaleFactor){


  if (scaleFactor != 1.0){ 
    _tMargin       = 0.0001;
    _bMargin       = _bMargin       * scaleFactor;
    _xAxisLabelSize = _xAxisLabelSize * scaleFactor;
    _yAxisLabelSize = _yAxisLabelSize * scaleFactor;
    _xAxisTitleSize = _xAxisTitleSize * scaleFactor;
    _yAxisTitleSize = _yAxisTitleSize * scaleFactor;
  }
  this->setCanvasDefaults(pad);
  gStyle->SetHistTopMargin(0.0);

  THStack* histogramStack = new THStack("stackedHistograms","stackedHistograms");


  for (unsigned int i = 0; i < _histograms.size(); i++) {
    this->setHistogramStyle( getHistogram(i) , 0);
    getHistogram(i)->SetFillColor(getColor(i));
    histogramStack->Add(getHistogram(i));
  }

  _canvas->cd();
  histogramStack->Draw("hist");
  pad->Draw(); 

  histogramStack->GetYaxis()->SetNdivisions( getHistogram(0)->GetYaxis()->GetNdivisions () );
  histogramStack->GetYaxis()->CenterTitle  ( getHistogram(0)->GetYaxis()->GetCenterTitle() );
  histogramStack->GetXaxis()->SetNdivisions( getHistogram(0)->GetXaxis()->GetNdivisions () );
  histogramStack->GetXaxis()->CenterTitle  ( getHistogram(0)->GetXaxis()->GetCenterTitle() );

  //The histogram is not drawn until the pad is drawn on a canvas - we must therefore 
  // draw the pad before the GetHistogram() returns something. Might run into problems
  // if we start using pads within pads

  this->setHistogramStyle( histogramStack->GetHistogram(), 0 );

  if (scaleFactor != 1.0){   
    //histogramStack->GetHistogram()->GetXaxis()->SetTickLength(histogramStack->GetHistogram()->GetXaxis()->GetTickLength()*scaleFactor);
    histogramStack->GetHistogram()->GetYaxis()->SetNdivisions(502);
  }

  pad->cd();
  histogramStack->DrawClone("hist");

  delete histogramStack;
}


/** plot with the option SAME  */
void RootPlotter::plotSame(TPad* pad, TString plotOptions, double scaleFactor){
  
  if (pad == 0) pad = _canvas;
  this->setCanvasDefaults(pad);
  pad->cd();
  
  if (scaleFactor != 1.0){ 
    //std::cout << "Scaling histogram in pad " << pad->GetName() << " by " << scaleFactor << std::endl;
    _yAxisLabelSize   = _yAxisLabelSize   * scaleFactor;
    _yAxisTitleSize   = _yAxisTitleSize   * scaleFactor;
    _xAxisLabelOffset = _xAxisLabelOffset * scaleFactor;
    _xAxisTickLength  = _xAxisTickLength  * scaleFactor;
    _yAxisTitleOffset = _yAxisTitleOffset * (1.0/scaleFactor);
    //getHistogram(0)->GetXaxis()->SetTickLength( getHistogram(0)->GetXaxis()->GetTickLength()*scaleFactor );
    //getHistogram(0)->GetYaxis()->SetTickLength( getHistogram(0)->GetYaxis()->GetTickLength()*scaleFactor );
    //getHistogram(0)->GetYaxis()->SetNdivisions(502);
  }
  

  
  VERBOSE_LOG << "Setting histogram style";
  this->setHistogramStyle( getHistogram(0) );
  if (_usePresetColours) getHistogram(0)->SetLineColor(getColor(0));
  if (_usePresetColours) getHistogram(0)->SetMarkerColor(getColor(0));
  VERBOSE_LOG << "Draw axis";

  getHistogram(0)->DrawCopy(plotOptions + "AXIS");
  this->setHistogramStyle( getHistogram(0) );
  getHistogram(0)->DrawCopy(plotOptions + "AXIS");


  for (unsigned int i = 0; i < _histograms.size(); i++){
    VERBOSE_LOG << "Set style" << i;
    this->setHistogramStyle( getHistogram(i) );   
    if (_usePresetColours)getHistogram(i)->SetLineColor(getColor(i));
    //getHistogram(i)->SetFillColor(getColor(i));
    if (_usePresetColours)getHistogram(i)->SetMarkerColor(getColor(i));
    VERBOSE_LOG << "Draw" << i;
    getHistogram(i)->DrawCopy(plotOptions + drawOptions(i) + "SAME");
  }



  for (unsigned int i = 0; i < _objToPlot.size(); i++){
    if ( TString(_objToPlot.at(i)->ClassName()).Contains("TGraph") ){
      if (_usePresetColours) ((TGraph*)_objToPlot.at(i))->SetLineColor  (getColor(i));
      if (_usePresetColours) ((TGraph*)_objToPlot.at(i))->SetMarkerColor(getColor(i));      
    }
    _objToPlot.at(i)->DrawClone(plotOptions + "SAME" + objDrawOptions(i) );
    //INFO_LOG << "Drawing object" << i << " with option " << objDrawOptions(i);
  }
  VERBOSE_LOG << "Draw axis";
  getHistogram(0)->DrawCopy(plotOptions + "SAME AXIS");

  drawLegend();

}

/** Call this to save the canvas to a specific directory. */
void RootPlotter::plot(TString plotDirectory, TString plotOptions, TPad* pad, double scaleFactor){
  
  if (pad == 0) pad = _canvas;

  this->setCanvasDefaults(pad);
  pad->cd();
  //set the histogram options on the first histogram   

  if (plotOptions == "STACKED")   plotStacked(pad, scaleFactor);
  else                            plotSame   (pad, plotOptions, scaleFactor);

  pad->Update();
  if (plotDirectory != "") {
    pad->Print(plotDirectory + s_imageformat);
    if (s_imageformat2 != "") pad->Print(plotDirectory + s_imageformat2);
    if (_allImageFormats == true){
      pad->Print(plotDirectory + ".pdf");
      pad->Print(plotDirectory + ".eps");
      pad->Print(plotDirectory + ".png");
      pad->Print(plotDirectory + ".C"  );
    }
  }
}


/** Destructor. */
RootPlotter::~RootPlotter(){

}
