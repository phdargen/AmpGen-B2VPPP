/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class to make plotting root histograms a little easier.
 *
 **/

 
#ifndef ROOT_PLOTTER_HH
#define ROOT_PLOTTER_HH

#include "TH1.h"
#include "THStack.h"
#include "Plotter.h"
#include "TASImage.h"
#include "TSystem.h"
#include "AmpGen/HyperPlot/MessageService.h"
#include "TLatex.h"
#include "TGraph.h"

class RootPlotter : public Plotter{

  protected:

  TString          _xAxisName;     /**< xAxis name - by default taken from the first histogram added */
  TString          _yAxisName;     /**< yAxis name - by default taken from the first histogram added */
  
  std::vector<TString> _drawOptions;    /**< options passed to Root on how to draw each histogram added */
  std::vector<TString> _objDrawOptions; /**< options passed to Root on how to draw each object added */
  
  /** purely virtual function that sets the correct histogram style (different for TH1 and TH2) */
  virtual void setHistogramStyle(TH1* histogram, bool setMinMax = 1) = 0;


  TH1* getHistogram(int i);

  public:
  
  void setXaxisName(TString name){_xAxisName = name;} /**< Set the x-axis title name */
  void setYaxisName(TString name){_yAxisName = name;} /**< Set the y-axis title name */

  virtual double getGlobalMax(){return 0.0;} /**< Get the minimum of all histograms added (or the forced min if set) */
  virtual double getGlobalMin(){return 0.0;} /**< Get the maximum of all histograms added (or the forced max if set) */

  RootPlotter(TH1* histogram, double width = 300, double height = 200);

  TString& drawOptions(int i);
  TString& objDrawOptions(int i);
  
  //void plotDivisons(TPad* pad);
  void plotStacked(TPad* pad, double scaleFactor);
  void plotSame   (TPad* pad, TString plotOptions, double scaleFactor = 1.0);
  virtual void plot(TString plotDirectory, TString plotOptions = "", TPad* pad = 0, double scaleFactor = 1.0);
  
  void addText(TString text, double x, double y, int alignh = 1, int alignv = 2, double size = 0.06, int ndc = true, int color = kBlack);
  
  void drawLegend();

  virtual ~RootPlotter();
  
};

#endif
