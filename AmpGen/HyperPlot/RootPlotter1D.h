/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class to make plotting 1D root histograms a little easier.
 *
 **/
 
#ifndef ROOT_PLOTTER_1D_HH
#define ROOT_PLOTTER_1D_HH

#include "AmpGen/HyperPlot/RootPlotter.h"


#include "AmpGen/HyperPlot/MessageService.h"

#include "TLine.h"

class RootPlotter1D : public RootPlotter{

  private:
 
  protected:
  

  double getRatioMin(RootPlotter1D* ratioPlotter); /**< \todo this stuff never really worked */
  double getRatioMax(RootPlotter1D* ratioPlotter); /**< \todo this stuff never really worked */
  double _ratioMax;                                /**< \todo this stuff never really worked */
  double _ratioMin;                                /**< \todo this stuff never really worked */

  virtual void setHistogramStyle(TH1* histogram, bool setMinMax = 1);


  public:
  
  virtual double getGlobalMin();
  virtual double getGlobalMax();

  RootPlotter1D(TH1* histogram, double width = 400, double height = 300);

  void addVerticalLine(double xpos  , int style=1, int colour=1);
  void addHorizontalLine(double ypos, int style=1, int colour=1);
  void addHorizontalBox(double ypos, double width, int fillColour);
  void addVerticalBox(double xmin, double xmax, int fillColour, int fillstyle);

  virtual ~RootPlotter1D();


  static double s_ratioMax;
  /**< \todo this stuff never really worked */
 
  static double s_ratioMin;
  /**< \todo this stuff never really worked */

  void setRatioMax(double val){_ratioMax = val;}
  /**< \todo this stuff never really worked */

  void setRatioMin(double val){_ratioMin = val;}
  /**< \todo this stuff never really worked */

  void plotRatio(TString plotDirectory, TString plotOptions = "", TPad* pad = 0, double scaleFactor = 1.0, double* returnMin = 0, double* returnMax = 0);
  /**< \todo this stuff never really worked */
  void plotPulls(TString plotDirectory, TString plotOptions = "", TPad* pad = 0, double scaleFactor = 1.0);
  /**< \todo this stuff never really worked */
  void plotWithRatio(TString plotDirectory, TString plotOptions = "", TPad* pad = 0);
  /**< \todo this stuff never really worked */
  void plotWithPulls(TString plotDirectory, TString plotOptions = "", TPad* pad = 0);
  /**< \todo this stuff never really worked */


  
};

#endif
