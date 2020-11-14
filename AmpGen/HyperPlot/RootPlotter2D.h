/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class to make plotting 2D root histograms a little easier.
 *
 **/
 
#ifndef ROOT_PLOTTER_2D_HH
#define ROOT_PLOTTER_2D_HH

#include "AmpGen/HyperPlot/RootPlotter.h"

#include "AmpGen/HyperPlot/MessageService.h"

#include "TLine.h"
#include "TArc.h"

class RootPlotter2D : public RootPlotter{

  private:
  
  protected:



  virtual void setHistogramStyle(TH1* histogram, bool setMinMax);

  public:

  virtual double getGlobalMax();
  virtual double getGlobalMin();  

  RootPlotter2D(TH1* histogram, double width = 350, double height = 300);

  void addVerticalLine(double xpos, int style=1);
  void addHorizontalLine(double ypos, int style=1);
  
  void addVerticalBox(double xmin, double xmax, int fillColour, int fillstyle);

  virtual ~RootPlotter2D();
  
};

#endif
