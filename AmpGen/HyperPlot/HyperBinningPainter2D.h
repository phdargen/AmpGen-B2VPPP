/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class plotting 2D HyperBinningHistograms 
 *
 **/

 
#ifndef HYPERBINNINGPAINTER2D_HH
#define HYPERBINNINGPAINTER2D_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/RootPlotter1D.h"
#include "AmpGen/HyperPlot/RootPlotter2D.h"
#include "AmpGen/HyperPlot/HyperBinningPainter.h"


// Root includes
#include "TH1D.h"
#include "TH2D.h"

// std includes



class HyperBinningPainter2D : public HyperBinningPainter {

  private:
  
  BinningBase* _binning; 
  /**< This gets filled in an alternate constructor (usually the HyperBinningHistogram
  get taken from the inhereted HyperBinningPainter class). This allows a HyperPointSet
  to be plotted on top of the HyperVolumeBinning */ 
  HyperPointSet* _hyperPoints;  
  /**< This gets filled in an alternate constructor (usually the HyperBinningHistogram
  get taken from the inhereted HyperBinningPainter class). This allows a HyperPointSet
  to be plotted on top of the HyperVolumeBinning */ 

  const BinningBase& getBinning();

  int getFillColour(double binContents);
  
  void addHyperPoints(TH2D* histogram);

  void drawFilledBins(RootPlotter2D* plotter, bool hashNeg = false);
  void drawFilledBin(RootPlotter2D* plotter, int bin, bool hashNeg = false);
  void drawFilledBin(RootPlotter2D* plotter, HyperCuboid* bin, double binContents);
  void drawFilledBin(RootPlotter2D* plotter, HyperCuboid* bin, int fillColor, int fillStyle);


  void drawBinEdges(RootPlotter2D* plotter);
  void drawBinEdge(RootPlotter2D* plotter, int bin);
  void drawBinEdge(RootPlotter2D* plotter, HyperCuboid* bin);


  void drawBinEdges2(RootPlotter2D* plotter);
  void drawBinEdge2(RootPlotter2D* plotter, int bin, double minWidX, double minWidY);
  void drawBinEdge2(RootPlotter2D* plotter, HyperCuboid* bin, double minWidX, double minWidY);
  void drawBinEdge2(RootPlotter2D* plotter, HyperCuboid* bin, int edge, double minWidX, double minWidY);


  void drawBinNumbers(RootPlotter2D* plotter);
  void drawBinNumbers(RootPlotter2D* plotter, int bin);

  void drawBinCont(RootPlotter2D* plotter);
  void drawBinCont(RootPlotter2D* plotter, int bin);


  public:

  HyperBinningPainter2D(BinningBase* binning, HyperPointSet* hyperPoints = 0);
  HyperBinningPainter2D(HyperHistogram* histogram);

  virtual void draw(TString path = "", TString option = "");

  ~HyperBinningPainter2D();


};



#endif

