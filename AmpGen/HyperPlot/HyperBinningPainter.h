/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class plotting HyperBinningHistograms of any dimension
 *
 **/
 
#ifndef HYPERBINNINGPAINTER_HH
#define HYPERBINNINGPAINTER_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/RootPlotter1D.h"
#include "AmpGen/HyperPlot/RootPlotter2D.h"

#include "AmpGen/HyperPlot/BinningBase.h"
#include "AmpGen/HyperPlot/HyperHistogram.h"


// Root includes
#include "TH1D.h"
#include "TH2D.h"

// std includes



class HyperBinningPainter {
  
  protected:

  HyperHistogram* _histogram; /**< Pointer to the HyperBinningHistogram I'm going to plot */
  bool _density; /**< Am I drawing the bin contents or the density? */ 

  public:
  
  HyperBinningPainter(HyperHistogram* histogram);
  
  void useDensity(bool val){_density = val;}
  /**< Am I drawing the bin contents or the density? */ 

  virtual void draw(TString path = "", TString option = "");

  virtual ~HyperBinningPainter();

};


#endif

