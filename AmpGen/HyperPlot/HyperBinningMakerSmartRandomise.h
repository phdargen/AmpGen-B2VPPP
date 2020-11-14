#ifndef HYPERBINNINGMAKERSMARTRANDOMISE_HH
#define HYPERBINNINGMAKERSMARTRANDOMISE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperBinningMaker.h"
#include "AmpGen/HyperPlot/LoadingBar.h"

// Root includes
#include "TMath.h"

// std includes
#include <complex>

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Algorithm to adaptively create a HyperVolumeBinning for a
 *  specific HyperPointSet, taking limits from a HyperCuboid. 
 *
 *  It first performs the same algorithm as HyperBinningMakerSmart, but when
 *  each bin is split, the splitting dimension is chosen at random
 **/
class HyperBinningMakerSmartRandomise : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerSmartRandomise(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerSmartRandomise();
  /**< Destructor  */  
};

#endif
