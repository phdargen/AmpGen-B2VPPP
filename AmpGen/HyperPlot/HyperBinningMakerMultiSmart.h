#ifndef HYPERBINNINGMAKERMULTISMART_HH
#define HYPERBINNINGMAKERMULTISMART_HH

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
 *  \todo describe how this works... 
 *
 **/
class HyperBinningMakerMultiSmart : public HyperBinningMaker {
  
  private:
  
  int _startingDim;  /**< the dimension to start splitting from  */
  
  public:
  
  HyperBinningMakerMultiSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMultiSmart();
  /**< Destructor  */  
  
};

#endif

