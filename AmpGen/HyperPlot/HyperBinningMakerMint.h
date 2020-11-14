#ifndef HYPERBINNINGMAKERMINT_HH
#define HYPERBINNINGMAKERMINT_HH

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
 *  It first splits all bins in the 'startingDim' so that each bin is half its size. 
 *  It then splits all the resulting bins in the dimension 'startingDim + 1' 
 *  using the same method.
 *  This process iterates until the minimum bin content or minimum bin widths 
 *  have been reached.
 *
 **/
class HyperBinningMakerMint : public HyperBinningMaker{

  private:

  int _startingDim; /**< the dimension to start splitting from  */

  public:

  HyperBinningMakerMint(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMint();
  /**< Destructor  */  
};

#endif

