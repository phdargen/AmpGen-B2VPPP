#ifndef HYPERBINNINGMAKERMINTSMART_HH
#define HYPERBINNINGMAKERMINTSMART_HH

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
 *  It first performs the same algorithm as HyperBinningMakerMint, and then
 *  follows this with HyperBinningMakerSmart
 *
 **/
class HyperBinningMakerMintSmart : public HyperBinningMaker{

  private:

  int _startingDim; /**< the dimension to start splitting from  */

  public:

  HyperBinningMakerMintSmart(const HyperCuboid& binningRange,const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerMintSmart();
  /**< Destructor  */  
};

#endif
