#ifndef HYPERBINNINGMAKERSMART_HH
#define HYPERBINNINGMAKERSMART_HH

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
 *  It first splits all bins in the 'startingDim' so that each bin contains 50%
 *  of the HyperPoints. 
 *  It then splits all the resulting bins in the dimension 'startingDim + 1' 
 *  using the same method.
 *  This process iterates until the minimum bin content or minimum bin widths 
 *  have been reached.
 *
 **/
class HyperBinningMakerSmart : public HyperBinningMaker {
  
  private:
  
  int _startingDim;  /**< the dimension to start splitting from  */
  
  public:
  
  HyperBinningMakerSmart(const HyperCuboid& binningRange, const HyperPointSet& data, int startingDim = 0);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerSmart();
  /**< Destructor  */  
  
};


#endif

