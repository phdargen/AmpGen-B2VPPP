#ifndef HYPERBINNINGMAKERSMARTLIKELIHOOD_HH
#define HYPERBINNINGMAKERSMARTLIKELIHOOD_HH

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
 *  \todo rememer how this algorithm works - its something along these
 *  lines... follows the same proceedure as HyperBinningMakerLikelihood
 *  but only the dimension is picked using the likelihood method. The
 *  splitting is done using the SMART method (HyperBinningMakerSmart) .
 *
 **/
class HyperBinningMakerSmartLikelihood : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerSmartLikelihood(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerSmartLikelihood();
  /**< Destructor  */  
};

#endif


