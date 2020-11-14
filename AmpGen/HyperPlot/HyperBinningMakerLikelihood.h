#ifndef HYPERBINNINGMAKERLIKELIHOOD_HH
#define HYPERBINNINGMAKERLIKELIHOOD_HH

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
 *  lines... Fits the points in the bin with a flat line, then fits it 
 *  with a step fuction (that has 3 free parameters, height before step
 *  height after step, and step point ). The step point that maximises
 *  the likelihood is chosen as the split point. 
 * 
 *  Additionally, you can compare the significance of splitting in different
 *  dimensions using the likelihood ratio of the flat line and the step function.
 *  choose to split in the dimension that gives the biggest significance.
 **/
class HyperBinningMakerLikelihood : public HyperBinningMaker{
  
  private:
  
  public:
  
  HyperBinningMakerLikelihood(const HyperCuboid& binningRange, const HyperPointSet& data);
  /**< Constructor that initiates the base class HyperBinningMaker  */
  virtual void makeBinning();
  /**< run the algorithm  */  
  ~HyperBinningMakerLikelihood();
  /**< Destructor  */  
};

#endif
