#ifndef HYPERBINNINGMAKERPHASEBINNING_HH
#define HYPERBINNINGMAKERPHASEBINNING_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperBinningMaker.h"
#include "AmpGen/HyperPlot/LoadingBar.h"
#include "AmpGen/HyperPlot/CyclicPhaseBins.h"

// Root includes
#include "TMath.h"

// std includes
#include <complex>




/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  This algorithm takes a HyperFunction which describes some kind of phase
 *  motion i.e. map from R^n -> [-pi, pi].
 *  In the final HyperVolumeBinning, each HyperVolume gets assigned a bin number depending on the phase
 *  at that HyperVolume center. The range  [ -pi, pi ] is split uniformly into
 *  a given number of bin pairs (so total bins = bin pairs x 2 ). The numbering
 *  of bins is like so...
 *
 *  [ -pi, pi ] -> [ -_numBinPairs, ..., -1, +1, ... , +_numBinPairs ]
 *  
 *  The algorithm tries to split HyperVolumes such that the bin number does
 *  not change within the HyperVolumes limits.
 *
 *
 **/
class HyperBinningMakerPhaseBinning : public HyperBinningMaker{
  
  private:
  
  int   _numBinPairs;
  CyclicPhaseBins _binEdges;

  int    _maximumRandWalks;
  int    _numWalkers;
  double     _walkSizeFrac;
  
  int _numberOfSystematicSplits;
  int _numberOfGradientSplits;
  
  private:

  int splitByCoord(int volumeNumber, int dimension, HyperPoint& coord);
  /**< Split a HyperVolume in a given dimension at a given coord */

  

  //Functions for determining the bin boundaries / bin number from a given bin number / phase.

  int getBinNumFromFunc(HyperPoint& point);
  /**< For a given HyperPoint, get the bin number */    
  int getBinNumFromFuncVal(double phase);
  /**< For a given phase, get the bin number */    
  double getLowBinBoundary(double phase);
  /**< For a given bin, get the lower phase boundary */  
  double getHighBinBoundary(double phase);
  /**< For a given bin, get the upper phase boundary */
  double closestBinBoundary(double val);
  /**< For a given phase, determine the closest bin boundary (as a phase)*/

  //Functions to evaluate derivatives and second derivatives of the phase motion

  HyperPoint getGrad(HyperPoint& point);
  /**< Evaluate the gradient at a given point. Step length is the minimium edge 
  length. Function is evaluated at (x+h and x-h) */

  HyperPoint getGradPos(HyperPoint& point, double funcValAtPoint);
  /**< Evaluate the gradient at a given point. Step length is the minimium edge 
  length. Function is evaluated at (x+h) only */

  double getSecondDerivative(HyperPoint& point, HyperPoint& vector, double funcValAtPoint, double& deriv);
  /**< Evaluate the second and first derivative in a given direction at a given point */  

  int splitDimFromGrad(int volumeNumber, HyperPoint gradient);
  /**< Use the graident at the center of the given volume number
  to decide what dimension to split in. */

  virtual int gradientSplit(int binNumber   , int& dimension);
  /**< The default split method. It uses the gradient of the 
  function to predict where the boundary wetween two bin numbers is. If
  it fails, it will call the systematic split instead */


  int systematicSplit      (int volumeNumber, int dimension, double valAtCenter, HyperPoint gradient);
  /**< This function looks at the function value at specific points within the HyperVolume. It uses
  the corners of the HyperVolume, a point in the middle of each face (NPlane), and a point in the middle
  of each edge. If the bin number changes at any of these points, it will result in the bin being split.
  To speed things up it uses the gradient at the HyperVolume center to order the points, on which are most
  liekly to result in a successful split. While looping through the points it estimates the uncertainty
  of the function estimate (from the graient). If the estimate indicates that the function is more than
  5 sigma from a bin split, it will skip that point (and all remaining points since they are ordered) */

  HyperPointSet getSplitCorners( int volumeNumber );
  /**< Get a list of all the corners of the volume number given */
  HyperPointSet getSplitFaces  ( int volumeNumber );
  /**< Get a list of all the face centers of the volume number given */
  HyperPointSet getSplitEdges  ( int volumeNumber );
  /**< Get a list of all the edge centers of the volume number given */

  HyperPoint orderAndTestSplitPoints(HyperPointSet& points, HyperPoint& point, double valAtPoint, HyperPoint gradient);
  /**< See systematicSplit for a full description. */


  int randomWalkSplit(int volumeNumber, int dimension);
  /**< Do random walk within the HyperVolume to search for a point where the
  bin number changes. Not really used anymore. */

  void walkOrthogonal(HyperPoint& point, HyperCuboid& walkLimits);
  /**< Do a random walk in one dimension (within the walk limits) */

  void walk(HyperPoint& point, HyperCuboid& walkLimits);
  /**< Do a random walk in a random dimension (within the walk limits) */



  public:
  
  HyperBinningMakerPhaseBinning(const HyperCuboid& binningRange, HyperFunction* func);
  /**< Constructor that initiates the base class HyperBinningMaker  */

  virtual void makeBinning();
  /**< run the algorithm  */  

  void setNumBinPairs(int binpairs); 
  /**< set the number of bin pairs (see class description of details) */  
  
  void setBinEdges(std::vector<double> binEdges);
  /**< set the bin edges - the highest bin edge is the lowest + pi, so only nBinPair edges need to be passed */  

  virtual int gradientSplitAll();


  ~HyperBinningMakerPhaseBinning();
  /**< Destructor  */  
};

#endif

