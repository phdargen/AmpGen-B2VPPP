/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 **/
 
#ifndef HYPERBINNINGMAKER_HH
#define HYPERBINNINGMAKER_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperCuboid.h"
#include "AmpGen/HyperPlot/HyperPointSet.h"
#include "AmpGen/HyperPlot/HyperBinningMemRes.h"
#include "AmpGen/HyperPlot/RootPlotter1D.h"
#include "AmpGen/HyperPlot/RootPlotter2D.h"
#include "AmpGen/HyperPlot/HyperFunction.h"

// Root includes
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"

// std includes
#include <algorithm>
#include <sstream>

class HyperHistogram;


/**
 * HyperBinningMaker is used to adaptively create a HyperVolumeBinning 
 * from a HyperPointSet. This class just contains a suite of useful tools
 * for splitting bins, and creating the bin hierarchy described in HyperVolumeBinning.
 * Specific types of adaptive binning inheret from this class, and are in HyperBinningMakers.h
 *
 * It is also possible to give it a 'shadow HyperPointSet' - this is useful if one is
 * binning the ratio of two HyperPointSet's and require a minium number of events
 * in both the numerator and denominator.
*/
class HyperBinningMaker {

  private:
  
  HyperBinningMaker();

  protected:


  static bool s_printBinning; /**< print out the status of the binning algorithm */

  //each bin has:
  //   o  a hyper cuboid defining it
  //   o  a set of data which is contained within it
  //   o  a status

  std::vector<HyperCuboid>        _hyperCuboids;          
  /**< vector of HyperCuboid's that defines the binning hierarchy */
  std::vector< std::vector<int> > _linkedBins;           
  /**< each HyperCuboid in the binning hierarchy has linked bins. 
  If these are empty (no linked bins) then this is a true bin. If it 
  contains links, this HyperVolume is just part of the binning hierarchy
  which is used to speed up the binning of events */

  std::vector<HyperPointSet>      _hyperPointSets;        
  /**< records how many HyperPoints in the inital HyperPointSet (that we're going to adaptively bin) 
  fall into each HyperVolume */

  std::vector<HyperPointSet>      _shadowHyperPointSets;  
  /**< records how many HyperPoints in the inital shadow HyperPointSet (that we're going to adaptively bin) 
  fall into each HyperVolume */
  
  /**  Is the volume split as many times as possible, or can I continue to split it  */
  enum VolumeStatus { DONE, CONTINUE };

  std::vector<int>                _status;                
  /**< the status of each HyperVolume i.e. can we continue splitting it into more bins VolumeStatus::CONTINUE 
  or has it been split as many times as possible VolumeStatus::DONE */

  std::vector< std::vector<int> >  _dimSpecificStatus;                
  /**< the status of each HyperVolume dimension i.e. can we continue splitting it into more bins VolumeStatus::CONTINUE 
  or has it been split as many times as possible VolumeStatus::DONE. This only applies to a specific dinemsion.
  i.e. if the bin is < twice the minimum egde length, can safely say that it cannot be resplit in that dimension.
  Note that if all dimensions are VolumeStatus::DONE, then _status should also be VolumeStatus::DONE!! */

  std::vector<int>                _binningDimensions;     
  /**< what dimensions are we allowed to bin in */
  
  bool _shadowAdded;
  /**< has a shadow HyperPointSet been provided  */

  bool _useEventWeights;
  /**< should event weights be considered when creating the binning  */

  double _minimumBinContent;
  /**< the minimum number of events allowed in a bin - if _useEventWeights is true this 
  corresponds to the sum of weights in each bin  */

  double _shadowMinimumBinContent;
  /**< the minimum number of shadow events allowed in a bin - if _useEventWeights is true this 
  corresponds to the sum of weights in each bin  */

  HyperPoint _minimumEdgeLength;
  /**< the minimum bin width in each dimension */

  TRandom* _random;
  /**< random number generator (some binning adaptive binning schemes may require this) */

  bool _drawAlgorithm;
  /**< if this is true, the binning will be drawn after every interation of the algorithm */

  TString _drawAlgorithmDir;
  /**< directory used to draw the binning */

  int _iterationNum;
  /**< what iteration are we on */

  HyperName _names;
  /**< used for the axis titles on any of the HyperBinningHistograms I create */
  
  HyperFunction* _func;
  /**< Some binning algorithms are based on a function rather than a dataset*/

  bool _snapToGrid;
  /**< if this is true, any split point needs to be located on a grid */

  HyperPoint _gridMultiplier;
  /**< if the  */
  

  public:
    
  HyperBinningMaker(const HyperCuboid& binningRange, const HyperPointSet& data);
  
  HyperBinningMaker(const HyperBinning& binning    , const HyperPointSet& data);


  /* ----------------------------------------------------------------*/

  virtual void makeBinning() = 0;
  /**< This function needs to be defined in the derrived classes, and
  runs the entire binning algorithm */


  /* ----------------------------------------------------------------*/

  //The following functions should be called before any binning 
  //binning algorithms commence.

  static void setOutputLevel(bool val){s_printBinning = val;}
  /**< set the verbosity of the output - by default this is on */

  void setBinningDimensions(std::vector<int> dims){_binningDimensions = dims;}
  /**< select which dimensions should be binned */

  void addShadowHyperPointSet(const HyperPointSet& data);

  void setSeed(int seed);

  void useSnapToGrid(bool val);
  void setGridMultiplier(HyperPoint& multipliers);
  void setGridMultiplier(double multiplier);


  void useEventWeights(bool val = true){_useEventWeights = val;}
  /**< select if weighted event should be used - by default this is off */

  void setMinimumBinContent(double val);
  void setShadowMinimumBinContent(double val);
  void setMinimumEdgeLength(double val);     
  void setMinimumEdgeLength(HyperPoint val);  
  
  void setHyperFunction(HyperFunction* fnc);  
  
  void drawAfterEachIteration(TString path);

  void updateFromExistingHyperBinning( const HyperBinning& binning );

  void setNames(HyperName names){_names = names;}
  /**< used for the axis titles on any of the HyperBinningHistograms I create */

  /*----------------------------------------------------------------*/


  int& getGlobalVolumeStatus(int volumeNumber){return _status.at(volumeNumber);}
  int& getDimensionSpecificVolumeStatus(int volumeNumber, int dimension){return _dimSpecificStatus.at(volumeNumber).at(dimension);}

  const int& getGlobalVolumeStatus(int volumeNumber) const {return _status.at(volumeNumber);}
  const int& getDimensionSpecificVolumeStatus(int volumeNumber, int dimension) const {return _dimSpecificStatus.at(volumeNumber).at(dimension);}

  int getNumContinueBins(int dimension = -1) const;
  int getNumBins        () const;
  int getNumHyperVolumes() const;


  /*----------------------------------------------------------------*/
  
  //These functions are common to all the HyperBinning makers
  //and used to divide bins and update their status, and check that
  //the resulting bins are OK.

  int  split(int volumeNumber, int dimension, double splitPoint);

  HyperCuboid splitBelowPoint(int dim, double splitPoint, const HyperCuboid& original, bool noSnapToGrid = false) const;
  HyperCuboid splitAbovePoint(int dim, double splitPoint, const HyperCuboid& original, bool noSnapToGrid = false) const;

  double getSumOfWeights(const HyperPointSet& hyperPointSet) const;
  double getWeight(const HyperPoint& hyperPoint) const;              //Will return 1 if _useEventWeights is false. If not, get event weight 0.
  bool isValidBinningDimension(int dimension);
  virtual bool passFunctionCriteria(HyperCuboid& cuboid1, HyperCuboid& cuboid2);

  HyperPointSet filterHyperPointSet(const HyperPointSet& hyperPointSet, const HyperCuboid& hyperCuboid, bool print = false) const;

  void addBin(const HyperCuboid& hyperCuboid, const HyperPointSet& hyperPointSet, const HyperPointSet& shadowHyperPointSet, int status);

  void setDimSpecStatusFromMinBinWidths (int volumeNumber);
  void updateGlobalStatusFromDimSpecific(int volumeNumber);
  
  bool snapToGrid(const HyperCuboid& cuboid, int dimension, double& splitCoord) const;
  
  /* ----------------------------------------------------------------*/
  
  virtual void startedAlgorithm();
  virtual void startedIteration();  

  void drawCurrentState(TString path) const;

  virtual void finishedIteration();  
  virtual void finishedAlgorithm();

  
  /* ----------------------------------------------------------------*/
  
  // These functions can be used to get a 
  // HyperBinningHistogram of the current state 

  HyperBinningMemRes getHyperVolumeBinning() const;
  HyperHistogram* getHyperBinningHistogram() const;
  HyperHistogram* getShadowHyperBinningHistogram() const;
  HyperHistogram* getRatioHyperBinningHistogram() const;

  virtual ~HyperBinningMaker();


  /********************* NON ESSENTIAL ************************/


  double countEventsBelowSplitPoint(int binNumber, int dimension, double splitPoint) const;
  double countEventsInHyperCuboid(const HyperPointSet& hyperPointSet, const HyperCuboid& hyperCuboid) const;
  double countShadowEventsBelowSplitPoint(int binNumber, int dimension, double splitPoint) const;
  
  //split the bin in a chosen dimension to give a chosen fraction of events in the resulting bin
  double findSmartSplitPoint(int binNumber, int dimension, double dataFraction) const;

  int smartSplit   (int binNumber, int dimension, double dataFraction);
  int smartSplitAll(int dimension, double dataFraction); 
  int smartSplitAllRandomise(double dataFraction = 0.5);


  int smartMultiSplit(int binNumber, int dimension, int parts);
  int smartMultiSplit(int binNumber, int dimension);
  int smartMultiSplitAll(int dimension);


  double findSmartSplitPointInt(int binNumber, int dimension, double dataFraction) const;

  int smartSplitInt(int binNumber, int dimension, double dataFraction);
  int smartSplitAllInt(int dimension, double dataFraction);

  int splitAll(int dimension, double splitPoint);
  int splitAllRandomise(double splitPoint = 0.5);




  //Likelihood shit

  int likelihoodSplit(int binNumber);
  int likelihoodSplitAll();

  int smartLikelihoodSplit(int binNumber);
  int smartLikelihoodSplitAll();

  void  getSplitToMinNeg2LLH(double& split, double& sig, int binNumber, int dimension, bool useConstraints = true);

  void  getDimWithLargestSplitSignificance(int& dim, double& split, int binNumber, bool useConstraints = true);

  TH1D* scanSig(int binNumber, int dimension, int nbins, bool useConstraints = true);

  double neg2LLH(int binNumber, int dimension, double splitPoint, bool useConstraints = true);
  double nullNeg2LLH(int binNumber);





};



#endif

