#include "AmpGen/HyperPlot/HyperBinningMaker.h"

#include "AmpGen/HyperPlot/HyperHistogram.h"

///Just an empty contructor to make it compile (private, will never be used)
///
HyperBinningMaker::HyperBinningMaker() :
  _minimumEdgeLength(0),
  _drawAlgorithm          (false),  
  _iterationNum           (0    ),
  _names( 0 ),
  _func( 0 ),
  _snapToGrid(false),  
  _gridMultiplier( 0 )
{

}

///Constuctor which enables the user to create a binning scheme 
///for the HyperPointSet given, and using limits provided by a HyperCuboid.
///
HyperBinningMaker::HyperBinningMaker(const HyperCuboid& binningRange, const HyperPointSet& data) :
  _hyperCuboids           (1, binningRange),
  _linkedBins             (1, std::vector<int>(0,0) ),
  _hyperPointSets         (1, filterHyperPointSet( data, binningRange, true )       ),
  _shadowHyperPointSets   (1, HyperPointSet( binningRange.getDimension() ) ), 
  _status                 (1, VolumeStatus::CONTINUE ),
  _dimSpecificStatus      (1, std::vector<int>( binningRange.getDimension(), VolumeStatus::CONTINUE )  ),  
  _shadowAdded            (false),
  _useEventWeights        (false),
  _minimumBinContent      (15),
  _shadowMinimumBinContent(15),
  _minimumEdgeLength      ( HyperPoint(binningRange.getDimension(), 0.1) ),
  _random                 (new TRandom3(0)),
  _drawAlgorithm          (false),
  _iterationNum           (0    ),
  _names( binningRange.getDimension() ),
  _func (0),
  _snapToGrid(false),
  _gridMultiplier( HyperPoint(binningRange.getDimension(), 1) )

{  
  for (int i = 0; i < binningRange.getDimension(); i++) _binningDimensions.push_back(i);
  WELCOME_LOG << "Good day from the HyperBinningMaker() Constructor"<<std::endl;

}
 

///Takes an existing HyperBinning and updates the current state of the HyperBinningMaker
///
void HyperBinningMaker::updateFromExistingHyperBinning( const HyperBinning& binning ){
  
  if (_hyperPointSets.size() != 1) {
    ERROR_LOG << "You can only call the updateFromExistingHyperBinning if the HyperBinningMaker has" << std::endl;
    ERROR_LOG << "a single HyperVolume";
  }
  else{
    INFO_LOG << "Congratulations, you've taken the bold step of calling the updateFromExistingHyperBinning function" << std::endl;
  }
  
  HyperPointSet points     = _hyperPointSets      .at(0);
  HyperPointSet shadPoints = _shadowHyperPointSets.at(0); 

  _hyperCuboids         .clear();   
  _linkedBins           .clear();
  _status               .clear();  
  _dimSpecificStatus    .clear();  
  

  int nvols = binning.getNumHyperVolumes();
  
  for (int i = 0; i < nvols; i++){
    _hyperCuboids.push_back( binning.getHyperVolume(i).at(0)  );
    _linkedBins  .push_back( binning.getLinkedHyperVolumes(i) );

    if ( _linkedBins.at(i).size() == 0 ){
      _status           .push_back( VolumeStatus::CONTINUE );
      _dimSpecificStatus.push_back( std::vector<int>( binning.getDimension(), VolumeStatus::CONTINUE ) );
    }
    else {
      _status           .push_back( VolumeStatus::DONE );
      _dimSpecificStatus.push_back( std::vector<int>( binning.getDimension(), VolumeStatus::DONE ) );
    }

  }
  
  _hyperPointSets.clear();  
  _hyperPointSets.resize(nvols, HyperPointSet(binning.getDimension()) );

  std::vector<int> binNums = binning.getBinNum(points);

  for (unsigned i = 0; i < points.size(); i++){
    int binNum = binNums.at(i);
    if (binNum == -1) continue;
    int volNum = binning.getHyperVolumeNumber( binNum );
    _hyperPointSets.at(volNum).push_back( points.at(i) );
  }

  _shadowHyperPointSets.clear();  
  _shadowHyperPointSets.resize(nvols, HyperPointSet(binning.getDimension()) );


  std::vector<int> shadBinNums = binning.getBinNum(shadPoints);
  
  for (unsigned i = 0; i < shadPoints.size(); i++){
    int binNum = shadBinNums.at(i);
    if (binNum == -1) continue;
    int volNum = binning.getHyperVolumeNumber( binNum );
    _hyperPointSets.at(volNum).push_back( shadPoints.at(i) );
  }  
  
  
}




///Add a shadow HyperPointSet to the HyperBinningMaker. This must be done
///before any adaptive binning algorithms commence.
///
void HyperBinningMaker::addShadowHyperPointSet(const HyperPointSet& data){
  if (_hyperCuboids.size() != 1) ERROR_LOG << "You have to add the shadow set before you make the binning"<<std::endl;
  else{
    _shadowHyperPointSets.at(0) = filterHyperPointSet(data, _hyperCuboids.at(0), true);
    _shadowAdded = true;
  }
}

///Split a HyperCuboid in a specific dimension, at point x_dim where: 
///
///x_dim = lowEdge_dim + (highEdge_dim - lowEdge_dim)*splitPoint
///
///i.e. splitPoint = 0.5 would split the bin into two equal pieces.
///
///This fuction returns the HyperCuboid above the split point 
HyperCuboid HyperBinningMaker::splitAbovePoint(int dim, double splitPoint, const HyperCuboid& original, bool noSnapToGrid) const{

  HyperPoint lowCorner   ( original.getLowCorner () );
  HyperPoint highCorner  ( original.getHighCorner() );

  double splitCoord = lowCorner.at(dim) + (highCorner.at(dim) - lowCorner.at(dim))*splitPoint;
  if ( _snapToGrid == true && noSnapToGrid == false ) {
    if ( snapToGrid(original, dim, splitCoord) == false ){
      return HyperCuboid(0);
    }  
  }

  HyperPoint newLowCorner( lowCorner );
  newLowCorner.at(dim) = splitCoord;
  
  HyperCuboid temp(newLowCorner, highCorner);
  return temp;

} 

///Split a HyperCuboid in a specific dimension, at point x_dim where: 
///
///x_dim = lowEdge_dim + (highEdge_dim - lowEdge_dim)*splitPoint
///
///i.e. splitPoint = 0.5 would split the bin into two equal pieces.
///
///This fuction returns the HyperCuboid below the split point 
HyperCuboid HyperBinningMaker::splitBelowPoint(int dim, double splitPoint, const HyperCuboid& original, bool noSnapToGrid) const{

  HyperPoint lowCorner   ( original.getLowCorner () );
  HyperPoint highCorner  ( original.getHighCorner() );

  double splitCoord = lowCorner.at(dim) + (highCorner.at(dim) - lowCorner.at(dim))*splitPoint;
  if ( _snapToGrid == true && noSnapToGrid == false ) {
    if ( snapToGrid(original, dim, splitCoord) == false ){
      return HyperCuboid(0);
    }
  }

  HyperPoint newHighCorner( highCorner );
  newHighCorner.at(dim) = splitCoord;

  HyperCuboid temp(lowCorner, newHighCorner);
  return temp;

}

///Return any HyperPoint's that are within the HyperCuboid
///
HyperPointSet HyperBinningMaker::filterHyperPointSet(const HyperPointSet& hyperPointSet, const HyperCuboid& hyperCuboid, bool print) const{
  
  if (print == true && s_printBinning == true) INFO_LOG << "Before filtering there are " << hyperPointSet.size() << " events" << std::endl;

  HyperPointSet temp(hyperPointSet.getDimension());
  for(unsigned int i = 0; i < hyperPointSet.size(); i++){
    if( hyperCuboid.inVolume(hyperPointSet.at(i)) ) temp.push_back(hyperPointSet.at(i));
  }

  if (print == true && s_printBinning == true) INFO_LOG << "After filtering there are " << temp.size() << " events" << std::endl;

  return temp;

}

///Set the random seed that is required for some binning schemes
///
void HyperBinningMaker::setSeed(int seed){

  delete _random;
  _random = new TRandom3(seed);

}

///Do you want to use the snap to grid feature i.e. force all the bin egdes 
///to lie on a well defined grid. Grid is has N bins in each direction between
///min and max (from the binning limits given)
///where N = floor( (max - min) / minBinWidth) * gridMultiplier
void HyperBinningMaker::useSnapToGrid(bool val){
  _snapToGrid = val;
}

///Set the grid mulipliers independently for each dimension. See 'useSnapToGrid'
///function for description
void HyperBinningMaker::setGridMultiplier(HyperPoint& multipliers){
  _gridMultiplier = multipliers;
}

///Set the grid muliplier to be the same in each dimesnion. See 'useSnapToGrid'
///function for description
void HyperBinningMaker::setGridMultiplier(double multiplier){
  HyperPoint point( _hyperCuboids.at(0).getDimension(), multiplier );
  setGridMultiplier( point );
}

///Get the sum of weights in a HyperPointSet. If _useEventWeights isn't true,
///just returns the number of events in the HyperPointSet.
double HyperBinningMaker::getSumOfWeights(const HyperPointSet& hyperPointSet) const{

  if (_useEventWeights == true) return hyperPointSet.getSumW();
  return hyperPointSet.size();
  
}


///Get the weight of an event - if _useEventWeights isn't true,
///this is just 1.0
double HyperBinningMaker::getWeight(const HyperPoint& hyperPoint) const{

  if (_useEventWeights == true) return hyperPoint.getWeight();
  return 1.0;

}


///Add a bin to the binning scheme 
///
void HyperBinningMaker::addBin(const HyperCuboid& hyperCuboid, const HyperPointSet& hyperPointSet, const HyperPointSet& shadowHyperPointSet, int status){
  _hyperPointSets.push_back(hyperPointSet);
  _shadowHyperPointSets.push_back(shadowHyperPointSet);
  _hyperCuboids  .push_back(hyperCuboid);
  _status        .push_back(status);
  _linkedBins    .push_back( std::vector<int>(0,0) );

  std::vector<int> dimStatus( hyperCuboid.getDimension(), status);
  _dimSpecificStatus.push_back(dimStatus);

}

///Check that we are allowed to bin in this dimension
///
bool HyperBinningMaker::isValidBinningDimension(int dimension){

  bool isBinningDim = false;
  for ( unsigned i = 0; i < _binningDimensions.size(); i++ ) { 
    if (_binningDimensions.at(i) == dimension) {
      isBinningDim = true;
      break;
    }
  }
  if (isBinningDim == false) return false;
  return true;

}

/**
  This checks the bin width in each dimension agaisnt the minimum bin width for
  that dimension. This is used to set the dimension speicific bin status. 
  If every status is DONE then the global status can also be set to DONE.

  Dimensions that are not in the list of binning dimensions are also, automatically
  set to DONE.  
*/
void HyperBinningMaker::setDimSpecStatusFromMinBinWidths(int volumeNumber){

  HyperCuboid&   chosenHyperCuboid          = _hyperCuboids        .at(volumeNumber);
  
  int dim = _hyperCuboids.at(0).getDimension();
  
  bool allDone = true;

  for (int i = 0; i < dim; i++){
    double low  =  chosenHyperCuboid.getLowCorner ().at(i);
    double high =  chosenHyperCuboid.getHighCorner().at(i);
    double width    = high - low;
    double minwidth = _minimumEdgeLength.at(i);

    if ( width < minwidth*2.0 || isValidBinningDimension(i) == false ) {
      getDimensionSpecificVolumeStatus(volumeNumber, i) = VolumeStatus::DONE;
    }
    else{
      allDone = false;
    }

  }  

  if (allDone == true){
    getGlobalVolumeStatus(volumeNumber) = VolumeStatus::DONE;
  }

}

/**
  For a given volume, see if all the dimension specific statuses are DONE.
  If so, set the global status to DONE.
*/
void HyperBinningMaker::updateGlobalStatusFromDimSpecific(int volumeNumber){

  int dim = _hyperCuboids.at(0).getDimension();
  
  bool allDone = true;

  for (int i = 0; i < dim; i++){

    int status = getDimensionSpecificVolumeStatus(volumeNumber, i);
    
    if (status != VolumeStatus::DONE){
      allDone = false;
      break;
    }

  }  

  if (allDone == true){
    getGlobalVolumeStatus(volumeNumber) = VolumeStatus::DONE;
  }

}

/**
This is the most important function in the class. It is used to split a specific
hypervolume in a specific dimension, at point x_dim where: 

x_dim = lowEdge_dim + (highEdge_dim - lowEdge_dim)*splitPoint

i.e. splitPoint = 0.5 would split the bin into two equal peices.

While splitting the bin it checks if the resulting bins follow the nessesary 
requirements. These are:

  - Number of events (or sum of weights) in each of the resulting bins 
    is greater than the minimum bin content (_minimumBinContent)
  - Number of shadow events (or sum of weights) in each of the resulting bins 
    is greater than the minimum bin content (_shadowMinimumBinContent)
  - The width of the resulting bins are larger than the minimum bin width
    ( _minimumEdgeLength.at(dimension) )
  - The dimension that is being split is a valid binning direction
    (i.e an element of _binningDimensions)

*/
int HyperBinningMaker::split(int volumeNumber, int dimension, double splitPoint){
  
  //check if we're allowed to bin in this dimension

  VERBOSE_LOG << "Calling the HyperBinningMaker::split function that gets used by all the algorithms" << std::endl;    


  if (isValidBinningDimension(dimension) == false) {
    VERBOSE_LOG << "This isn't a valid binning dimension, so not splitting" << std::endl;    
    return 0;
  }
  
  //get the HyperPointSet, ShadowHyperPointSet, and HyperCuboid associated to
  //this volume number
  HyperPointSet& chosenHyperPointSet        = _hyperPointSets      .at(volumeNumber);
  HyperPointSet& chosenShadowHyperPointSet  = _shadowHyperPointSets.at(volumeNumber);
  HyperCuboid&   chosenHyperCuboid          = _hyperCuboids        .at(volumeNumber);
  
  //create the two HyperCuboid's that are been split
  HyperCuboid cuboid1 = splitBelowPoint(dimension, splitPoint, chosenHyperCuboid);
  HyperCuboid cuboid2 = splitAbovePoint(dimension, splitPoint, chosenHyperCuboid);
  
  if (cuboid1.getDimension() == 0 || cuboid2.getDimension() == 0){
    VERBOSE_LOG << "It looks like the snap to grid option means that this bin cannot be split. Returning 0."<<std::endl;
    return 0;    
  }

  //calcuate the edge length of each HyperCuboid in the binning dimension
  double edgeLength1   = cuboid1.getHighCorner().at(dimension) - cuboid1.getLowCorner().at(dimension);
  double edgeLength2   = cuboid2.getHighCorner().at(dimension) - cuboid2.getLowCorner().at(dimension);
  double minEdgeLength = _minimumEdgeLength.at(dimension);

  //first check if the new bins are too small - if they are, return 0
  if (edgeLength1 < minEdgeLength || edgeLength2 < minEdgeLength){
    VERBOSE_LOG << "Tired to split bin but one of the resulting bins is too small... hopefully spliting in another dim will help."<<std::endl;
    return 0;
  }  

  //find the HyperPoint's that fall into each of the new HyperCuboid's
  HyperPointSet hyperPointSet1 = filterHyperPointSet(chosenHyperPointSet, cuboid1);
  HyperPointSet hyperPointSet2 = filterHyperPointSet(chosenHyperPointSet, cuboid2);
  
  //find the shadow HyperPoint's that fall into each of the new HyperCuboid's
  HyperPointSet shadowHyperPointSet1 = filterHyperPointSet(chosenShadowHyperPointSet, cuboid1);
  HyperPointSet shadowHyperPointSet2 = filterHyperPointSet(chosenShadowHyperPointSet, cuboid2);
  
  double evts1       = getSumOfWeights(hyperPointSet1      );
  double evts2       = getSumOfWeights(hyperPointSet2      );
  double shadowEvts1 = getSumOfWeights(shadowHyperPointSet1);
  double shadowEvts2 = getSumOfWeights(shadowHyperPointSet2);

  //check if there are enough HyperPoint's in the split bins - if not, return 0
  if ( evts1 < _minimumBinContent || evts2 < _minimumBinContent){
    VERBOSE_LOG << "Tired to split bin but one half has too little events... hopefully spliting in another dim will help."<<std::endl;
    VERBOSE_LOG << "It contained " << chosenHyperPointSet.size() << " events, and was split into " << hyperPointSet1.size() << " and " << hyperPointSet2.size()<<std::endl;
    return 0;
  }

  //check if there are enough shadow HyperPoint's in the split bins - if not, return 0
  if (_shadowAdded == true){
    if ( shadowEvts1 < _shadowMinimumBinContent || shadowEvts2 < _shadowMinimumBinContent){
      VERBOSE_LOG << "Tired to split bin but one half has too little events... hopefully spliting in another dim will help."<<std::endl;
      VERBOSE_LOG << "It contained " << chosenHyperPointSet.size() << " events, and was split into " << hyperPointSet1.size() << " and " << hyperPointSet2.size()<<std::endl;
      return 0;
    }
  }

  //There is an option to pass a function - if so check that the new
  //bins pass the 'function criteria' (in passFunctionCriteria).
  //Default behaviour returns true, but some Algs may override this
  if (_func != 0){
    if ( passFunctionCriteria(cuboid1, cuboid2) == 0 ) {
      VERBOSE_LOG << "I tired to split this bin but the Function criteria failed." <<std::endl;
      return 0;
    }
  }

  // Add bins to the HyperVolume vector 
  // if the bin content is less than double the _minimumBinContent,
  // there is no way to split it any further. In this case, mark
  // the bin as DONE. In not mark the bin as CONTINUE

  if ( evts1 >= 2.0*_minimumBinContent && (shadowEvts1 >= 2.0*_shadowMinimumBinContent || _shadowAdded == false) ){
    VERBOSE_LOG << "Adding HyperVolume 1 with status 1";  
    addBin(cuboid1, hyperPointSet1, shadowHyperPointSet1, VolumeStatus::CONTINUE);
  } 
  else{
    VERBOSE_LOG << "Adding HyperVolume 1 with status 0"<<std::endl;  
    addBin(cuboid1, hyperPointSet1, shadowHyperPointSet1, VolumeStatus::DONE    );
  }
  
  if ( evts2 >= 2.0*_minimumBinContent && (shadowEvts2 >= 2.0*_shadowMinimumBinContent || _shadowAdded == false) ){
    VERBOSE_LOG << "Adding HyperVolume 2 with status 1"<<std::endl; 
    addBin(cuboid2, hyperPointSet2, shadowHyperPointSet2, VolumeStatus::CONTINUE);
  } 
  else{
    VERBOSE_LOG << "Adding HyperVolume 2 with status 0"<<std::endl; 
    addBin(cuboid2, hyperPointSet2, shadowHyperPointSet2, VolumeStatus::DONE    );
  }
  
  //Link the old bin to the new bins

  int newVolumeNum1 = _hyperPointSets.size() - 2;
  int newVolumeNum2 = _hyperPointSets.size() - 1; 
  
  VERBOSE_LOG << "Linking old bin to new bins"<<std::endl;
  _linkedBins.at( volumeNumber ).push_back( newVolumeNum1 );
  _linkedBins.at( volumeNumber ).push_back( newVolumeNum2 );
  
  //clear the HyperPoints and Shadow HyperPoints associated
  //to the original volume.

  int dim = _hyperPointSets.at(volumeNumber).getDimension();

  VERBOSE_LOG << "Removing data associated with old bin..." << std::endl;
  _hyperPointSets      .at(volumeNumber) = HyperPointSet(dim);
  _shadowHyperPointSets.at(volumeNumber) = HyperPointSet(dim);
  VERBOSE_LOG << "and setting it's status to 0" << std::endl;

  //finally, set the status of the original volume to DONE

  _status.at(volumeNumber) = VolumeStatus::DONE;

  setDimSpecStatusFromMinBinWidths(newVolumeNum1);
  setDimSpecStatusFromMinBinWidths(newVolumeNum2);

  return 1;


  //if (_verbose){
  //  std::stringstream ss;
  //  chosenHyperCuboid.print(ss, 0);
  //  VERBOSE_LOG << "Original Cuboid:  " << ss.str();
  //  std::stringstream ss2;
  //  cuboid1.print(ss2, 0);
  //  VERBOSE_LOG << "Split Cuboid 1:  " << ss2.str();
  //  std::stringstream ss3;
  //  cuboid2.print(ss3, 0);
  //  VERBOSE_LOG << "Split Cuboid 2:  " << ss3.str();  
  //  VERBOSE_LOG << "HyperCuboid 1 contains " << getSumOfWeights( hyperPointSet1 ) << "events";
  //  VERBOSE_LOG << "HyperCuboid 2 contains " << getSumOfWeights( hyperPointSet2 ) << "events";  
  //}

}

/// Sometimes it's nice to have all the bin edges along some kind of grid structure.
/// My hope is that this also allows better compression of the binning scheme.
bool HyperBinningMaker::snapToGrid(const HyperCuboid& cuboid, int dimension, double& splitCoord) const{
  
  double absMax = _hyperCuboids.at(0).getHighCorner ().at(dimension);
  double absMin = _hyperCuboids.at(0).getLowCorner  ().at(dimension);
  
  double minEdgeLength = _minimumEdgeLength.at(dimension);
  
  int gridMultiplier = floor( _gridMultiplier.at(dimension) );
  int nbins = floor((absMax - absMin)/minEdgeLength)*gridMultiplier;

  double gridWidth = (absMax - absMin)/double(nbins);
  
  
  int binNumLow  = floor( (splitCoord - absMin)/gridWidth );
  int binNumHigh = ceil ( (splitCoord - absMin)/gridWidth );

  double lowCoord  = absMin + double(binNumLow )*gridWidth;
  double highCoord = absMin + double(binNumHigh)*gridWidth;
  

  double closestGridPoint  = 0.0;
  double furthestGridPoint = 0.0;


  if ( fabs(lowCoord - splitCoord) <= fabs(highCoord - splitCoord) ){
    closestGridPoint  = lowCoord ;
    furthestGridPoint = highCoord;
  }
  else{
    closestGridPoint  = highCoord;
    furthestGridPoint = lowCoord ;    
  }

  double min = cuboid.getLowCorner ().at(dimension);
  double max = cuboid.getHighCorner().at(dimension);
  
  //std::cout << min << "    " << lowCoord << "   " << splitCoord << "    " << highCoord << "    " << max << std::endl;


  if ( closestGridPoint - min >= minEdgeLength && max - closestGridPoint >= minEdgeLength ){
    splitCoord = closestGridPoint;
    return true;
  }
  if ( furthestGridPoint - min >= minEdgeLength && max - furthestGridPoint >= minEdgeLength ){
    splitCoord = furthestGridPoint;
    return true;    
  }   
  
  VERBOSE_LOG << "I failed to snap to grid - very sad" << std::endl;

  splitCoord = closestGridPoint;
  return false;

}




/// If a HyperFunction is passed through the setHyperFunction function,
/// this function is called, and must return true for any successful bin split.
/// Default behavour is to return true, but other derived classes can override
bool HyperBinningMaker::passFunctionCriteria(HyperCuboid& cuboid1, HyperCuboid& cuboid2){
  
  cuboid1.getLowCorner();
  cuboid2.getHighCorner();

  return true; 
}



/// \todo remember how this works
///
///
double HyperBinningMaker::neg2LLH(int binNumber, int dimension, double splitPoint, bool useConstraints){
  
  double binLength = _hyperCuboids.at(binNumber).getHighCorner().at(dimension) - _hyperCuboids.at(binNumber).getLowCorner().at(dimension);

  double total  = getSumOfWeights( _hyperPointSets.at(binNumber) );
  double nBelow = countEventsBelowSplitPoint(binNumber, dimension, splitPoint);
  double nAbove = total - nBelow;
  
  if (useConstraints == true){
    if (nBelow < _minimumBinContent || nAbove < _minimumBinContent) return nullNeg2LLH(binNumber);
    if (binLength*splitPoint < _minimumEdgeLength.at(dimension) || binLength*(1.0 - splitPoint) < _minimumEdgeLength.at(dimension)) return nullNeg2LLH(binNumber);
  }

  if (_shadowAdded == false){

    double fracBelow = nBelow/total;
    double fracAbove = nAbove/total;
    
    double PDFabove = fracAbove/(1.0 - splitPoint);
    double PDFbelow = fracBelow/splitPoint;
    
    if (PDFabove  == 0.0) PDFabove = 1e-16;
    if (PDFbelow  == 0.0) PDFbelow = 1e-16;

    return -2.0*(nAbove*log(PDFabove) + nBelow*log(PDFbelow));
  }
  
  double shadowTotal = getSumOfWeights( _shadowHyperPointSets.at(binNumber) );
  double shadowNBelow = countShadowEventsBelowSplitPoint(binNumber, dimension, splitPoint);
  double shadowNAbove = shadowTotal - shadowNBelow;
  
  double allTotal = shadowTotal + total;

  double PDFNumeratorAbove   = (nAbove       / allTotal) / (1.0 -splitPoint);
  double PDFDenominatorAbove = (shadowNAbove / allTotal) / (1.0 -splitPoint);
  double PDFNumeratorBelow   = (nBelow       / allTotal) / (splitPoint);
  double PDFDenominatorBelow = (shadowNBelow / allTotal) / (splitPoint);
  
  if (PDFNumeratorAbove    <= 0.0) PDFNumeratorAbove   = 1e-16;
  if (PDFDenominatorAbove  <= 0.0) PDFDenominatorAbove = 1e-16;
  if (PDFNumeratorBelow    <= 0.0) PDFNumeratorBelow   = 1e-16;
  if (PDFDenominatorBelow  <= 0.0) PDFDenominatorBelow = 1e-16;

  return -2.0*( shadowNBelow*log(PDFDenominatorBelow) + shadowNAbove*log(PDFDenominatorAbove) + nBelow*log(PDFNumeratorBelow) + nAbove*log(PDFNumeratorAbove) );

}



/// \todo remember how this works
///
///
double HyperBinningMaker::nullNeg2LLH(int binNumber){

  if (_shadowAdded == false) return 0.0;
  
  double num  = getSumOfWeights( _hyperPointSets.at(binNumber) );
  double sha  = getSumOfWeights( _shadowHyperPointSets.at(binNumber) );
  double tot = num + sha;

  double probNum = num/tot;
  double probDen = sha/tot;

  return -2.0*(num*log(probNum) + sha*log(probDen));

}

/// \todo remember how this works
///
///
TH1D* HyperBinningMaker::scanSig(int binNumber, int dimension, int nbins, bool useConstraints){

  TH1D* scan = new TH1D("scan", "scan", nbins, 0.0, 1.0);
  double nullHypothesis = nullNeg2LLH(binNumber);

  for (int i = 1; i < nbins; i++){
    double neg2LLHval = neg2LLH(binNumber, dimension, scan->GetXaxis()->GetBinCenter(i), useConstraints );
    neg2LLHval = nullHypothesis - neg2LLHval;
    if (neg2LLHval < 0.0) neg2LLHval = 0.0;
    double sig = sqrt( neg2LLHval );
    scan->SetBinContent( i, sig );
  }
  
  return scan;

}

/// \todo remember how this works
///
///
void HyperBinningMaker::getSplitToMinNeg2LLH(double& split, double& sig, int binNumber, int dimension, bool useConstraints){

  bool plot = false;

  int nBins = 25;
  
  //std::cout << "OK here?" << std::endl;

  TH1D* splitHist = scanSig(binNumber, dimension, nBins, useConstraints);

  int splitBin = splitHist->GetMaximumBin();

  split    = splitHist->GetXaxis()->GetBinCenter(splitBin);
  sig      = splitHist->GetBinContent(splitBin);
  
  //std::cout << "OK here?" << sig << std::endl;

  if (plot == true){
    RootPlotter1D plotter(splitHist);
    TString name = "splitHist";
    name += binNumber;
    name += "_";
    name += dimension;
    plotter.plot(name);
    INFO_LOG << "Min split point = " << split << " in dimension " << dimension<<std::endl;
  }

  delete splitHist;

}

/// \todo remember how this works
///
///
void HyperBinningMaker::getDimWithLargestSplitSignificance(int& dim, double& split, int binNumber, bool useConstraints){

  int nDim = _binningDimensions.size();

  double bigestSig = -2e40;

  for (int i = 0; i < nDim; i++){

    double tmpSplit = 0.0; double tmpSig = 0.0;
    getSplitToMinNeg2LLH(tmpSplit, tmpSig, binNumber, _binningDimensions.at(i), useConstraints);

    if ( tmpSig > bigestSig ) {
      bigestSig = tmpSig;
      dim       = _binningDimensions.at(i);
      split     = tmpSplit;
    }

  }

}

/// \todo remember how this works
///
///
int HyperBinningMaker::likelihoodSplit(int binNumber){
  
  int      dim = -1;
  double splitPoint = -1.0;

  getDimWithLargestSplitSignificance(dim, splitPoint, binNumber, true);
  
  //INFO_LOG << "likelihoodSplit(" << binNumber << ") gives dim " << dim << " and split point " <<  splitPoint;

  splitPoint = _random->Gaus(splitPoint, 0.2);  //I think this is needed to make the binning unbiased
  if (splitPoint <= 0.0 || splitPoint >= 1.0) return likelihoodSplit(binNumber);  
  
  if ( split(binNumber, dim, splitPoint) == 1) return 1;
  if ( smartSplit(binNumber, dim, 0.5)   == 1) return 1;

  int ndim = _binningDimensions.size();
  for (int i = 0; i < ndim; i++){
    if ( smartSplit(binNumber, _binningDimensions.at(i), 0.5) == 1 ) return 1;
  }

  return 0;

}

/// \todo remember how this works
///
///
int HyperBinningMaker::smartLikelihoodSplit(int binNumber){
  
  int      dim = -1;
  double splitPoint = -1.0;

  getDimWithLargestSplitSignificance(dim, splitPoint, binNumber, false);

  if ( smartSplit(binNumber, dim, 0.5)   == 1) return 1;

  int ndim = _binningDimensions.size();
  for (int i = 0; i < ndim; i++){
    if ( smartSplit(binNumber, _binningDimensions.at(i), 0.5) == 1 ) return 1;
  }

  return 0;

}

/// \todo remember how this works
///
///
int HyperBinningMaker::likelihoodSplitAll(){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == 1) {
      nSplits += likelihoodSplit(i);
    }
  }

  return nSplits;
}

/// \todo remember how this works
///
///
int HyperBinningMaker::smartLikelihoodSplitAll(){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == 1) {
      smartLikelihoodSplit(i);
      nSplits++;
    }
  }

  return nSplits;
}


///for a specific volume number count how many HyperPoint's (or the sum of those
///HyperPoint's weights) are below the split value x_dim where
///
///x_dim = lowEdge_dim + (highEdge_dim - lowEdge_dim)*splitPoint
///
double HyperBinningMaker::countEventsBelowSplitPoint(int binNumber, int dimension, double splitPoint ) const{
  
  bool forceNoSnapToGrid = true;

  const HyperCuboid& cuboid = _hyperCuboids.at(binNumber);
  HyperCuboid newCuboid = splitBelowPoint(dimension, splitPoint, cuboid, forceNoSnapToGrid);
  
  const HyperPointSet& chosenHyperPointSet = _hyperPointSets.at(binNumber);

  double events = countEventsInHyperCuboid(chosenHyperPointSet, newCuboid);
  
  return events;
}

///for a specific volume number count how many HyperPoint's (or the sum of those
///HyperPoint's weights) are above the split value x_dim where
///
///x_dim = lowEdge_dim + (highEdge_dim - lowEdge_dim)*splitPoint
///
double HyperBinningMaker::countShadowEventsBelowSplitPoint(int binNumber, int dimension, double splitPoint) const{

  bool forceNoSnapToGrid = true;

  const HyperCuboid& cuboid = _hyperCuboids.at(binNumber);
  HyperCuboid newCuboid = splitBelowPoint(dimension, splitPoint, cuboid, forceNoSnapToGrid);
  
  const HyperPointSet& chosenHyperPointSet = _shadowHyperPointSets.at(binNumber);

  double events = countEventsInHyperCuboid(chosenHyperPointSet, newCuboid);
  
  return events;
}

///Count the number of HyperPoints (or the sum of their weights) within the HyperCuboid
///
double HyperBinningMaker::countEventsInHyperCuboid(const HyperPointSet& hyperPointSet, const HyperCuboid& hyperCuboid) const{

  double count = 0.0;
  for(unsigned int i = 0; i < hyperPointSet.size(); i++){
    if( hyperCuboid.inVolume(hyperPointSet.at(i)) ) {
       count = count + getWeight( hyperPointSet.at(i) ) ;
    }
  }
  return count;

}

///Find the split point such that a specified fraction of events falls within the
///resulting two bins.
double HyperBinningMaker::findSmartSplitPoint(int binNumber, int dimension, double dataFraction) const{
  
  double splitPoint = 0.5;
  double shift      = 0.25;

  double dataBefore = getSumOfWeights( _hyperPointSets.at(binNumber) );

  while (shift > 0.001){
    double dataNow = countEventsBelowSplitPoint(binNumber, dimension, splitPoint);
    double dataFrac = dataNow/dataBefore;

    if (dataFrac < dataFraction) splitPoint += shift;
    else if (dataFrac > dataFraction) splitPoint -= shift;
    else break;

    shift = shift*0.5;
  }
  
  return splitPoint;
}

///Same as findSmartSplitPoint, but only allow splits to be made halfway between 
///two integers. This assumes that all the HyperPoint elements of this dimension
///are also integers.
double HyperBinningMaker::findSmartSplitPointInt(int binNumber, int dimension, double dataFraction) const{
  
  const HyperCuboid& cuboid = _hyperCuboids.at(binNumber);

  HyperPoint lowCorner   ( cuboid.getLowCorner () );
  HyperPoint highCorner  ( cuboid.getHighCorner() );  

  double lowEdge   = lowCorner .at(dimension);
  double highEdge  = highCorner.at(dimension);

  int lowInt  = ceil (lowEdge );
  int highInt = floor(highEdge);
  
  if (lowInt == highInt) return 0.0; // everything is an integer - so if this is true, there is only one integer between the limits already

  VERBOSE_LOG << "------------------> bin num " << binNumber << ", dimension " << dimension << std::endl;
  VERBOSE_LOG << "------------------> [ " << lowEdge << ", " << highEdge << "]" << std::endl;
  VERBOSE_LOG << "------------------> [ " << lowInt  << ", " << highInt  << "]" << std::endl;

  double dataBefore = getSumOfWeights( _hyperPointSets.at(binNumber) );
  
  double splitPoint = 0.0;

  for (int i = lowInt; i < highInt; i++){

    double val = double(i) + 0.5; 
    splitPoint = (val - lowEdge)/(highEdge - lowEdge);

    double dataNow = countEventsBelowSplitPoint(binNumber, dimension, splitPoint);
    double dataFrac = dataNow/dataBefore;
    
    if (dataFrac > dataFraction) break;

  }
  
  VERBOSE_LOG << "------------------> split point " << splitPoint << std::endl;

  return splitPoint;
}


///split the bin in a chosen dimension to give a chosen fraction of events in the resulting bin
///
int HyperBinningMaker::smartSplit(int binNumber, int dimension, double dataFraction){

  double splitPoint = findSmartSplitPoint(binNumber, dimension, dataFraction);
  return split(binNumber, dimension, splitPoint);

}
 



///Split a bin into N parts, each which contain the same number of events
///
int HyperBinningMaker::smartMultiSplit(int binNumber, int dimension, int parts){
  
  int binToSplit = binNumber;
  
  int nSplits = 0;

  for (int i = 0; i < (parts - 1); i++){
    double fraction = 1.0/double(parts - i);
    double splitPoint = findSmartSplitPoint(binToSplit, dimension, fraction);
    bool ok = split(binToSplit, dimension, splitPoint);
    nSplits += ok;
    if (!ok) return nSplits;
  
    binToSplit = _linkedBins.at(binToSplit).at(1);
  }

  return nSplits;  
}

///Split a bin into N parts, each which contain the same number of events
///
int HyperBinningMaker::smartMultiSplit(int binNumber, int dimension){
  
  const HyperPointSet& points = _hyperPointSets.at(binNumber);
  
  double nEvents = getSumOfWeights(points);
  double ratio   = nEvents/_minimumBinContent;
  
  if (ratio >= 0.0 && ratio < 3.0) return smartMultiSplit(binNumber, dimension, 2);
  if (ratio >= 3.0 && ratio < 4.0) return smartMultiSplit(binNumber, dimension, 3);
  if (ratio >= 4.0 && ratio < 5.0) return smartMultiSplit(binNumber, dimension, 2);
  if (ratio >= 5.0 && ratio < 6.0) return smartMultiSplit(binNumber, dimension, 5);
  
  return smartMultiSplit(binNumber, dimension, 2);

}


///same as smartSplit but assumes all HyperPoint elements in the chosen dimension
///are integers. Split points are therefore always an integer + 0.5
int HyperBinningMaker::smartSplitInt(int binNumber, int dimension, double dataFraction){

  double splitPoint = findSmartSplitPointInt(binNumber, dimension, dataFraction);
  return split(binNumber, dimension, splitPoint);

}

///Split every volume with the CONTINUE status at the split point
///given (split point of 0.5 would split the bin into two equal parts)
///
int HyperBinningMaker::splitAll(int dimension, double splitPoint){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == VolumeStatus::CONTINUE) {
      nSplits += split(i, dimension, splitPoint);
    }
  }

  return nSplits;
}

/// Split every volume with the CONTINUE status using 
/// the smartSplit function. This makes one of the resulting bins
/// have a fraction 'dataFraction' of the original events.
int HyperBinningMaker::smartSplitAll(int dimension, double dataFraction){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == VolumeStatus::CONTINUE) {
      nSplits += smartSplit(i, dimension, dataFraction);
    }
  }


  return nSplits;
}

/// Split every volume with the CONTINUE status using 
/// the smartSplit function. This makes one of the resulting bins
/// have a fraction 'dataFraction' of the original events.
int HyperBinningMaker::smartMultiSplitAll(int dimension){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == VolumeStatus::CONTINUE) {
      nSplits += smartMultiSplit(i, dimension);
    }
  }


  return nSplits;
}



/// Split every volume with the CONTINUE status at the using 
/// the smartSplit function. This makes one of the resulting bins
/// have a fraction 'dataFraction' of the original events.
///
/// This assumes that the HyperPoint elements of this dimension
/// are integers, so splits are only made at integers + 0.5
int HyperBinningMaker::smartSplitAllInt(int dimension, double dataFraction){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == 1) {
      nSplits += smartSplitInt(i, dimension, dataFraction);
    }
  }

  return nSplits;
}

/// Split every volume with the CONTINUE status using 
/// the smartSplit function. This makes one of the resulting bins
/// have a fraction 'dataFraction' of the original events.
/// The dimension to split in is chosen at random.
int HyperBinningMaker::smartSplitAllRandomise(double dataFraction){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;
  
  int ndim = _binningDimensions.size();

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == VolumeStatus::CONTINUE) {
      int dimension = floor(_random->Uniform(0, ndim));
      nSplits += smartSplit(i, _binningDimensions.at(dimension) , dataFraction);
    }
  }

  return nSplits;

}

///Split every volume with the CONTINUE status at the split point
///given (split point of 0.5 would split the bin into two equal parts)
/// The dimension to split in is chosen at random.
int HyperBinningMaker::splitAllRandomise(double splitPoint){
  
  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  int ndim = _binningDimensions.size();

  for (int i = 0; i < initialSize; i++){
    if (_status.at(i) == VolumeStatus::CONTINUE) {
      int dimension = floor(_random->Uniform(0, ndim));
      nSplits += split(i, _binningDimensions.at(dimension), splitPoint);
    }
  }

  return nSplits;

}

///Set the minimum bin content allowed in any bin 
/// (calcualted from the HyperPoints)
void HyperBinningMaker::setMinimumBinContent(double val){
  _minimumBinContent = val;
}


///Set the HyperFunction - only used by some binning Algs
void HyperBinningMaker::setHyperFunction(HyperFunction* fnc){
  _func = fnc;
} 


///Set the minimum bin content allowed in any bin 
/// (calcualted from the shadow HyperPoints)
void HyperBinningMaker::setShadowMinimumBinContent(double val){
  _shadowMinimumBinContent = val;
}

///Set the minimum edge length for all dimensions
///
void HyperBinningMaker::setMinimumEdgeLength(double val){

  _minimumEdgeLength = HyperPoint( _hyperCuboids.at(0).getDimension(), val );

}

///Set the minimum edge length for all dimensions.
///
void HyperBinningMaker::setMinimumEdgeLength(HyperPoint val)  { 

  _minimumEdgeLength = val; 

}

///Find the number of bins - this is NOT the number of
///HyperVolumes. Some HyperVolumes are stored as part of the binning
///hierarchy which help event to be binned efficiently - this is
///descirbed in detail in the HyperVolumeBinning class.
int HyperBinningMaker::getNumBins() const{
  int count = 0;
  for ( int i = 0; i < getNumHyperVolumes(); i++){
    if (_linkedBins.at(i).size() == 0) count++;
  }
  return count;
}

///return the number of HyperVolumes - this is NOT the number
///of bins
int HyperBinningMaker::getNumHyperVolumes() const{
  
  return _hyperCuboids.size();

}

///Find the number of bins with the global status of
///CONTINUE, and the dimension specific status CONTINUE.
///Note if the dimension given is -1, then only consider
///Global status
int HyperBinningMaker::getNumContinueBins(int dimension) const{
  
  int size  = getNumHyperVolumes();
  int count = 0;
  
  //std::cout << "getNumContinueBins " << dimension << std::endl;

  for (int i = 0; i < size; i++){
    if ( getGlobalVolumeStatus(i) == VolumeStatus::CONTINUE ) {
      if ( dimension == -1 ){
        count++;
      } 
      else{
        if ( getDimensionSpecificVolumeStatus(i, dimension) == VolumeStatus::CONTINUE ){
          count++;
        }
      }
    }
  }  
  
  return count;
}



///Use the current state of the HyperBinningMaker to create
///a HyperVolumeBinning
HyperBinningMemRes HyperBinningMaker::getHyperVolumeBinning() const{
  
  int dimension = _hyperCuboids.at(0).getDimension();

  HyperBinningMemRes temp;
  
  for (unsigned int i = 0; i < _hyperCuboids.size(); i++){
    HyperVolume hyperVolume(dimension);
    hyperVolume.addHyperCuboid(_hyperCuboids.at(i));
    temp.addHyperVolume(hyperVolume, _linkedBins.at(i));
  }

  temp.addPrimaryVolumeNumber(0);

  if (s_printBinning == true) INFO_LOG << "Made HyperVolumeBinning"<<std::endl;

  return temp;
 
}


///Use the current state of the HyperBinningMaker to create
///a HyperBinningHistogram
HyperHistogram* HyperBinningMaker::getHyperBinningHistogram() const{

  HyperBinningMemRes binning = getHyperVolumeBinning();
  
  HyperHistogram* histogram = new HyperHistogram( binning );

  for ( int i = 0; i < binning.getNumBins(); i++){
    int volumeNumber = binning.getHyperVolumeNumber(i);
    histogram->setBinContent(i, double( _hyperPointSets.at(volumeNumber).getSumW() ) );
    histogram->setBinError  (i, double( sqrt(_hyperPointSets.at(volumeNumber).getSumW2() ) ) );
  }
  
  if (s_printBinning == true) {
    INFO_LOG << "Made HyperHistogram"<<std::endl;
  }
  
  histogram->setNames(_names);

  return histogram;

}

///Use the current state of the HyperBinningMaker to create
///a HyperBinningHistogram for the shadow events
HyperHistogram* HyperBinningMaker::getShadowHyperBinningHistogram() const{

  HyperBinningMemRes binning = getHyperVolumeBinning();
  
  HyperHistogram* histogram = new HyperHistogram( binning );

  for ( int i = 0; i < binning.getNumBins(); i++){
    int volumeNumber = binning.getHyperVolumeNumber(i);
    histogram->setBinContent(i, double( _shadowHyperPointSets.at(volumeNumber).getSumW() ) );
    histogram->setBinError  (i, double( sqrt(_shadowHyperPointSets.at(volumeNumber).getSumW2() ) ) );
  }
  
  if (s_printBinning == true) INFO_LOG << "Made HyperHistogram"<<std::endl;
  
  histogram->setNames(_names);
  
  return histogram;

}

///Use the current state of the HyperBinningMaker to create
///a HyperBinningHistogram for the ratio of events to shadow events
HyperHistogram* HyperBinningMaker::getRatioHyperBinningHistogram() const{

  HyperHistogram* numerator = new HyperHistogram( *getHyperBinningHistogram() );
  HyperHistogram denominator( *getShadowHyperBinningHistogram() );

  numerator->divide(denominator);

  return numerator;

}

/** 
Use the current status of the HyperBinningMaker to make a 
HyperVolumeBinning, and draw it.
*/
void HyperBinningMaker::drawCurrentState(TString path) const{

  HyperHistogram* hist = getHyperBinningHistogram();
  
  hist->drawDensity(path);

}


/** 
Call this at the beginning of the algorithm
*/
void HyperBinningMaker::startedAlgorithm(){


}


/** 
Call this at the end of the algorithm

*/
void HyperBinningMaker::finishedAlgorithm(){

}


/** 
When making an algorithm, call this at the beginning of each iteration
if possible
*/
void HyperBinningMaker::startedIteration(){



}

/** 
When making an algorithm, call this after each iteration
if possible
*/
void HyperBinningMaker::finishedIteration(){

  if (_drawAlgorithm){

    TString outputdir = _drawAlgorithmDir + "_iteration";
    outputdir += _iterationNum;

    drawCurrentState(outputdir);

  } 

  _iterationNum++;

}

/** 
Use this option if you want the binning to be drawn after every 
iteration of the binning algorithm.
*/
void HyperBinningMaker::drawAfterEachIteration(TString path){
  _drawAlgorithm = true;
  _drawAlgorithmDir = path;
}


bool HyperBinningMaker::s_printBinning = true;

///destructor
HyperBinningMaker::~HyperBinningMaker(){
  delete _random;
}



