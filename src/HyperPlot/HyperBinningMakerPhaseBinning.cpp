#include "AmpGen/HyperPlot/HyperBinningMakerPhaseBinning.h"







HyperBinningMakerPhaseBinning::HyperBinningMakerPhaseBinning(const HyperCuboid& binningRange, HyperFunction* func) :
  HyperBinningMaker(binningRange, HyperPointSet( binningRange.getDimension() )),
  _numBinPairs     (0),
  _maximumRandWalks(30),
  _numWalkers      (5),
  _walkSizeFrac  (0.12),
  _numberOfSystematicSplits(0),
  _numberOfGradientSplits  (0)
{
  setNumBinPairs(3);
  setHyperFunction(func);
  WELCOME_LOG << "Good day from the HyperBinningMakerPhaseBinning() Constructor"<<std::endl; 
}

void HyperBinningMakerPhaseBinning::makeBinning(){
   
  int dimension = _binningDimensions.size();
  
  int splitDim = 0;
  if (splitDim >= dimension) splitDim = 0;

  int nBins = 0;
  int unchanged = 0;

 
  while ( 1==1 ){
    finishedIteration();

    if (s_printBinning == true) {
      INFO_LOG << "--------------------------------------------------" <<std::endl;
      INFO_LOG << "Trying to split all bins" << std::endl;
      INFO_LOG << "There are " << getNumContinueBins() << " splittable bins" << std::endl;
    }

    gradientSplitAll();

    if (nBins == getNumBins()) unchanged++;
    else unchanged = 0;
    if (unchanged >= dimension) break;
    nBins = getNumBins();
    if (s_printBinning == true) INFO_LOG << "There is now a total of " << nBins << " bins" << std::endl;

  }

  if (s_printBinning == true) INFO_LOG << "Mint binning algorithm complete " <<std::endl;

} 


int HyperBinningMakerPhaseBinning::gradientSplitAll(){
  
  VERBOSE_LOG <<  "HyperBinningMakerPhaseBinning::gradientSplitAll( )" << std::endl;

  int initialSize = _hyperCuboids.size();
  int nSplits = 0;

  int splittable = getNumContinueBins();
  LoadingBar loadingbar(splittable);
  int splittableDone = 0;

  for (int i = 0; i < initialSize; i++){
    if (getGlobalVolumeStatus(i) == VolumeStatus::CONTINUE) {
           
      loadingbar.update(splittableDone);
      splittableDone++;
      
      int dimension = -1;
      int split = gradientSplit( i, dimension );
      
      nSplits += split;

      if (split == 0){
        getDimensionSpecificVolumeStatus(i, dimension) = VolumeStatus::DONE;
        updateGlobalStatusFromDimSpecific(i);
      }

    }
  }

  return nSplits;
}


void HyperBinningMakerPhaseBinning::walkOrthogonal(HyperPoint& point, HyperCuboid& walkLimits){

  //randomly select what direction we're going to walk in
  int walkDimNum = floor(_random->Uniform(0.0, _binningDimensions.size()));
  int walkDim    = _binningDimensions.at(walkDimNum);
  
  double walkSize = _walkSizeFrac*( walkLimits.getHighCorner().at(walkDim) - walkLimits.getLowCorner().at(walkDim) );
  //if (walkSize < _minimumEdgeLength.at(walkDim)) walkSize = _minimumEdgeLength.at(walkDim);
  
  //randomly move along this axis
  double direction = _random->Gaus(0.0, 1.0);
  point.at(walkDim) += walkSize*direction;
  
  //If this walk went out of the walk limits, move back, and try again.
  if ( walkLimits.inVolume(point) == false ){
    point.at(walkDim) -= walkSize*direction;
    walk( point, walkLimits);
  }
  
}

void HyperBinningMakerPhaseBinning::walk(HyperPoint& point, HyperCuboid& walkLimits){
  
  HyperPoint walkWidth = (walkLimits.getHighCorner() - walkLimits.getLowCorner())*_walkSizeFrac;
  
  for (int i = 0; i < point.getDimension(); i++){
    if ( isValidBinningDimension(i) == false ) {
      walkWidth.at(i) = 0.0;
    }
    else{
      walkWidth.at(i) = walkWidth.at(i)*_random->Gaus(0.0, 1.0);
    }
  } 

  point = point + walkWidth;
  
  //If this walk went out of the walk limits, move back, and try again.
  if ( walkLimits.inVolume(point) == false ){
    point = point - walkWidth;
    walk( point, walkLimits);
  }
  
}

int HyperBinningMakerPhaseBinning::splitByCoord(int volumeNumber, int dimension, HyperPoint& point){

  double low  = _hyperCuboids.at(volumeNumber).getLowCorner ().at(dimension);
  double high = _hyperCuboids.at(volumeNumber).getHighCorner().at(dimension);      
  double splitPoint = (point.at(dimension) - low)/(high-low);
  return split(volumeNumber, dimension, splitPoint);

}


void HyperBinningMakerPhaseBinning::setNumBinPairs(int binpairs){
  _numBinPairs = binpairs;
  std::vector<double> binEdges;

  for (int i = 0; i < 2*_numBinPairs; i++){
    binEdges.push_back( - TMath::Pi() + (TMath::Pi()/double(_numBinPairs))*i );
  }

  _binEdges = CyclicPhaseBins(binEdges);
} 
  
void HyperBinningMakerPhaseBinning::setBinEdges(std::vector<double> binEdges){ 
  _binEdges    = CyclicPhaseBins(binEdges); 
  _numBinPairs = binEdges.size(); 
} 


int HyperBinningMakerPhaseBinning::getBinNumFromFuncVal(double phase){

  if (phase < -10.0) return -1;
  
  return _binEdges.getBinNumber(phase);

}

double HyperBinningMakerPhaseBinning::getLowBinBoundary(double phase){
  
  return _binEdges.getLowBinBoundary(phase);

}

double HyperBinningMakerPhaseBinning::getHighBinBoundary(double phase){
  
  return _binEdges.getHighBinBoundary(phase);

}

double HyperBinningMakerPhaseBinning::closestBinBoundary(double val){

  double low  = getLowBinBoundary (val);
  double high = getHighBinBoundary(val);

  double distLow  = fabs(val - low );
  double distHigh = fabs(val - high);

  if (distLow < distHigh) return low;

  return high;

}

int HyperBinningMakerPhaseBinning::getBinNumFromFunc(HyperPoint& point) {
  
  return getBinNumFromFuncVal(_func->getVal(point));
  
}

HyperPoint HyperBinningMakerPhaseBinning::getGrad(HyperPoint& point){

  int nBinningDims = _binningDimensions.size();
  int nDims = _hyperCuboids.at(0).getDimension();
  
  HyperPoint grad(nDims, 0.0);

  for (int i = 0; i < nBinningDims; i++){
    int dim = _binningDimensions.at(i);
    double stepsize = _minimumEdgeLength.at(dim)*0.5;
    
    HyperPoint low (point);
    HyperPoint high(point);
    low .at(i) -= stepsize;
    high.at(i) += stepsize;
    
    double valLow  = _func->getVal(low );
    double valHigh = _func->getVal(high);
    
    std::complex<double> compLow (cos(valLow ), sin(valLow ));
    std::complex<double> compHigh(cos(valHigh), sin(valHigh));

    grad.at(i) = arg( compHigh*conj(compLow) )/(2.0*stepsize);
  }
  return grad;
}

HyperPoint HyperBinningMakerPhaseBinning::getGradPos(HyperPoint& point, double funcValAtPoint){

  int nBinningDims = _binningDimensions.size();
  int nDims = _hyperCuboids.at(0).getDimension();
  
  HyperPoint grad(nDims, 0.0);

  for (int i = 0; i < nBinningDims; i++){
    int dim = _binningDimensions.at(i);
    double stepsize = _minimumEdgeLength.at(dim)*0.5;
    
    HyperPoint high(point);
    high.at(i) += stepsize;
    
    double valHigh = _func->getVal(high);
    
    std::complex<double> compVal (cos(funcValAtPoint ), sin(funcValAtPoint ));
    std::complex<double> compHigh(cos(valHigh), sin(valHigh));

    grad.at(i) = arg( compHigh*conj(compVal) )/(stepsize);
  }

  return grad;
}

double HyperBinningMakerPhaseBinning::getSecondDerivative(HyperPoint& point, HyperPoint& vector, double funcValAtPoint, double& firstDeriv){
  
  bool centeral = true;
  
  HyperPoint normVector = vector/vector.norm();

  double stepLength = fabs(_minimumEdgeLength.dotProduct(normVector)*0.5);
  //double stepLength = vector.norm()*0.25;

  if (centeral){
    HyperPoint low  = point - normVector*stepLength; 
    HyperPoint high = point + normVector*stepLength; 
  
    double val     = funcValAtPoint;
    double lowVal  = _func->getVal(low  );
    double highVal = _func->getVal(high );
  
    double h = (high - low).norm()*0.5;
    
    firstDeriv = (highVal - lowVal)/(2.0*h);
    double secondDeriv = (lowVal + highVal - 2.0*val)/(h*h);
    return secondDeriv;
  }

  HyperPoint high1 = point + normVector*stepLength; 
  HyperPoint high2 = point + normVector*stepLength*2.0; 
  
  double val      = _func->getVal(point  );
  double high1Val = _func->getVal(high1  );
  double high2Val = _func->getVal(high2  );
  
  double h = (high1 - high2).norm();
    
  double secondDeriv = (val + high2Val - 2.0*high1Val)/(h*h);
  return secondDeriv;
  

}


int HyperBinningMakerPhaseBinning::splitDimFromGrad(int volumeNumber, HyperPoint gradient){
  
  double minGrad = gradient.norm()*0.25;


  double maxGrad        = -1.0;
  int    dimWithMaxGrad = -1;

  for (int i = 0; i < gradient.getDimension(); i++){
    
    if ( getDimensionSpecificVolumeStatus(volumeNumber, i) == VolumeStatus::DONE ) continue;

    double gradInThisDim = fabs( gradient.at(i) ) + minGrad;
    
    gradInThisDim *= _hyperCuboids.at(volumeNumber).getWidth(i);
    gradInThisDim /= _minimumEdgeLength.at(i);

    if ( gradInThisDim > maxGrad ) {
      maxGrad        = gradInThisDim;
      dimWithMaxGrad = i;
    }

  }
  
  if (dimWithMaxGrad == -1){
    ERROR_LOG << "Something has gone horribly wrong - if no dimensions can be split, I shouldn't be arriving here" << std::endl;
  }

  //std::cout << gradient << "  " << dimWithMaxGrad <<std::endl;

  return dimWithMaxGrad;
}

HyperPointSet HyperBinningMakerPhaseBinning::getSplitCorners( int volumeNumber ){
  
  //Get the cuboid we are trying to split  
  HyperCuboid&   chosenHyperCuboid = _hyperCuboids.at(volumeNumber);

  //Figure out what dimesions we need
  int dimensions     = _hyperCuboids.at(0).getDimension();

  std::vector<int> binningDimensions;
  for (int i = 0; i < dimensions; i++){
    if ( chosenHyperCuboid.getWidth(i) < 2.0*_minimumEdgeLength.at(i) ) continue;
    if ( isValidBinningDimension(i) == false ) continue;
    binningDimensions.push_back(i);
  }

  int nSplittingDims = binningDimensions.size();

  //Get the center of the bin and get the function / bin value
  HyperPoint point = chosenHyperCuboid.getCenter();
  
  //Project the cuboid into the splitting dimensions
  HyperCuboid projectedCuboid = chosenHyperCuboid.project(binningDimensions);
  
  VERBOSE_LOG << "Geting vector of all my verticies" << std::endl;

  //Get the verticies of the projected cuboid
  HyperPointSet projectedVerticies = projectedCuboid.getVertices();

  //Turn the `projected verticies` into `verticies` that use the coordinates
  //of the bin centre for all dimensions that are not in binningDimensions

  HyperPointSet verticies(dimensions);

  for (unsigned i = 0; i < projectedVerticies.size(); i++){
    HyperPoint vertex(point);
    
    for (int j = 0; j < nSplittingDims; j++){
      vertex.at( binningDimensions.at(j) ) = projectedVerticies.at(i).at(j);
    }
    
    //move the verticies a very small amount towards the bin centre.
    //can cause issues if the verticies are on the limit of the 
    //kinematically allowed region

    for (int j = 0; j < nSplittingDims; j++){
       
      int thisDim = binningDimensions.at(j);

      double addOrSubtract = 0.0;
      double pointVal = point .at( thisDim );
      double vertVal  = vertex.at( thisDim );
      if ( vertVal > pointVal ) { addOrSubtract = -1.00; }
      else                      { addOrSubtract = +1.00; }
      
      //if (dimension != thisDim  ) {
      //  addOrSubtract *= 0.99;
      //}

      vertex.at( thisDim ) += addOrSubtract*_minimumEdgeLength.at( thisDim );
    }

    verticies.push_back(vertex);
  }

  return verticies;
}

HyperPointSet HyperBinningMakerPhaseBinning::getSplitFaces( int volumeNumber ){
  

  //Get the cuboid we are trying to split  
  HyperCuboid&   chosenHyperCuboid = _hyperCuboids.at(volumeNumber);

  //Figure out what dimesions we need
  int dimensions     = _hyperCuboids.at(0).getDimension();

  std::vector<int> binningDimensions;
  for (int i = 0; i < dimensions; i++){
    if ( chosenHyperCuboid.getWidth(i) < 2.0*_minimumEdgeLength.at(i) ) continue;
    if ( isValidBinningDimension(i) == false ) continue;
    binningDimensions.push_back(i);
  }

  int nSplittingDims = binningDimensions.size();


  //Get the center of the bin and get the function / bin value
  HyperPoint point = chosenHyperCuboid.getCenter();

  //Project the cuboid into the splitting dimensions
  HyperCuboid projectedCuboid = chosenHyperCuboid.project(binningDimensions);
  
  VERBOSE_LOG << "Geting vector of all my faces" << std::endl;



  HyperPointSet verticies(dimensions);

  for (int i = 0; i < nSplittingDims; i++){

    int thisDim = binningDimensions.at(i);
    HyperPoint faceA(point);
    HyperPoint faceB(point);

    faceA.at(thisDim) = faceA.at(thisDim) - 0.5*chosenHyperCuboid.getWidth(thisDim) + _minimumEdgeLength.at(thisDim);
    faceB.at(thisDim) = faceB.at(thisDim) + 0.5*chosenHyperCuboid.getWidth(thisDim) - _minimumEdgeLength.at(thisDim);
    
    verticies.push_back(faceA);
    verticies.push_back(faceB);

  }

  return verticies;

}

HyperPointSet HyperBinningMakerPhaseBinning::getSplitEdges( int volumeNumber ){
  

  //Get the cuboid we are trying to split  
  HyperCuboid&   chosenHyperCuboid = _hyperCuboids.at(volumeNumber);

  //Figure out what dimesions we need
  int dimensions     = _hyperCuboids.at(0).getDimension();

  std::vector<int> binningDimensions;
  for (int i = 0; i < dimensions; i++){
    if ( chosenHyperCuboid.getWidth(i) < 2.0*_minimumEdgeLength.at(i) ) continue;
    if ( isValidBinningDimension(i) == false ) continue;
    binningDimensions.push_back(i);
  }

  int nSplittingDims = binningDimensions.size();

  if (nSplittingDims == 0){
    ERROR_LOG << "There should always be at least 1 splittable dimension!" << std::endl;
  }

  //Get the center of the bin and get the function / bin value
  HyperPoint point = chosenHyperCuboid.getCenter();
  
  //Project the cuboid into the splitting dimensions
  HyperCuboid projectedCuboid = chosenHyperCuboid.project(binningDimensions);
  
  VERBOSE_LOG << "Geting vector of all my edges" << std::endl;
  VERBOSE_LOG << "There are " << nSplittingDims << " split dimensions" << std::endl;

  //Get the verticies of the projected cuboid
  HyperPointSet projectedVerticies = projectedCuboid.getEdgeCenters();

  VERBOSE_LOG << "There are " << projectedVerticies.size() << " verticies" << std::endl;

  //Turn the `projected verticies` into `verticies` that use the coordinates
  //of the bin centre for all dimensions that are not in binningDimensions

  HyperPointSet verticies(dimensions);

  for (unsigned i = 0; i < projectedVerticies.size(); i++){
    HyperPoint vertex(point);
    
    for (int j = 0; j < nSplittingDims; j++){
      vertex.at( binningDimensions.at(j) ) = projectedVerticies.at(i).at(j);
    }
    
    //move the verticies a very small amount towards the bin centre.
    //can cause issues if the verticies are on the limit of the 
    //kinematically allowed region

    for (int j = 0; j < nSplittingDims; j++){
       
      int thisDim = binningDimensions.at(j);

      double addOrSubtract = 0.0;
      double pointVal = point .at( thisDim );
      double vertVal  = vertex.at( thisDim );
      if      ( vertVal == pointVal ) { addOrSubtract =  0.00; }
      else if ( vertVal >  pointVal ) { addOrSubtract = -1.00; }
      else                            { addOrSubtract = +1.00; }
      
      //if (dimension != thisDim  ) {
      //  addOrSubtract *= 0.99;
      //}

      vertex.at( thisDim ) += addOrSubtract*_minimumEdgeLength.at( thisDim );
    }

    verticies.push_back(vertex);
  }

  VERBOSE_LOG << "About to return my verticies" << std::endl;


  return verticies;

}



HyperPoint HyperBinningMakerPhaseBinning::orderAndTestSplitPoints(HyperPointSet& points, HyperPoint& point, double valAtPoint, HyperPoint gradient){

  int    bin = getBinNumFromFuncVal(valAtPoint);

  //To speed things up, we want to order the corners by 
  //the most likely to split. Can use gradient to see
  
  //The FOM is essentially just the distance to the nearest bin edge.
  //It's negative if we are still inside the bin.

  std::vector<double> functionEstimateFom;
  std::vector<double> functionEstimate;
  std::vector<int>    index;
  
  double lowBoundary  = getLowBinBoundary (valAtPoint);
  double highBoundary = getHighBinBoundary(valAtPoint);
  
  VERBOSE_LOG << "Sorting points to speed my life up" << std::endl;

  for (unsigned i = 0; i < points.size(); i++){

    index.push_back(i);

    double estimate = valAtPoint + gradient.dotProduct(points.at(i) - point);
    double fom = 0.0;
    if ( estimate >= lowBoundary && estimate <= highBoundary ){
      double fom1 = lowBoundary - estimate;
      double fom2 = estimate - highBoundary; 

      if (fom1 > fom2)  {fom = fom1;}
      else              {fom = fom2;}  
    }
    if (estimate < lowBoundary){
      fom = lowBoundary - estimate;
    }
    if (estimate > highBoundary){
      fom = estimate - highBoundary;
    }
    
    functionEstimate   .push_back(estimate);
    functionEstimateFom.push_back(fom     );

  }
  

  TMath::Sort( (int)points.size(), &(functionEstimateFom.at(0)), &(index.at(0)) );
  //for (unsigned i = 0; i < points.size(); i++){
  //  std::cout << functionEstimateFom.at( index.at(i) ) << std::endl;
  //}

  //Evaluate the function for each vertex and see if the bin number changes. 
  //If so, choose a split point that is somewhere between the vetex and the
  //bin center.
  
  //Each time we evaluate the function at a corner, we can see how good
  //our predictions are, and set an error.

  double uncertaintySum = 0.0;
  int uncertainiesAdded = 0;
    
  for (unsigned i = 0; i < points.size(); i++){

    double vtxNum = index.at(i);
    double estimate = functionEstimate   .at(vtxNum);
    double fom      = functionEstimateFom.at(vtxNum);    
    double estimateUncert = TMath::Pi()*2.0;

    if (uncertainiesAdded >= 2){
      estimateUncert = uncertaintySum/double(uncertainiesAdded - 1.0);
    }
    
    // i.e. if our estimate indicates that this bin is
    // more than 3.5 sigma from a bin boundary, we give up.
    if ( fom < -estimateUncert*5.0 ) {
      //std::cout << "Breaking loop at " << i << "; fom est = " <<  fom << " ± " << estimateUncert << std::endl;
      break; 
    }
    double valHere = _func->getVal(points.at(vtxNum));
    int binHere    = getBinNumFromFuncVal(valHere);
    
    std::complex<double> cEst(cos(estimate), sin(estimate));
    std::complex<double> cValHere(cos(valHere) , sin(valHere ));
    
    double diff = arg( cEst*conj(cValHere) );
    
    uncertaintySum += fabs(diff);
    uncertainiesAdded++;

    //std::cout << i << "   "<< estimate << " : [ " << lowBoundary << ", " << highBoundary <<  "] "<< valHere << std::endl;
    
    if (binHere != bin){
      HyperPoint returnPoint(points.at(vtxNum));
      returnPoint.setWeight(0, valHere);
      return returnPoint;
    }

  }
  
  HyperPoint returnPoint(0);
  if (returnPoint.getDimension() !=0 ){
    ERROR_LOG << "This should be zero" << std::endl;
  }

  return returnPoint;

}

int HyperBinningMakerPhaseBinning::systematicSplit(int volumeNumber, int dimension, double valAtCenter, HyperPoint gradient){
  
  _numberOfSystematicSplits++;
  
  //int binNumAtCenter = getBinNumFromFuncVal(valAtCenter);

  HyperCuboid&   chosenHyperCuboid = _hyperCuboids.at(volumeNumber);
  HyperPoint centerPoint    = chosenHyperCuboid.getCenter();


  
  HyperPointSet pointsToTest( _hyperCuboids.at(0).getDimension() );

  HyperPointSet corners = getSplitCorners(volumeNumber);
  HyperPointSet faces   = getSplitFaces  (volumeNumber);
  HyperPointSet edges   = getSplitEdges  (volumeNumber);
  
  pointsToTest.addHyperPointSet(corners);
  pointsToTest.addHyperPointSet(faces  );
  pointsToTest.addHyperPointSet(edges  );
  
  VERBOSE_LOG << "Got a big list of " << pointsToTest.size() << " to test " << std::endl;

  //pointsToTest.print();

  //for (int i = 0; i < corners.size(); i++){
  //  pointsToTest.push_back( (corners.at(i) + centerPoint)*0.5 );
  //}


  HyperPoint point          = orderAndTestSplitPoints(pointsToTest, centerPoint, valAtCenter, gradient);

  if ( point.getDimension() == 0 ) return 0;
  
  double valAtPoint = point.getWeight();

  std::complex<double> cValAtPoint (cos(valAtPoint)  , sin(valAtPoint ));
  std::complex<double> cValAtCenter(cos(valAtCenter) , sin(valAtCenter));

  double dy = arg( cValAtCenter*conj(cValAtPoint) );
      
  double binEdge = 0.0;
  if (dy > 0.0){
    binEdge = getLowBinBoundary(valAtCenter);
  }
  else{
    binEdge = getHighBinBoundary(valAtCenter);
  }

  double scale = fabs( (binEdge - valAtCenter)/dy );
  if (scale > 1.00){
    //ERROR_LOG << "THIS HAS TO BE LESS THAN 1.0 - It's " << scale << std::endl;
    scale = 1.0;
  }

  HyperPoint splitPoint = centerPoint + ( point - centerPoint )*scale;
  return splitByCoord(volumeNumber, dimension, splitPoint );


}






int HyperBinningMakerPhaseBinning::gradientSplit(int volumeNumber, int& dimension){
  
  _numberOfGradientSplits++;

  VERBOSE_LOG <<"Calling gradientSplit(" << volumeNumber << ", " << dimension << ")" << std::endl;

  //first get the Hypercube that we want to split

  HyperCuboid&   chosenHyperCuboid = _hyperCuboids.at(volumeNumber);


  //Get the center of the bin. Ideally we want
  //to find a split point as colse as possible to this
  HyperPoint binCenter = chosenHyperCuboid.getCenter();
  
  //Evaluate the function at the bin center and find the
  //'bin number' from this. Also find out which bin boundary
  //is closest to us. 

  double valCenter    = _func->getVal(binCenter);
  int    binNumber    = getBinNumFromFuncVal(valCenter);
  double closestBound = closestBinBoundary(valCenter);
  
  //Gradient given in the same coordinates as the binning  
  HyperPoint gradient = getGradPos(binCenter, valCenter);
  
  dimension = splitDimFromGrad(volumeNumber, gradient);
  
  //if (fabs(gradient.at(dimension))*10.0 < gradient.norm()){
  //  return randomWalkSplit(volumeNumber, dimension);
  //}

  //When we are splitting the bin, there is no point in splitting in a
  //region that will result in bins that are too narrow. Therefore define the "splitLimits" 

  HyperPoint     lowPoint          = chosenHyperCuboid.getLowCorner ();
  HyperPoint     highPoint         = chosenHyperCuboid.getHighCorner();
  
  lowPoint  = lowPoint  + _minimumEdgeLength*0.499;
  highPoint = highPoint - _minimumEdgeLength*0.499;

  lowPoint .at(dimension) += _minimumEdgeLength.at(dimension)*0.5;
  highPoint.at(dimension) -= _minimumEdgeLength.at(dimension)*0.5;
  
  HyperCuboid splitLimits(lowPoint, highPoint);
  HyperPoint  splitLimitsWidth = highPoint - lowPoint;
  VERBOSE_LOG << "I have a nice set of split limits" << std::endl;


  //tranform gradient to a coordinate system where the bin 
  //is limited by [-1, 1] in all dims
  
  HyperPoint gradientTrans = gradient;

  int nDims = _hyperCuboids.at(0).getDimension();

  for (int i = 0; i < nDims; i++){
    gradientTrans.at(i) *= (splitLimitsWidth.at(i)*2.0);
  }
  
  //if I move along in space by 'gradient' the change in the function will be

  double rateOfInc = gradientTrans.norm2();
  
  //We need the function to change by the distance to the nearest boundary: 

  double scale = (closestBound - valCenter)/rateOfInc;

  HyperPoint splitVector = gradientTrans*scale;

  //Now we need to transform the SplitVector back to 
  //the real coordinate system. 
  //The fact that we multiply by the same thing again
  //I found counter-intuative, but I think it's correct.

  for (int i = 0; i < nDims; i++){
    splitVector.at(i) *= (splitLimitsWidth.at(i)*2.0);
  }
  
  //Now we have a direction to look along, we can find the 
  //second derivate in this direction to improve the estimate

  //double firstDeriv  = gradientOrig.dotProduct(splitVector)/splitVector.norm();
  double firstDeriv  = 0.0;
  double secondDeriv = getSecondDerivative(binCenter, splitVector, valCenter, firstDeriv);
  
  //   f(x) = f(x0) + f'(x0) (x-x0) + f''(x0) (x-x0)^2 / 2  
  //solve the quadratic eqn

  double a = secondDeriv/2.0;
  double b = firstDeriv;
  double c = valCenter - closestBound;
  
  double quadraticSol1 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
  double quadraticSol2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
  
  //choose the solution that is closest to the bin center

  double quadraticSol = quadraticSol1;
  if ( fabs(quadraticSol1) > fabs(quadraticSol2) ){
    quadraticSol   = quadraticSol2;
  }
  
  //If the solution to the quadratic is NaN, just go with the
  //more simple cornerSplit function.

  if (quadraticSol != quadraticSol){
    return systematicSplit(volumeNumber, dimension, valCenter, gradient);
  }
  
  //Get a refined split vector using improved first deriv and second deriv
  HyperPoint refinedSplitVector = (splitVector * quadraticSol  ) / splitVector.norm();

  //Use the split vectors to get split points 

  HyperPoint splitPoint        = splitVector        + binCenter;
  HyperPoint refinedSplitPoint = refinedSplitVector + binCenter;
  
  //Want to add a small amount extra to make sure we cross a bin boundary.

  HyperPoint refinedSplitPointExtended = refinedSplitVector*1.05 + binCenter;

  //As a cross check, can make a function in the original coordinate
  //system that models the function using f(x0) and grad(x0). We then
  //evaluate this at the split point that was determined. Hopefully
  //this returns the nearest bin boundary!!! (It does)

  bool crossCheck1 = false;
  bool crossCheck2 = false;
  bool crossCheck4 = false;

  if (crossCheck1 == true){

    double predictedVal = valCenter + gradient.dotProduct( splitPoint - binCenter );
    std::cout << predictedVal << ", " << closestBound << std::endl;

  }
  
  if (crossCheck2 == true){

    double directionVal = splitVector.dotProduct(refinedSplitPoint - binCenter);
    
    if (directionVal < 0.0) directionVal = -1.0;
    else                    directionVal =  1.0;

    double predictedVal = 0.5*secondDeriv*(refinedSplitPoint - binCenter).norm2() + directionVal*firstDeriv*(refinedSplitPoint - binCenter).norm() + valCenter;

    std::cout << predictedVal << ", " << closestBound << std::endl;

  }
  
  if (crossCheck4 == true){

    std::cout << gradient.dotProduct( splitVector )/splitVector.norm() << ",  " << firstDeriv << std::endl;
  }


  //See if the predicted split point is within the bin limits, if not
  //return try a less sophisticated method.

  if (splitLimits.inVolume(refinedSplitPointExtended) == false) {
    return systematicSplit(volumeNumber, dimension, valCenter, gradient);
  }

  //See what the function and the bin number is at the new point
  double valNew    = _func->getVal(refinedSplitPointExtended);
  int    binNew    = getBinNumFromFuncVal(valNew);    
 
  //As a cross check we can see how the extrapolation agrees with the real value.

  bool crossCheck3 = false;
  
  if (crossCheck3 == true){

    double valReal     = _func->getVal(       splitPoint);
    double valRealRef  = _func->getVal(refinedSplitPoint);
    double valExtrap   = closestBound;

    std::complex<double> cValReal   (cos(valReal)   ,  sin(valReal)  );
    std::complex<double> cValRealRef(cos(valRealRef),  sin(valRealRef)  );    
    std::complex<double> cValExtrap (cos(valExtrap) ,  sin(valExtrap));
    double diff    = arg(cValReal   *conj(cValExtrap));
    double diffRef = arg(cValRealRef*conj(cValExtrap));

    std::cout << valCenter << " -> (" << valReal << " ?= " << valExtrap << " ) -> " << diff / fabs(valCenter-valExtrap) << std::endl;
    std::cout << diff << " -> " << diffRef << std::endl;

  }

  //If we end up in the same bin, try a less sophisticated approach

  if ( binNew == binNumber ) {
    return systematicSplit(volumeNumber, dimension, valCenter, gradient);
  }

  //Finally, if everything else goes OK, split bin at the chosen split point

  return splitByCoord(volumeNumber, dimension, refinedSplitPointExtended );

}

int HyperBinningMakerPhaseBinning::randomWalkSplit(int volumeNumber, int dimension){
  
  VERBOSE_LOG <<"Calling gradientSplit(" << volumeNumber << ", " << dimension << ")" << std::endl;

  //first get the Hypercube that we want to split

  HyperCuboid&   chosenHyperCuboid = _hyperCuboids.at(volumeNumber);

  HyperPoint     lowPoint          = chosenHyperCuboid.getLowCorner ();
  HyperPoint     highPoint         = chosenHyperCuboid.getHighCorner();
  

  //get the minumum bin width in the dimension we're splitting in
  double minBinWidth = _minimumEdgeLength.at(dimension);

  //if the volume we're splitting is narrower than minBinWidth*2 we bail now
  //if ( binWidth < minBinWidth*2.0 ) return 0;


  
  //When we do the random walk, there is no point in walking to regions that will
  //result in bins that are too narrow. Therefore define the "walkLimits" 

  
  HyperPoint binCenter = chosenHyperCuboid.getCenter();
  
  lowPoint  = lowPoint  + _minimumEdgeLength*0.499;
  highPoint = highPoint - _minimumEdgeLength*0.499;

  //for ( int i = 0; i < lowPoint.getDimension(); i++ ){
  //  if ( highPoint.at(i) - lowPoint.at(i) < 2.0*_minimumEdgeLength.at(i) ){
  //    lowPoint .at(i) = binCenter.at(i) - 0.01*_minimumEdgeLength.at(i);
  //    highPoint.at(i) = binCenter.at(i) + 0.01*_minimumEdgeLength.at(i);      
  //  }
  //  else{
  //    lowPoint .at(i) += _minimumEdgeLength.at(i);
  //    highPoint.at(i) -= _minimumEdgeLength.at(i);
  //  }    
  //}

  lowPoint .at(dimension) += minBinWidth*0.5;
  highPoint.at(dimension) -= minBinWidth*0.5;

  VERBOSE_LOG <<"Walk width in split dim = " << highPoint.at(dimension) - lowPoint.at(dimension) << std::endl;

  //INFO_LOG << "I bet this happens first" << std::endl;
  HyperCuboid walkLimits(lowPoint, highPoint);


  //Random walk starting point
  HyperPoint point = chosenHyperCuboid.getCenter();
  //point.at(2) = 0.3;
  //point.at(3) = -0.2;
  //point.at(4) = 2.0;

  HyperPointSet points(_numWalkers, point);

  double valCenter = getBinNumFromFunc(point);


  for ( int i = 0; i < _maximumRandWalks; i++ ){
    
    VERBOSE_LOG <<"Random walk " << i << " of " << _maximumRandWalks << std::endl;


    for (int j = 0; j < _numWalkers; j++){
      VERBOSE_LOG <<"-- walker " << j << " fancies exploring this space" << std::endl;

      walk( points.at(j), walkLimits );
      VERBOSE_LOG <<"-- walker " << j;
      double val = getBinNumFromFunc( points.at(j) );
      VERBOSE_LOG << " has fnc val = " << val << std::endl;

      if ( fabs(val - valCenter) > 0.5 ) {

        return splitByCoord(volumeNumber, dimension, points.at(j) );
        
      }

    }
    
    //make two pointx which are a short distance either size of the split

    //HyperPoint point1(point);
    //HyperPoint point2(point);
    //point1.at(dimension) += walkSize*0.5;
    //point2.at(dimension) -= walkSize*0.5;
    //
    //if ( walkLimits.inVolume(point1) == false ){
    //  point1.at(dimension) -= walkSize*0.5;
    //  point1.at(dimension) += minBinWidth*0.5;      
    //}
    //
    //if ( walkLimits.inVolume(point2) == false ){
    //  point2.at(dimension) += walkSize*0.5;
    //  point2.at(dimension) -= minBinWidth*0.5;      
    //}

    //Evaluate the function for each point
    //double val1 = _func->getVal(point1);
    //double val2 = _func->getVal(point2);

    //double val = _func->getVal(point);


    //If the function values are sufficiently different - SPLIT!!
    //if ( fabs(val1 - val2) > 0.5 ) {
    //if ( fabs(val - valCenter) > 0.5 ) {

    //  splitByCoord();

    //}

  }
  
  return 0;

}



HyperBinningMakerPhaseBinning::~HyperBinningMakerPhaseBinning(){
  
  INFO_LOG << "Number of gradient   split attempts: " << _numberOfGradientSplits   << std::endl;
  INFO_LOG << "Number of systematic split attempts: " << _numberOfSystematicSplits << std::endl;


  GOODBYE_LOG << "Goodbye from the HyperBinningMakerPhaseBinning() Constructor" <<std::endl; 
}
