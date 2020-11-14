#include "AmpGen/HyperPlot/HyperCuboid.h"

///Most basic constructor where only the dimension
///of the cuboid is specified.
HyperCuboid::HyperCuboid(int dimension) : 
  _dimension ( dimension ),
  _lowCorner ( dimension ),
  _highCorner( dimension )
{ 
}

///Construct the HyperCuboid with two HyperPoints, one in the low corner, 
///and one in the high corner
HyperCuboid::HyperCuboid(const HyperPoint& lowCorner, const HyperPoint& highCorner) :
  _dimension( lowCorner.size() ),
  _lowCorner ( _dimension ),  //note this assignment sets the dimension of this variable forever, so if 
  _highCorner( _dimension )   //low and high corners are a different size, an error will be thrown inside HyperPoint
{

  if (lowCorner.allLT(highCorner)){
    _lowCorner  = lowCorner;
    _highCorner = highCorner;
  }
  else{
    ERROR_LOG << "Your lowCorner isn't lower than your highCorner"<< std::endl;
    _lowCorner .print();
    _highCorner.print();
  }

}

///Construct the HyperCuboid with two HyperPoints (x, x, x ....) 
///and (y, y, y ....)
HyperCuboid::HyperCuboid(int dimension, double low, double high) :
  _dimension ( dimension ),
  _lowCorner ( dimension ),  //note this assignment sets the dimension of this variable forever, so if 
  _highCorner( dimension )   //low and high corners are a different size, an error will be thrown inside HyperPoint
{
  
  for (int i = 0; i < _dimension; i++){
    _lowCorner.at(i)  = low ;
    _highCorner.at(i) = high;
  }

  if (low > high) {
    ERROR_LOG << "Your lowCorner isn't lower than your highCorner" << std::endl;
  }
}

///See if two HyperCuboids are the same
///
bool HyperCuboid::operator ==(const HyperCuboid& other) const{
  if (other.getDimension() != this->getDimension()) {
    ERROR_LOG << "Trying to compare HyperCuboids of different dimensions, returning false";
    return false; 
  }
  if (getLowCorner()  != other.getLowCorner() ) return false;
  if (getHighCorner() != other.getHighCorner()) return false;
  return true;
}

///See if two HyperCuboids are not the same
///
bool HyperCuboid::operator !=(const HyperCuboid& other) const{
  if (other.getDimension() != this->getDimension()) {
    ERROR_LOG << "Trying to compare HyperCuboids of different dimensions, returning false";
    return false; 
  }
  if (getLowCorner()  != other.getLowCorner() ) return true;
  if (getHighCorner() != other.getHighCorner()) return true;
  return false;
}

///Inflate the boundaries of the HyperCuboid by some percentage.
///
/// x_low_0' = x_low_0 - (x_high_0 - x_low_0) * percent
///
/// x_high_0' = x_high_0 + (x_high_0 - x_low_0) * percent
const HyperCuboid& HyperCuboid::inflateCuboid(double percent){

  for (int i = 0; i < getDimension(); i++){
    double diff = getHighCorner().at(i) - getLowCorner().at(i);
    double increase = diff*percent;
    getHighCorner().at(i) = getHighCorner().at(i) + increase;
    getLowCorner() .at(i) = getLowCorner() .at(i) - increase;
  }  
  return *this;
}

///Get the HyperPoint at the center of the HyperCuboid
///
HyperPoint HyperCuboid::getCenter() const{
  return (_lowCorner + _highCorner)*0.5;
}

///Print the HyperCuboid to given std::ostream
///
void HyperCuboid::print(std::ostream& os, int endline) const{

  os << "Lower Corner: ";
  _lowCorner.print(os, 0);
  os << "     High Corner: ";
  _highCorner.print(os, 0);
  if (endline) os << std::endl;

}

///Project the HyperCuboid into a lower dimensional space.
///E.g. start with HyperCuboid defined by (0,-1,-2) , (0,+1,+2)
///and project over dimensions 0 and 2. This would return a
///HyperCuboid defined by (0,-2) , (0,+2)
HyperCuboid HyperCuboid::project(std::vector<int> dims) const{
  
  int newdim = (int)dims.size();
  HyperPoint lowCorner (newdim);
  HyperPoint highCorner(newdim);

  for (int i = 0; i < newdim; i++){
    lowCorner .at(i) = getLowCorner ().at( dims.at(i) );
    highCorner.at(i) = getHighCorner().at( dims.at(i) );
  }

  return HyperCuboid(lowCorner, highCorner);

}

///Project the HyperCuboid into a lower dimensional space.
///E.g. start with HyperCuboid defined by (0,-1,-2) , (0,+1,+2)
///and projectOpposite over dimensions 1. This would return a
///HyperCuboid defined by (0,-2) , (0,+2)
HyperCuboid HyperCuboid::projectOpposite(std::vector<int> dims) const{
  
  int currentDim = getDimension();
  int newdim     = (int)dims.size();
  
  std::vector<int> projdims;

  for (int i = 0; i < currentDim; i++){
    bool doesExist = false;
    for (int j = 0; j < newdim; j++){
      int dim = dims.at(j);
      if (i == dim) doesExist = true;
    }

    if (doesExist == false) projdims.push_back(i);
  }  

  return project(projdims);

}

///Get point on the opposite side of the HyperCuboid 
///
///Imagine a vector V going from the center of the
///HyperCuboid, C, to the given point. This returns
/// C - V
HyperPoint HyperCuboid::getOppositePoint(const HyperPoint& point) const{

  HyperPoint shift = point - _lowCorner;
  HyperPoint opposite = _highCorner - shift;
  
  return opposite;
}

///Get a random point in the HyperCuboid (uses gRandom)
///
HyperPoint HyperCuboid::getRandomPoint(TRandom* random) const{
  
  HyperPoint point(getDimension());
  for (int i = 0; i < getDimension(); i++){
    point.at(i) = random->Uniform(getLowCorner().at(i), getHighCorner().at(i));
  }
  return point;

}

///Get a random point in the HyperCuboid (uses gRandom)
///
HyperPointSet HyperCuboid::getRandomPoints(int nPoints, TRandom* random) const{
  
  HyperPointSet points(getDimension());
  for (int i = 0; i < nPoints; i++){
    points.push_back( getRandomPoint(random) );
  }
  return points;
  
}



///Give it a corner of the HyperCuboid and it will return all faces (HyperPlanes) connected 
///to that corner. Note that nDim points are needed to define such a HyperPlane. 
///
/// e.g. for a square the 'faces' are lines  defined by 2 points 
///      for a cube   the 'faces' are planes defined by 3 points 
std::vector<HyperPlane> HyperCuboid::getConnectedHyperPlanes(const HyperPoint& point) const{
  
  std::vector<HyperPlane> connectedHyperPlanes;
  
  //First find all the verticies that are connected to this one by an edge.
  //There will be nDim such verticies
  int dim = getDimension();
  HyperPointSet connectedVerticies = getConnectedVerticies(point);
  
  //We now have the vertex given, and an additional nDim connected 
  //verticies. There is therefore nDim possible Hyperplanes one can form.
  for (int i = 0; i < dim; i++){
    HyperPointSet pointsInPlane(dim);
    pointsInPlane.push_back(point);
    for(int j = 0; j < dim; j++){
      if (i != j) pointsInPlane.push_back(connectedVerticies.at(j));
    }
    pointsInPlane.sort();
    HyperPlane plane(pointsInPlane);
    connectedHyperPlanes.push_back(plane);
  }
  
  return connectedHyperPlanes;

}


///Give it a corners of the HyperCuboid and it will return all faces (HyperPlanes) connected 
///to that corner. Note that nDim points are needed to define such a HyperPlane. 
///
/// e.g. for a square the 'faces' are lines  defined by 2 points 
///      for a cube   the 'faces' are planes defined by 3 points 
///
///This is done for all corners given in the HyperPointSet
std::vector<HyperPlane> HyperCuboid::getConnectedHyperPlanes(const HyperPointSet& point) const{
   
  std::vector<HyperPlane> connectedHyperPlanes;

  for(unsigned i = 0; i < point.size(); i++){
    std::vector<HyperPlane> someHyperPlanes = getConnectedHyperPlanes(point.at(i));
    for(unsigned j = 0; j < someHyperPlanes.size(); j++){
      connectedHyperPlanes.push_back(someHyperPlanes.at(j));
    }
  }
  
  return connectedHyperPlanes;

}

///Get all the HyperPlanes that define the `faces' of the HyperCuboid.
///Note that nDim points are needed to define such HyperPlanes. 
///
/// e.g. for a square the 'faces' are lines  defined by 2 points 
///      for a cube   the 'faces' are planes defined by 3 points 
///
std::vector<HyperPlane> HyperCuboid::getBoundaryHyperPlanes() const{
  
  //Get all the verticies of the HyperCuboid, then find all the HyperPlanes
  //that are conected to these verticies.
  std::vector<HyperPlane> planes = getConnectedHyperPlanes(getVertices());
  std::vector<int> isUnique(planes.size(), 1);
  

  //remove all the duplicates.

  VERBOSE_LOG << "There are " << planes.size() << " planes before removing doubles"; 
  

  for(unsigned j = 0; j < planes.size(); j++){

    if (isUnique.at(j) == 0){ VERBOSE_LOG << "    Skippied"; continue; }

    for(unsigned i = j+1; i < planes.size(); i++){

      if (isUnique.at(i) == 0) continue;
      if (planes.at(i) == planes.at(j)) isUnique.at(i) = 0;

    }
  }

  std::vector<HyperPlane> faces;
  for(unsigned j = 0; j < planes.size(); j++){
    if (isUnique.at(j) == 1) faces.push_back(planes.at(j));
  }


  VERBOSE_LOG << "There are " << faces.size() << " planes after removing doubles"; 


  return faces;

}

///It takes a while to find all the HyperPlanes (that define the 
/// faces of the HyperCuboid), so they are cashed using this
/// function. Beware if you later change the corners of the HyperCuboid
void HyperCuboid::updateFaceCash() const{
  
  if(_faces.size() == 0) _faces = getBoundaryHyperPlanes();

}

///A HyperLine is parameterised by (v_1 + x.v_2). This finds one of the two points, x_+,
///where the HyperLine intersects a face of the HyperCuboid
double HyperCuboid::getPositiveIntersectionParameter(const HyperLine& line) const{

  updateFaceCash();
  
  int foundSolution = 0;
  double negativeLineParam = 0.0;

  for(unsigned i = 0; i < _faces.size(); i++){
    HyperPoint LineParam  = _faces.at(i).findLineIntersectionParameter(line);
    HyperPoint PlaneParam = _faces.at(i).findPlaneIntersectionParameters(line);

    VERBOSE_LOG << "-----------------------" << i << "--------------------------";
    VERBOSE_LOG << "Get the parameters for the line and the plane";  

    if (LineParam.getDimension()==0){
      VERBOSE_LOG << "Matrix is singular...";
      continue;
    }

    if (PlaneParam <= 1.0 && PlaneParam >= 0.0) {
      VERBOSE_LOG << "This intersection is within the cuboid";
      if(LineParam > 0 && foundSolution == 1){
        double otherPoint = LineParam.at(0);
        double diff = otherPoint - negativeLineParam;
        if (fabs(diff)>10e-2) ERROR_LOG << "Somehow there are two distinct negative intercepts";
        else ERROR_LOG << "Already have this solution... weird";
      }
      else if(LineParam > 0) {
        VERBOSE_LOG << "This positive intersection point " << LineParam.at(0);
        
        negativeLineParam = LineParam.at(0);
        foundSolution = 1;
      }
      else{
        VERBOSE_LOG << "This positive intersection point" << LineParam.at(0);
      }
    }
    else{
      VERBOSE_LOG << "This intersection is not within the cuboid";
    }

  }

  if (foundSolution == 0) ERROR_LOG << "Found no intercept in cuboid.";

  return negativeLineParam;

}

///A HyperLine is parameterised by (v_1 + x.v_2). This finds one of the two points, x_-,
///where the HyperLine intersects a face of the HyperCuboid.
double HyperCuboid::getNegativeIntersectionParameter(const HyperLine& line) const{

  updateFaceCash();
  
  int foundSolution = 0;
  double negativeLineParam = 0.0;

  for(unsigned i = 0; i < _faces.size(); i++){
    HyperPoint LineParam  = _faces.at(i).findLineIntersectionParameter(line);
    HyperPoint PlaneParam = _faces.at(i).findPlaneIntersectionParameters(line);

    VERBOSE_LOG << "-----------------------" << i << "--------------------------";
    VERBOSE_LOG << "Get the parameters for the line and the plane";  

    if (LineParam.getDimension()==0){
      VERBOSE_LOG << "Matrix is singular...";
      continue;
    }

    if (PlaneParam <= 1.0 && PlaneParam >= 0.0) {
      VERBOSE_LOG << "This intersection is within the cuboid";
      if(LineParam < 0 && foundSolution == 1){
        double otherPoint = LineParam.at(0);
        double diff = otherPoint - negativeLineParam;
        if (fabs(diff)>10e-2) ERROR_LOG << "Somehow there are two distinct negative intercepts";
        else ERROR_LOG << "Already have this solution... weird";
      }
      else if(LineParam < 0) {
        VERBOSE_LOG << "This positive intersection point " << LineParam.at(0);

        negativeLineParam = LineParam.at(0);
        foundSolution = 1;
      }
      else{
        VERBOSE_LOG << "This positive intersection point" << LineParam.at(0);
      }
    }
    else{
      VERBOSE_LOG << "This intersection is not within the cuboid";
    }

  }

  if (foundSolution == 0) ERROR_LOG << "Found no intercept in cuboid.";

  return negativeLineParam;
} 

///A HyperLine is parameterised by (v_1 + x.v_2). This finds one of the two points, (v_1 + x_+.v_2),
///where the HyperLine intersects a face of the HyperCuboid
HyperPoint HyperCuboid::getPositiveIntersectionPoint(const HyperLine& line) const{

  return line.getParametricPoint(getPositiveIntersectionParameter(line));

}

///A HyperLine is parameterised by (v_1 + x.v_2). This finds one of the two points, (v_1 + x_-.v_2),
///where the HyperLine intersects a face of the HyperCuboid
HyperPoint HyperCuboid::getNegativeIntersectionPoint(const HyperLine& line) const{

  return line.getParametricPoint(getNegativeIntersectionParameter(line));

}

HyperPointSet HyperCuboid::getEdgeCenters() const{

  int ndim = getDimension();

  std::bitset< 15 > binary( 0 );

  HyperPoint center = getCenter();
  HyperPoint vectorToHighCorner = _highCorner - center;

  int nVertices = pow(2, ndim);
  
  HyperPointSet verticies(ndim);
  
  for (int dimToIgnore = 0; dimToIgnore < ndim; dimToIgnore++){

    for (int i = 0; i < nVertices; i++){
      HyperPoint vectorToVertex(vectorToHighCorner);
      binary = i;
      for (int dim = 0; dim < ndim; dim++){
        if (binary[dim] == 0){
          vectorToVertex.at(dim) = -vectorToVertex.at(dim);
        }
        if (dimToIgnore == dim){
          vectorToVertex.at(dim) = 0.0;
        }        
      }
      verticies.push_back(vectorToVertex + center);
    }  
  }
  
  verticies.sort();
  verticies.removeDuplicates();

  return verticies;

}




///Get all verticies of the HyperCuboid with no repeats.
///
///This is done by starting from one corner and following
///the edges to the other corners. This has to be done nDim
///times in order to get all the corners.
HyperPointSet HyperCuboid::getVertices() const{
  
  //New faster method. We know the number of verticies, so just
  //loop from 0 -> nVertex-1 and convert this number to binary.
  //Imagine the binary number is a vector V. 
  //Also have the center of the cube, C, and a vector to the top right
  //corner T.
  //
  //

  int ndim = getDimension();

  std::bitset< 15 > binary( 0 );

  HyperPoint center = getCenter();
  HyperPoint vectorToHighCorner = _highCorner - center;

  int nVertices = pow(2, ndim);
  
  HyperPointSet verticies(ndim);

  for (int i = 0; i < nVertices; i++){
    HyperPoint vectorToVertex(vectorToHighCorner);
    binary = i;
    for (int dim = 0; dim < ndim; dim++){
      if (binary[dim] == 0){
        vectorToVertex.at(dim) = -vectorToVertex.at(dim);
      }
    }
    verticies.push_back(vectorToVertex + center);
  }
  
  return verticies;

  //OLD method
  //HyperPointSet hyperPointSet(getDimension());
  //hyperPointSet.push_back(_lowCorner);
  //
  //HyperPointSet previousHyperPointSet = hyperPointSet;

  //for(int i = 1; i < getDimension(); i++){
  //  previousHyperPointSet = getConnectedVerticies(previousHyperPointSet);
  //  hyperPointSet.addHyperPointSet(previousHyperPointSet);
  //  hyperPointSet.sort();
  //  hyperPointSet.removeDuplicates();
  //}
  //
  //hyperPointSet.push_back(_highCorner);


  //return hyperPointSet;

}

///If a vertex of the HyperCuboid is given, this finds all
///other verticies that are connected to that one by an edge.
///It does this for all HyperPoints in the HyperPointSet given.
///
///For example, if you give it the top right corner of a square,
/// it will return the top left and bottom right corners.
HyperPointSet HyperCuboid::getConnectedVerticies(const HyperPointSet& pointSet) const{
  
  HyperPointSet hyperPointSet(getDimension());
  
  for (unsigned int i = 0; i < pointSet.size(); i++){
    hyperPointSet.addHyperPointSet ( getConnectedVerticies(pointSet.at(i)) );
    VERBOSE_LOG << "Sorting....";
    hyperPointSet.sort();
    VERBOSE_LOG << "Removing Duplicates....";
    hyperPointSet.removeDuplicates();
  }
  
  return hyperPointSet;
}

///If a vertex of the HyperCuboid is given, this finds all
///other verticies that are connected to that one by an edge
///
///For example, if you give it the top right corner of a square,
/// it will return the top left and bottom right corners.
HyperPointSet HyperCuboid::getConnectedVerticies(const HyperPoint& point) const{
  
  int dim = getDimension();

  HyperPointSet hyperPointSet(dim);

  HyperPoint oppositePoint    = getOppositePoint(point);
  HyperPoint vectorToOpposite = oppositePoint - point;

  for (int i = 0; i < dim; i++){
    HyperPoint translation(dim);
    translation.at(i) = vectorToOpposite.at(i);
    HyperPoint newCorner = point + translation;
    hyperPointSet.push_back(newCorner);
  }
  
  return hyperPointSet;

}

///Set the upper and lower corners of the HyperCuboid
///
bool HyperCuboid::setCorners(const HyperPoint& lowCorner, const HyperPoint& highCorner){

  if (lowCorner.allLTOE(highCorner)){
    _lowCorner  = lowCorner;
    _highCorner = highCorner;
  }
  else{
    ERROR_LOG << "Your lowCorner isn't lower than your highCorner" <<std::endl;
    return 0;
  }
  return 1;
}

///See if a HyperPoint is within the HyperCuboid volume
///
bool HyperCuboid::inVolume(const HyperPoint& coords) const{

  if (_lowCorner.allLT(coords) &&  _highCorner.allGTOE(coords)) return 1;
  return 0;

}

///See if a HyperPoint is within the HyperCuboid volume, but only
///for selected dimensions.
bool HyperCuboid::inVolume(const HyperPoint& coords, std::vector<int> dims) const{

  for (unsigned i = 0; i < dims.size(); i++){
    int dim = dims.at(i);
    double minEdge = getLowCorner ().at(dim);
    double maxEdge = getHighCorner().at(dim);
    double val     = coords.at(dim);
    if ( ( minEdge < val && val <= maxEdge ) == false ) return false;
  }

  return true;

}

///Find the volume of the HyperCuboid
///
double HyperCuboid::volume() const{
  
  return (_highCorner - _lowCorner).multiplyElements();

}

///Find width of the HyperCuboid as a HyperPoint
///
HyperPoint HyperCuboid::getWidth() const{
  return _highCorner - _lowCorner;  
}

///Find the width of the HyperCuboid in a given dimension
///
double HyperCuboid::getWidth(int dim) const{
  return _highCorner.at(dim) - _lowCorner.at(dim);
}


///split the cuboid into two along the dimension given and 
///return the resulting HyperCuboid. The fractional split 
///point [0,1] decides where (in that dimesnion) the split is
HyperCuboid HyperCuboid::splitBelow(int dimension, double fractionalSplitPoint) const{
  
  if (dimension < 0 || dimension >= getDimension()){
    ERROR_LOG << "You cannot split in a dimension that doesn't exist!!" << std::endl;
    return HyperCuboid(*this);
  }
  if (fractionalSplitPoint <= 0.0 || fractionalSplitPoint >= 1.0){
    ERROR_LOG << "Your fractional split point must be between 0.0 and 1.0" << std::endl;
    return HyperCuboid(*this);
  }

  HyperPoint lowCorner1 (_lowCorner );
  HyperPoint highCorner1(_highCorner);
  
  double low   = _lowCorner .at(dimension);
  double high  = _highCorner.at(dimension);
  double width = high - low;

  double splitCoord = low + width*fractionalSplitPoint;

  highCorner1.at(dimension) = splitCoord;

  HyperCuboid lowCuboid (lowCorner1, highCorner1); 
  
  return lowCuboid;

} 

///split the cuboid into two along the dimension given and 
///return the resulting HyperCuboid. The fractional split 
///point [0,1] decides where (in that dimesnion) the split is
HyperCuboid HyperCuboid::splitAbove(int dimension, double fractionalSplitPoint) const{
  
  if (dimension < 0 || dimension >= getDimension()){
    ERROR_LOG << "You cannot split in a dimension that doesn't exist!!" << std::endl;
    return HyperCuboid(*this);
  }
  if (fractionalSplitPoint <= 0.0 || fractionalSplitPoint >= 1.0){
    ERROR_LOG << "Your fractional split point must be between 0.0 and 1.0" << std::endl;
    return HyperCuboid(*this);
  }

  HyperPoint lowCorner2 (_lowCorner );
  HyperPoint highCorner2(_highCorner);
  
  double low   = _lowCorner .at(dimension);
  double high  = _highCorner.at(dimension);
  double width = high - low;

  double splitCoord = low + width*fractionalSplitPoint;

  lowCorner2 .at(dimension) = splitCoord;

  HyperCuboid highCuboid(lowCorner2, highCorner2); 
  
  return highCuboid;

}

///split the cuboid into two along the dimension given and 
///return the resulting HyperCuboids in a HyperCuboidSet. The fractional split 
///point [0,1] decides where (in that dimesnion) the split is
HyperVolume HyperCuboid::split(int dimension, double fractionalSplitPoint) const{

  return HyperVolume(
             splitBelow(dimension,fractionalSplitPoint), 
             splitAbove(dimension,fractionalSplitPoint)
             );

}

///split the cuboid into two along the dimension given and 
///append the resulting HyperCuboids to the HyperCuboidSet given. The fractional split 
///point [0,1] decides where (in that dimesnion) the split is
void HyperCuboid::split(int dimension, double fractionalSplitPoint, HyperVolume& set) const{

  set.push_back( splitBelow(dimension,fractionalSplitPoint) );
  set.push_back( splitAbove(dimension,fractionalSplitPoint) );

}

///Destructor
///
HyperCuboid::~HyperCuboid() { 
}


