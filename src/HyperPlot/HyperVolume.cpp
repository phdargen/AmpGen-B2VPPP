#include "AmpGen/HyperPlot/HyperVolume.h"

///Simple constuctor that only takes the dimensionality of 
///the HyperVolume.
///
HyperVolume::HyperVolume(int dimension) : 
  _dimension ( dimension )
{ 
}

///Constuctor that adds just one HyperCuboid to the HyperVolume,
///defined by its upper and lower corners.
HyperVolume::HyperVolume(const HyperPoint& lowCorner, const HyperPoint& highCorner) : 
  _dimension ( lowCorner.size() )
{ 
  _hyperCuboids.push_back(HyperCuboid(lowCorner, highCorner));
}


///Constuctor that adds just one HyperCuboid to the HyperVolume.
HyperVolume::HyperVolume(const HyperCuboid& cuboid) :
  _dimension ( cuboid.getDimension() )
{

  _hyperCuboids.push_back(cuboid);

}

///Constuctor that adds just two HyperCuboids to the HyperVolume.
HyperVolume::HyperVolume(const HyperCuboid& cuboid1, const HyperCuboid& cuboid2) :
  _dimension ( cuboid1.getDimension() )
{
  push_back(cuboid1);
  push_back(cuboid2);
}

///Check to see if a HyperPoint is within the HyperVolume
///
bool HyperVolume::inVolume(const HyperPoint& coords) const{

  for(unsigned int i = 0; i < _hyperCuboids.size(); i++){
    if(_hyperCuboids.at(i).inVolume(coords)==1) return 1;
  }
  return 0;

}

///Check to see if all the HyperPoint in a HyperPointSet 
///are within the HyperVolume.
///
bool HyperVolume::inVolume(const HyperPointSet& coords) const{
  
  int ncoords = coords.size();
  for (int i = 0; i < ncoords; i++){
    if ( inVolume(coords.at(i)) == false ) return false;
  }
  return true;
}

///Make a HyperVolume that contains the HyperCuboids
///from both HyperVolumes.
///
HyperVolume HyperVolume::operator+ (const HyperVolume & other) const{
  
  if (getDimension() != other.getDimension()) {
    ERROR_LOG << "You are trying to add HyperVolumes of different dimensions.";
    return *this;
  }
  HyperVolume volume(getDimension());
  for(int i = 0; i < size(); i++){
    volume.addHyperCuboid(getHyperCuboid(i));
  }
  for(int i = 0; i < other.size(); i++){
    volume.addHyperCuboid(other.getHyperCuboid(i));
  }
  return volume;
}

///Overwrite this HyperVolume with another.
///
HyperVolume & HyperVolume::operator= (const HyperVolume & other){
  
  _hyperCuboids = other._hyperCuboids;
  _dimension =    other._dimension;

  return *this;
}

///add a HyperCuboid (defined by upper and lower HyperPoints) to the HyperVolume
///
void HyperVolume::addHyperCuboid(const HyperPoint& lowCorner,const HyperPoint& highCorner){

  addHyperCuboid( HyperCuboid(lowCorner, highCorner) );
}

///add a HyperCuboid to the HyperVolume
///
void HyperVolume::addHyperCuboid(const HyperCuboid& hyperCuboid){
  if (hyperCuboid.getDimension() == _dimension) _hyperCuboids.push_back(hyperCuboid);
  else ERROR_LOG << "The HyperCuboid you are adding to this HyperVolume has the wrong dimension";
}

///add a HyperCuboid to the HyperVolume
///
void HyperVolume::push_back(const HyperCuboid& hyperCuboid){
  addHyperCuboid(hyperCuboid);
}



///Find the center of mass of the HyperVolume, assuming the mass
///is proportional to volume.
///
HyperPoint HyperVolume::getAverageCenter() const{
  
  HyperPoint binCenter(_dimension);
  double sumW = 0.0;

  for(unsigned int i = 0; i < _hyperCuboids.size(); i++){
    HyperPoint ithBinCenter = _hyperCuboids.at(i).getCenter();
    double weight =           _hyperCuboids.at(i).volume();
    
    binCenter = binCenter + ithBinCenter*weight;
    sumW += weight;
  }

  binCenter = binCenter/sumW;

  return binCenter;

}

///Find the volume of the HyperVolume (sum of the volume of the 
/// HyperCuboids, doesn't consider overlap).
///
double HyperVolume::volume() const{

  double volume = 0.0;
  for(unsigned int i = 0; i < _hyperCuboids.size(); i++){
    volume += _hyperCuboids.at(i).volume();
  }
  return volume;

}

///Find the Minimum value in a given dimension.
///
double HyperVolume::getMin(int dimension) const{

  if (dimension < _dimension) {
    double min = _hyperCuboids.at(0).getLowCorner().at(dimension);
    for(unsigned int i = 1; i < _hyperCuboids.size(); i++){
      double temp = _hyperCuboids.at(i).getLowCorner().at(dimension);
      if (min > temp) min = temp;
    }
    return min;
  }

  ERROR_LOG << "You are requesting a dimensionality that does not exist in this HyperVolume";
  return -1.0;

}

///Find the Maximum value in a given dimension.
///
double HyperVolume::getMax(int dimension) const{

  if (dimension < _dimension) {
    double max = _hyperCuboids.at(0).getHighCorner().at(dimension);
    for(unsigned int i = 1; i < _hyperCuboids.size(); i++){
      double temp = _hyperCuboids.at(i).getHighCorner().at(dimension);
      if (max < temp) max = temp;
    }
    return max;
  }

  ERROR_LOG << "You are requesting a dimensionality that does not exist in this HyperVolume";
  return -1.0;

}

///get the limits of the HyperVolume
///
HyperCuboid HyperVolume::getLimits() const{
  HyperCuboid limits(getDimension());
  for (int i = 0; i < getDimension(); i++){
    limits.getLowCorner ().at(i) = getMin(i);
    limits.getHighCorner().at(i) = getMax(i);
  }
  return limits;
}


///Print the HyperCuboids that define the HyperVolume
///
void HyperVolume::print(std::ostream& os, int endline) const{

  for(unsigned int i = 0; i < _hyperCuboids.size(); i++){
    _hyperCuboids.at(i).print(os, 1);
  }  
  if (endline == 1) os << std::endl;
}

///Used for taking a slice through a HyperVolume. Pick a HyperPoint 'coords'
///e.g coords =(-1,5,4) then pick the dimensions you wish to slice through e.g. dims = (0).
///You can then imagine the plane defined by x_0 = -1, which intersects your HyperVolume.
///The 2D shape of this intersection is what will be returned.
///Further if you were to set dims = (0,1) you would get the 1D bin boundarys defined by the intersection
///of the line (x_0 = -1, x_1 = 5) and the HyperVolume.
HyperVolume HyperVolume::slice(const HyperPoint& coords, std::vector<int> dims) const{
  
  int currentDim = getDimension();
  int sliceDim   = (int)dims.size();  
  int newDim     = currentDim - sliceDim;

  HyperVolume slicedVolume(newDim);

  for(unsigned int i = 0; i < _hyperCuboids.size(); i++){
    if(_hyperCuboids.at(i).inVolume(coords, dims)==1){
      //std::cout << "Woohoo - I'm in volume, lets project" << std::endl;
      HyperCuboid cube = _hyperCuboids.at(i).projectOpposite(dims);
      //std::cout << " ---- wow, it worked" << std::endl;
      slicedVolume.addHyperCuboid(cube);
    } 
  }

  return slicedVolume;

}

///Split all of the HyperCuboids within the HyperVolume in the
///dimension given, and at the fractionalSplitPoint given.
HyperVolume HyperVolume::splitAll(int dimension, double fractionalSplitPoint){

  HyperVolume hyperCuboidSet( getDimension() );

  for (int i = 0; i < size(); i++){
    at(i).split(dimension, fractionalSplitPoint, hyperCuboidSet);
  }

  return hyperCuboidSet;

}



///Destructor
///
HyperVolume::~HyperVolume() { 
}




