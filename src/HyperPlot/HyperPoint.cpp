#include "AmpGen/HyperPlot/HyperPoint.h"

///Most basic constructor that makes a HyperPoint with a 
///specified dimension, and filled with 0.0 for each element.
HyperPoint::HyperPoint(int dimension) :
  _coords(dimension, 0.0)
{
}

///Constructor that makes a HyperPoint with a 
///specified dimension, and fills each element with
///a specified value.
HyperPoint::HyperPoint(int dimension, double val) :
  _coords(dimension, val)
{
}

///Constuctor that takes a std::vector of elements
///
HyperPoint::HyperPoint(std::vector<double> coords) :
  _coords(coords)
{
}

///Constuctor for a 1D space.
///
HyperPoint::HyperPoint(double x1) :
  _coords(1, 0.0)
{
  at(0) = x1;
}

///Constuctor for a 2D space.
///
HyperPoint::HyperPoint(double x1, double x2) :
  _coords(2, 0.0)
{
  at(0) = x1;
  at(1) = x2;
}

///Constuctor for a 3D space.
///
HyperPoint::HyperPoint(double x1, double x2, double x3) :
  _coords(3, 0.0)
{
  at(0) = x1;
  at(1) = x2;
  at(2) = x3;
}

///Constuctor for a 4D space.
///
HyperPoint::HyperPoint(double x1, double x2, double x3, double x4) :
  _coords(4, 0.0)
{
  at(0) = x1;
  at(1) = x2;
  at(2) = x3;
  at(3) = x4;
}

///Constuctor for a 5D space.
///
HyperPoint::HyperPoint(double x1, double x2, double x3, double x4, double x5) :
  _coords(5, 0.0)
{
  at(0) = x1;
  at(1) = x2;
  at(2) = x3;
  at(3) = x4;
  at(4) = x5;
}

HyperPoint::HyperPoint(double x1, double x2, double x3, double x4, double x5, double x6) :
_coords(6, 0.0)
{
    at(0) = x1;
    at(1) = x2;
    at(2) = x3;
    at(3) = x4;
    at(4) = x5;
    at(5) = x6;
}

HyperPoint::HyperPoint(double x1, double x2, double x3, double x4, double x5, double x6, double x7) :
_coords(7, 0.0)
{
    at(0) = x1;
    at(1) = x2;
    at(2) = x3;
    at(3) = x4;
    at(4) = x5;
    at(5) = x6;
    at(6) = x7;
}


///Multiply all the elements together.
///
double HyperPoint::multiplyElements() const{
  double total = 1.0;
  for (int i = 0; i < size(); i++) total *= this->at(i);
  return total;
}

///Check if two HyperPoints are compatible  (of the same dimension)
///
bool HyperPoint::compatible(const HyperPoint & other, bool printError) const{

  if (getDimension() != other.getDimension()) {
    if (printError == true) ERROR_LOG << "HyperPoints are of different dimension i.e. not compatible!!!"; 
    return 0;
  }
  return 1;

}

///Find the distance between this point and another
///
double HyperPoint::distanceTo(const HyperPoint & other) const{
  if (this->compatible(other) == 0) return 0.0; 
  HyperPoint temp = other - *this;
  return temp.norm();
}

///Find the norm = sqrt(x1^2 + x2^2 + ...)
///
double HyperPoint::norm() const{
  double normVal = 0.0;
  for (int i = 0; i < size(); i++) normVal += at(i)*at(i);
  return sqrt(normVal);
}

///Find the norm squared =(x1^2 + x2^2 + ...)
///
double HyperPoint::norm2() const{
  double normVal = 0.0;
  for (int i = 0; i < size(); i++) normVal += at(i)*at(i);
  return normVal;
}

///Find the dot product of two points
///
double     HyperPoint::dotProduct(const HyperPoint & other) const{
  if (this->compatible(other) == 0) return 0.0; 

  double val = 0.0;
  for (int i = 0; i < size(); i++) val += at(i)*other.at(i);
  return val;
}

///Fill the point with random numbers, uniformly distributed 
///between min and max. Uses gRandom for random number.
void       HyperPoint::fillRandom(double min, double max){

  for(int i = 0; i < size(); i++) at(i) = gRandom->Uniform(min,max);
}

///Project this vector onto another = v1 (v1.v2/v1.v1)
///
HyperPoint HyperPoint::project(const HyperPoint & other) const{
  if (this->compatible(other) == 0) return *this; 

  HyperPoint projection(*this);
  double scale = this->dotProduct(other) / this->dotProduct(*this);
  return (projection * scale);
}

///Copy all elements and weights from other HyperPoint to this.
///Only function that allows the dimensionality of a point to
///be changed after it is created. Should be careful with this,
///for instance a HyperPointSet only accepts HyperPoint's of the
///same dimension. If you change the dimensionality of a point
///after it has been added, this could cause big problems.
HyperPoint& HyperPoint::operator= (const HyperPoint & other){
  //if (this->compatible(other) == 0) return *this; 
  
  _coords  = other._coords;
  _weights = other._weights;
  
  return *this;
}

///Add elements of two HyperPoints together - leave weights unchanged.
///
HyperPoint HyperPoint::operator+ (const HyperPoint & other) const{

  if (this->compatible(other) == 0) return *this; 

  HyperPoint temp(size());

  for (int i = 0; i < size(); i++) temp.at(i) = this->at(i) + other.at(i);
  
  return temp;

}

///Subtract one HyperPoint from another - leave weights unchanged.
///
HyperPoint HyperPoint::operator- (const HyperPoint & other) const{

  if (this->compatible(other) == 0) return *this; 

  HyperPoint temp(size());

  for (int i = 0; i < size(); i++) temp.at(i) = this->at(i) - other.at(i);
  
  return temp;

}

///Subtract one value from all elements of the HyperPoint.
///Leave weights unchanged.
HyperPoint HyperPoint::operator- (const double & other) const{

  HyperPoint temp(size());

  for (int i = 0; i < size(); i++) temp.at(i) = this->at(i) - other;

  return temp;

}

///Add one value to all elements of the HyperPoint.
///Leave weights unchanged.
HyperPoint HyperPoint::operator+ (const double & other) const{

  HyperPoint temp(size());

  for (int i = 0; i < size(); i++) temp.at(i) = this->at(i) + other;

  return temp;

}

///Divide all elements of a HyperPoint by one value.
///Leave weights unchanged.
HyperPoint HyperPoint::operator/ (const double & other) const{

  HyperPoint temp(size());

  for (int i = 0; i < size(); i++) temp.at(i) = this->at(i) / other;

  return temp;

}

///Multiply all elements of a HyperPoint by one value.
///Leave weights unchanged.
HyperPoint HyperPoint::operator* (const double & other) const{

  HyperPoint temp(size());

  for (int i = 0; i < size(); i++) temp.at(i) = this->at(i) * other;

  return temp;

}

///Perform a linear transformation on the HyperPoint using a 
///specified matrix. v_new = v_old M  
///Leave weights unchanged.
HyperPoint HyperPoint::linearTransformation( const TMatrixD& matrix ){

  HyperPoint point(getDimension());

  for (int i = 0; i < getDimension(); i++){
    double val = 0.0;
    for (int j = 0; j < getDimension(); j++){
      val += at(j)*matrix(i,j);
    }
    point.at(i) = val;
  }

  return point;

} 

///Perform a linear transformation on the HyperPoint using a 
///specified matrix. v_new = v_old M  
///Leave weights unchanged.
bool HyperPoint::operator <(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) < other.at(i)) return true;
    if (this->at(i) > other.at(i)) return false;
  }
  return false;
}

///Compare two HyperPoints. First compares 0th element
///of each HyperPoint, if these are the same, compares the 
///next, etc.
bool HyperPoint::operator >(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) < other.at(i)) return false;
    if (this->at(i) > other.at(i)) return true;
  }
  return false;
}

///Compare two HyperPoints. First compares 0th element
///of each HyperPoint, if these are the same, compares the 
///next, etc.
bool HyperPoint::operator <=(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) < other.at(i)) return true;
    if (this->at(i) > other.at(i)) return false;
  }
  return true;
}

///Compare two HyperPoints. First compares 0th element
///of each HyperPoint, if these are the same, compares the 
///next, etc.
bool HyperPoint::operator >=(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) < other.at(i)) return false;
    if (this->at(i) > other.at(i)) return true;
  }
  return true;
}

///Checks if every element of v1 is less than the  
///same element in v2. i.e. v1_0 < v2_0 && v1_1 < v2_1 && ...
bool HyperPoint::allLT (const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) >= other.at(i)) return false;
  }
  return true;
}

///Checks if every element of v2 is greater than the  
///same element in v2. i.e. v1_0 > v2_0 && v1_1 > v2_1 && ...
bool HyperPoint::allGT  (const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) <= other.at(i)) return false;
  }
  return true;
}

///Checks if every element of v1 is less than or equal to the  
///same element in v2. i.e. v1_0 <= v2_0 && v1_1 <= v2_1 && ...
bool HyperPoint::allLTOE(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 

  for (int i = 0; i < size(); i++) {
    if (this->at(i) > other.at(i)) return false;
  }
  return true;
}

///Checks if every element of v1 is greater than or equal to the  
///same element in v2. i.e. v1_0 >= v2_0 && v1_1 >= v2_1 && ...
bool HyperPoint::allGTOE(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 
  
  for (int i = 0; i < size(); i++) {
    if (this->at(i) < other.at(i)) return false;
  }
  return true;
}

///Checks if two HyperPoints are identical
///
bool HyperPoint::operator ==(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 
  
  for (int i = 0; i < size(); i++) if (this->at(i) != other.at(i)) return false;
  return true;
}

///Checks if two HyperPoints are not identical
///
bool HyperPoint::operator !=(const HyperPoint& other) const{
  if (this->compatible(other) == 0) return false; 
  
  for (int i = 0; i < size(); i++) if (this->at(i) != other.at(i)) return true;
  return false;
}

///Checks if all elements of a HyperPoint are smaller than a given value
///
bool HyperPoint::operator <(const double& other) const{
  for (int i = 0; i < size(); i++) if ((this->at(i) < other)==0) return false;
  return true;
}

///Checks if all elements of a HyperPoint are greater than a given value
///
bool HyperPoint::operator >(const double& other) const{
  for (int i = 0; i < size(); i++) if ((this->at(i) > other)==0) return false;
  return true;
}

///Checks if all elements of a HyperPoint are smaller than or equal to a given value
///
bool HyperPoint::operator <=(const double& other) const{
  for (int i = 0; i < size(); i++) if ((this->at(i) <= other)==0) return false;
  return true;
}

///Checks if all elements of a HyperPoint are greater than or equal to a given value
///
bool HyperPoint::operator >=(const double& other) const{
  for (int i = 0; i < size(); i++) if ((this->at(i) >= other)==0) return false;
  return true;
}

///Checks if all elements of a HyperPoint are equal to a given value
///
bool HyperPoint::operator ==(const double& other) const{
  for (int i = 0; i < size(); i++) if ((this->at(i) == other)==0) return false;
  return true;
}

///returns a const reference to a choosen element of the HyperPoint.
///Can be used for getting the value, but not setting
const double& HyperPoint::at  (int i) const{
  if (i >= this->size()) {
    ERROR_LOG << "The component of this HyperPoint you have selected is out of range";
    return _coords.at(this->size() - 1);
  }
  return _coords.at(i);
}

///returns a reference to a choosen element of the HyperPoint.
///Can be used for getting or setting it's value.
double& HyperPoint::at  (int i){
  if (i >= size()) {
    ERROR_LOG << "The component of this HyperPoint you have selected is out of range";
    return _coords.at(size() - 1);
  }
  return _coords.at(i);
}

///Print the HyperPoint to the given std::ostream (default is std::cout).
///
void HyperPoint::print(std::ostream& os, int endline) const{
  
  std::ostringstream oss; 

  oss << "( ";
  for(int i = 0; i < (size() - 1); i++){
    oss << std::setw(10) << std::left << _coords.at(i) << ", ";
  }
  if (size() != 0) oss << std::setw(10) << std::left << _coords.at(size() - 1);
  oss << ")";    
  
  this->printWeight(oss, 0);

  if (endline) oss << std::endl;

  os << oss.str();
}

std::ostream& operator<<(std::ostream& os, const HyperPoint& point) {
  point.print(os, 0);
  return os;
} 

///Destructor
///
HyperPoint::~HyperPoint(){
}


