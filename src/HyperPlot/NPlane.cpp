#include "AmpGen/HyperPlot/NPlane.h"

///Standard constuctor that takes between 1 and nDim HyperPoints.
///Internally this is converted into an origin and between 0 and 
/// (nDim - 1) vectors.
NPlane::NPlane(const HyperPointSet& a) :
  _origin ( a.at(0)          ),
  _v      ( a.getDimension() )
{
  
  WELCOME_LOG << "Hello from the NPlane() Constructor";
  
  for(unsigned i = 1; i < a.size(); i++){
    _v.push_back( a.at(i) - a.at(0) );
  }

  if ( _v.getDimension() < (int)_v.size() ) ERROR_LOG << "You are providing too many directions in this dimension - your points must be L.D.";
  if ( 0                == (int)_v.size() ) ERROR_LOG << "You have not provided any directions";


}

///Constuctor that takes between an origin, and between 0 and 
/// (nDim - 1) vectors.
NPlane::NPlane(const HyperPoint& origin, const HyperPointSet& v) :
  _origin ( origin ),
  _v      ( v      )
{
  WELCOME_LOG << "Hello from the NPlane() Constructor";

  if ( _v.getDimension() < (int)_v.size() ) ERROR_LOG << "You are providing too many directions in this dimension - your points must be L.D.";
  if ( 0                == (int)_v.size() ) ERROR_LOG << "You have not provided any directions";

}


///The subspace is defined by (O + x.v_1 + y.v_2 + z.v_3 ...) - a point in this space
///is then defined by the `ParametricPoint' (x,y,z...). This function takes (x,y,z...)
/// and returns (O + x.v_1 + y.v_2 + z.v_3 ...)
HyperPoint NPlane::getParametricPoint(const HyperPoint& t) const{
  if ( t.size() != (int)_v.size() ) ERROR_LOG << "You have not provided the correct number of parameters";

  HyperPoint point = _origin;

  for(unsigned int i = 0; i < _v.size(); i++){
    point = point + _v.at(i)*t.at(i); //std::vector times a double
  }

  return point;

}

///Print out the origin and vectors that define the NPlane
///
///
void NPlane::print(std::ostream& os, int endline) const{
  os << "NPlane: v = ";
  _origin.print(os, 0);

  for (unsigned int i = 0; i < _v.size(); i++){
    os << " + ";
    _v.at(i).print(os, 0);
    os << "t_" << i;
  }
  if (endline) os << std::endl;
}

///Destructor
///
///
NPlane::~NPlane(){

}


