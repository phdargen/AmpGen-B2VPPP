#include "AmpGen/HyperPlot/HyperLine.h"

///Standard constructor that takes two HyperPoints
///
HyperLine::HyperLine(const HyperPoint& a, const HyperPoint& b) :
  NPlane( HyperPointSet(a, b) )
{ 

  WELCOME_LOG << "Hello from the HyperLine() Constructor";
}

///Print the parameters that define the HyperLine
///
void HyperLine::print(std::ostream& os, int endline) const{

  os << "HyperLine: v = ";
  _origin.print(os, 0);
  os << " + ";
  _v.at(0).print(os, 0);
  os << "t";
  if (endline) os << std::endl;

}

///Destructor
///
HyperLine::~HyperLine() { 
  GOODBYE_LOG << "Goodbye from the HyperLine() Constructor";
}


