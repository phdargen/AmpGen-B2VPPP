/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * A HyperLine that is defined by 2 points nDim space.
 *
 **/
 
#ifndef HYPERLINE_HH
#define HYPERLINE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/NPlane.h"

// Root includes

// std includes


class HyperLine : public NPlane {
  
  public:

  HyperLine(const HyperPoint& a, const HyperPoint& b);
  
  virtual void print(std::ostream& os=std::cout, int endline=1) const;

  //const HyperPoint& getDirection() const{ return this->getDirection(0); }

  ~HyperLine();

};



#endif

