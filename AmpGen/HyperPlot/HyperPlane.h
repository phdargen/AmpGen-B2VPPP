/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * A HyperPlane that is defined by nDim points in space.
 * i.e. a line in 2D, a plane in 3D ...
 *
 **/
 
#ifndef HYPERPLANE_HH
#define HYPERPLANE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/NPlane.h"
#include "AmpGen/HyperPlot/HyperLine.h"

// Root includes
#include "TDecompLU.h"

// std includes

class HyperPlane : public NPlane {

  public:

  //you need #dim points to define a hyperplane. These are passed via a HyperPointSet 

  HyperPlane(const HyperPointSet& a);

  virtual void print(std::ostream& os=std::cout, int endline=1) const;

  //These find plane - point interceptions
  bool pointInPlane(const HyperPoint& hyperPoint) const;
  HyperPoint findPlaneIntersectionParameters(const HyperPoint& hyperPoint, int dimToOmmit = 0) const;
  /**< \todo remember how this works */

  //These find plane - line interceptions

  HyperPoint findIntersectionPoint(const HyperLine& hyperLine) const;
  /**< \todo remember how this works */
  HyperPoint findLineIntersectionParameter  (const HyperLine& hyperLine) const;
  /**< \todo remember how this works */
  HyperPoint findPlaneIntersectionParameters(const HyperLine& hyperLine) const;
  /**< \todo remember how this works */
  TMatrixD   findIntersectionSolution       (const HyperLine& hyperLine) const;
  /**< \todo remember how this works */

  bool operator ==(const HyperPlane& other) const;

  ~HyperPlane();

};



#endif

