/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
**/
#ifndef NPLANE_HH
#define NPLANE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/HyperPointSet.h"

// Root includes

// std includes


class HyperPointSet;

/**
 *
 * A NPlane is some subspace of nDimensional space
 * defined by between 1 and nDim HyperPoints.  
 * e.g. to 3D space, you could have a point, a line, or a plane 
 * (defined by 1, 2 and 3 points respectively)
 *
 **/

class NPlane {

  protected:
  
  HyperPoint _origin; /**< The origin of the subspace */

  HyperPointSet _v;   /**< The vectors that define the dimensionality of the subspace */

  public:

  //you need #dim points to define a hyperplane. These are passed via a HyperPointSet 

  NPlane(const HyperPointSet& a);
  NPlane(const HyperPoint& origin, const HyperPointSet& v);

  int getDimension () const {return _origin.getDimension();}
  /**< get the dimensionality */
  int getN         () const {return _v     .size();}
  /**< get the number of vectors that define the NPlane */


  const HyperPoint& getOrigin()    const{ return _origin; }
  /**< get the origin of the NPlane */

  virtual const HyperPoint& getDirection(int i = 0) const{ return _v.at(i); }
  /**< get one of the vectors that defines the NPlane */

  HyperPoint getParametricPoint(const HyperPoint& t) const;

  virtual void print(std::ostream& os=std::cout, int endline=1) const;

  virtual ~NPlane();

};


#endif

