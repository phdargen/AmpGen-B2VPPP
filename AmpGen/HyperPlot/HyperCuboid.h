/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * A multi-dimensional cuboid that is defined by 2 corners (low and high)
 *
 **/
 
#ifndef HYPERCUBOID_HH
#define HYPERCUBOID_HH

class HyperCuboid;

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/HyperPlane.h"
#include "AmpGen/HyperPlot/HyperLine.h"
#include "AmpGen/HyperPlot/HyperVolume.h"

// Root includes

// std includes
#include <bitset>


class HyperCuboid {

  private:
  
  int _dimension;         /**< The dimensionality of the cuboid */

  HyperPoint _lowCorner;  /**< The lower  corner of the cuboid */
  HyperPoint _highCorner; /**< The higher corner of the cuboid */
  
  /** The HyperPlanes that define the faces of the cuboid. Dy default these
  are not filled, but if needed it's useful to cashe them here */
  mutable std::vector<HyperPlane> _faces; 
  
  void updateFaceCash() const;

  public:
    
  const int& getDimension () const {return _dimension;}
  /**< get the dimensionality */

  HyperCuboid(int dimension);
  HyperCuboid(const HyperPoint& lowCorner, const HyperPoint& highCorner);
  HyperCuboid(int dimension, double low, double high);

  void print(std::ostream& os = std::cout, int endline=1) const;

  HyperPoint getCenter() const;
  bool inVolume(const HyperPoint& coords) const;

  const HyperCuboid& inflateCuboid(double percent);

  double volume() const;
  
  double getPositiveIntersectionParameter(const HyperLine& line) const;
  HyperPoint getPositiveIntersectionPoint(const HyperLine& line) const;
  
  double getNegativeIntersectionParameter(const HyperLine& line) const;
  HyperPoint getNegativeIntersectionPoint(const HyperLine& line) const;

  //pick a random point from within the hypercuboid
  HyperPoint getRandomPoint(TRandom* random = gRandom) const;
  HyperPointSet getRandomPoints(int nPoints = 100, TRandom* random = gRandom) const;

  std::vector<HyperPlane> getConnectedHyperPlanes(const HyperPoint& point) const;
  std::vector<HyperPlane> getConnectedHyperPlanes(const HyperPointSet& point) const;
  std::vector<HyperPlane> getBoundaryHyperPlanes() const;

  HyperPoint getOppositePoint(const HyperPoint& point) const;
  HyperPointSet getConnectedVerticies(const HyperPointSet& pointSet) const;
  HyperPointSet getConnectedVerticies(const HyperPoint& point) const;
  HyperPointSet getVertices() const;
  
  HyperPointSet getEdgeCenters() const;


  HyperCuboid project        (std::vector<int> dims) const;
  HyperCuboid projectOpposite(std::vector<int> dims) const;
  bool        inVolume(const HyperPoint& coords, std::vector<int> dims) const;

  const HyperPoint& getLowCorner () const{ return _lowCorner ; }
  /**< return the low HyperPoint corner */
  const HyperPoint& getHighCorner() const{ return _highCorner; }
  /**< return the high HyperPoint corner */
  HyperPoint& getLowCorner () { return _lowCorner ; }
  /**< return the low HyperPoint corner */
  HyperPoint& getHighCorner() { return _highCorner; }
  /**< return the high HyperPoint corner */
  

  HyperCuboid splitAbove(int dimension, double fractionalSplitPoint) const;
  HyperCuboid splitBelow(int dimension, double fractionalSplitPoint) const;
  HyperVolume split(int dimension, double fractionalSplitPoint) const;
  void split(int dimension, double fractionalSplitPoint, HyperVolume& set) const;

  HyperPoint getWidth() const;
  double getWidth(int dim) const;

  bool setCorners(const HyperPoint& lowCorner, const HyperPoint& highCorner);

  bool operator ==(const HyperCuboid& other) const;
  bool operator !=(const HyperCuboid& other) const;

  ~HyperCuboid();

};



#endif

