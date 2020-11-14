/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * HyperVolume is made up from a vector of HyperCuboids
 *
 **/

 
#ifndef HYPERVOLUME_HH
#define HYPERVOLUME_HH

class HyperVolume;

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/HyperPointSet.h"
#include "AmpGen/HyperPlot/HyperCuboid.h"

// Root includes

// std includes


class HyperVolume {
  
  private:

  int _dimension;                          /**< Dimensionality of the HyperVolume */

  std::vector<HyperCuboid> _hyperCuboids;  /**< vector containing HyperCuboids that define the HyperVolume */


  public:

  HyperVolume(int dimension); 
  HyperVolume(const HyperPoint& lowCorner, const HyperPoint& highCorner);
  HyperVolume(const HyperCuboid& cuboid);
  HyperVolume(const HyperCuboid& cuboid1, const HyperCuboid& cuboid2);


  const int& getDimension () const {return _dimension;}
  /**< get the dimensionality of the HyperVolume */

  void addHyperCuboid(const HyperPoint& lowCorner, const HyperPoint& highCorner);
  void addHyperCuboid(const HyperCuboid& hyperCuboid);
  void push_back(const HyperCuboid& hyperCuboid);
  
  const std::vector<HyperCuboid>& getHyperCuboids()     const{return _hyperCuboids;      }
  /**< return the std::vector containing the HyperCuboids */

  const HyperCuboid&              getHyperCuboid(int i) const{return _hyperCuboids.at(i);}
  const HyperCuboid&              at            (int i) const{return _hyperCuboids.at(i);}
  HyperCuboid&              at            (int i){return _hyperCuboids.at(i);}

  /**< return one of the HyperCuboids */

  int size() const { return _hyperCuboids.size();} 
  /**< return the number of HyperCuboids that make up the HyperVolume */

  HyperVolume & operator= (const HyperVolume & other);
  HyperVolume operator+ (const HyperVolume & other) const;

  HyperPoint getAverageCenter() const;
  
  void print(std::ostream& os = std::cout, int endline=1) const;

  bool inVolume(const HyperPoint& coords) const;
  bool inVolume(const HyperPointSet& coords) const;
  double volume() const;
  
  HyperVolume slice(const HyperPoint& coords, std::vector<int> dims) const;

  double getMin(int dimension) const;
  double getMax(int dimension) const;
  
  HyperCuboid getLimits() const;


  HyperVolume splitAll(int dimension, double fractionalSplitPoint);

  ~HyperVolume();

};



#endif

