/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * A point in multi-dimensional space  
 *
 **/
 
#ifndef HYPERPOINT_HH
#define HYPERPOINT_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h" 
#include "AmpGen/HyperPlot/Weights.h"

// Root includes
#include "TMatrixD.h"
#include "TRandom.h"

// std includes
#include <sstream> 
#include <iomanip>


class HyperPoint : public Weights {

  protected:
  
  std::vector<double> _coords; /**< The coordinates of the point in muli-dimensional spave */
  

  public:

  HyperPoint(int dimension);
  HyperPoint(int dimension, double val);
  HyperPoint(std::vector<double> coords);

  //some constuctors for a set number of dimensions
  HyperPoint(double x1);
  HyperPoint(double x1, double x2);
  HyperPoint(double x1, double x2, double x3);
  HyperPoint(double x1, double x2, double x3, double x4);
  HyperPoint(double x1, double x2, double x3, double x4, double x5);
  HyperPoint(double x1, double x2, double x3, double x4, double x5, double x6);
  HyperPoint(double x1, double x2, double x3, double x4, double x5, double x6, double x7);

  const std::vector<double>& getVector(){return _coords;} /**< Get the std::vector<double> that contains the coordinates */

  HyperPoint linearTransformation( const TMatrixD& matrix );

  virtual void print(std::ostream& os=std::cout, int endline=1) const;

  //std::vector compatibility
  const double& at(int i) const;
  double& at(int i);

  HyperPoint & operator= (const HyperPoint & other);
  HyperPoint operator+ (const HyperPoint & other) const;
  HyperPoint operator- (const HyperPoint & other) const;
  
  HyperPoint operator- (const double & other) const;
  HyperPoint operator+ (const double & other) const;
  HyperPoint operator/ (const double & other) const;
  HyperPoint operator* (const double & other) const;

  bool operator <(const HyperPoint& other) const;
  bool operator >(const HyperPoint& other) const;
  bool operator <=(const HyperPoint& other) const;
  bool operator >=(const HyperPoint& other) const;
  bool operator ==(const HyperPoint& other) const;
  bool operator !=(const HyperPoint& other) const;

  bool allLT  (const HyperPoint& other) const;
  bool allGT  (const HyperPoint& other) const;
  bool allLTOE(const HyperPoint& other) const;
  bool allGTOE(const HyperPoint& other) const;

  bool operator <(const double& other) const;
  bool operator >(const double& other) const;
  bool operator <=(const double& other) const;
  bool operator >=(const double& other) const;
  bool operator ==(const double& other) const;

  double multiplyElements() const;
  

  double dotProduct(const HyperPoint & other) const;
  void fillRandom(double min = -1.0, double max = 1.0);
  HyperPoint project(const HyperPoint & other) const;


  double distanceTo(const HyperPoint & other) const;
  double norm() const;
  double norm2() const;
  int size() const {return (int)_coords.size();}
  /**< get the dimensionality of the HyperPoint */

  int getDimension() const{return size();}
  /**< get the dimensionality of the HyperPoint */
  
  bool compatible(const HyperPoint& other, bool printError = true) const;

  friend std::ostream& operator<<(std::ostream& os, const HyperPoint& point);


  ~HyperPoint();

};



#endif

