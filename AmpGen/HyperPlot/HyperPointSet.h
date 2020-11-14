/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * A vector of HyperPoints 
 *
 **/

 
#ifndef HYPERPOINTSET_HH
#define HYPERPOINTSET_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/StatisticsFinder.h"

// Root includes
#include "TMatrixD.h"
#include "TFile.h"
#include "TTree.h"

// std includes
#include <algorithm>
#include <functional>
#include <array>
#include <iostream>


class HyperPointSet {
  
  int _dimension;                  /**< The dimensionality of HyperPoints in the HyperPointSet */
  std::vector<HyperPoint> _points; /**< std::vector containing the HyperPoints */

  void save();

  public:

  const int& getDimension() const{return _dimension;} /**< Get the dimensionality of HyperPoints in the HyperPointSet */
  
  bool compatible(const HyperPoint& other, bool printError = true) const;

  const HyperPoint& at(int i) const      ;
  HyperPoint& at(int i)                  ;
  unsigned int size() const              ;
  void push_back(const HyperPoint& point);

  void addHyperPointSet(const HyperPointSet& other);
  void removeDuplicates();
  void sort();
  
  double getSumW() const;
  double getSumW2() const;

  double   getCorrelation(int i, int j) const;
  TMatrixD getCorrelationMatrix() const;

  double   getCovarience(int i, int j) const;
  TMatrixD getCovarienceMatrix() const;
  
  HyperPoint getMin() const;
  HyperPoint getMax() const;

  HyperPoint mean()          const;
  HyperPoint geometricMean() const; 
  HyperPoint harmonicMean()  const;
  
  void save(TString path);
  void load(TString path);

  void print(std::ostream& os = std::cout) const;

  HyperPointSet(int dimension);
  HyperPointSet(const HyperPoint& point);
  HyperPointSet(const HyperPoint& point1, const HyperPoint& point2);
  HyperPointSet(const HyperPoint& point1, const HyperPoint& point2, const HyperPoint& point3);
  HyperPointSet(int npoints, const HyperPoint& point);

  HyperPointSet(TString path);

  virtual ~HyperPointSet();

  //These last functions are more std::vector oriented (i.e. treating it as a set of HyperVectors 
  //rather than a set of HyperPoints)

  bool linearlyIndependant() const;  //check that the std::vectors are L.I.

  // gramSchmidtProcess makes the hyperPoints into an orthanormal basis. 
  // Note that only the first dim points will be used.

  HyperPointSet gramSchmidtProcess() const; 

  void normalise();

};


#endif

