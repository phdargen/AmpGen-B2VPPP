/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Holds a n-dim vector of names that can be associtated
 *  to a n-dim Histogram (for axis titles etc). Also possible
 *  to store the units.
 **/

 
#ifndef HYPERNAME_HH
#define HYPERNAME_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"

// Root includes
#include "TString.h"

// std includes


class HyperName {

  std::vector<TString> _names; /**< vector of names */
  std::vector<TString> _units; /**< vector of units */

  public:
  
  HyperName(int dim); 
  
  HyperName(std::vector<TString> names);

  HyperName(TString name0);
  HyperName(TString name0, TString name1);
  HyperName(TString name0, TString name1, TString name2);
  HyperName(TString name0, TString name1, TString name2, TString name3);
  HyperName(TString name0, TString name1, TString name2, TString name3, TString name4);
  
  void setUnits(TString units);
  void setUnits(TString units, int dim);

  TString getAxisString(int dim);
  
  int getDimension() const;

  HyperName slice(std::vector<int> sliceDims);

  const TString& at(int dim) const{return _names.at(dim);} /**< access an element of the name vector */
  TString& at(int dim) {return _names.at(dim);}            /**< access an element of the name vector */

  ~HyperName(); 

}; 


#endif

