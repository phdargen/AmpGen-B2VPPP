#ifndef CYCLIC_PHASE_BINS_HH
#define CYCLIC_PHASE_BINS_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"

#include "TMath.h"


/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Dealing with cyclic bin boundaries is annoying. This class makes life easier.
 *
 *
 **/

class CyclicPhaseBins{

  std::vector<double> _binEdges;
  
  public:
  
  CyclicPhaseBins(std::vector<double> binEdges);
  /**< Constuct a phase binning with the given bin edges. Since the binning is
  cyclic you only need nBins edges. The edges must be given in asending order, 
  and within 2pi of each other */    

  CyclicPhaseBins();
  /**< Empty constuctor i.e. no bins */    
  
  void setUniformCiSiBins(int nBinPairs               );
  /**< set nBinPairs*2 uniformly spaces bins between -pi and pi */ 

  void setBinEdges       (std::vector<double> binEdges);
  /**< Set the bin edges. See constuctor for details */    
  
  int getBinNumber(double phase) const;
  /**< For a given phase, return the bin number. Note that the phase does not have to
  be mapped into any particular range i.e. [-pi, +pi] */

  double getLowBinBoundary (int bin) const;
  /**< Get the low bin edge for given bin number. This just uses the nmbers stored in _binEdges */

  double getHighBinBoundary(int bin) const;
  /**< Get the high bin edge for given bin number. This just uses the nmbers stored in _binEdges */

  double getLowBinBoundary (double phase) const;
  /**< Get the low bin edge for given phase. This appropriately maps the numbers stored in _binEdges
  so that the bin edge is below the phase given, but within 2pi */

  double getHighBinBoundary(double phase) const;
  /**< Get the high bin edge for given phase. This appropriately maps the numbers stored in _binEdges
  so that the bin edge is below the phase given, but within 2pi */

  int getNumBins() const;
  /**< Get number of bins */

  ~CyclicPhaseBins();
  /**< desctructor */
  
};


#endif

