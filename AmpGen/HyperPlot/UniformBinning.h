/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Uniform binning
 *
 **/

/** \class UniformBinning



*/

 
#ifndef UNIFORMBINNING_HH
#define UNIFORMBINNING_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/HyperPointSet.h"
#include "AmpGen/HyperPlot/HyperCuboid.h"
#include "AmpGen/HyperPlot/HyperVolume.h"
#include "AmpGen/HyperPlot/RootPlotter1D.h"
#include "AmpGen/HyperPlot/RootPlotter2D.h"
#include "AmpGen/HyperPlot/HyperName.h"
#include "AmpGen/HyperPlot/BinningBase.h"
#include "AmpGen/HyperPlot/LoadingBar.h"
#include "AmpGen/HyperPlot/CachedVar.h"


// Root includes
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

// std includes
#include <algorithm>
#include <sstream>

class UniformBinning : public BinningBase {

  private:
  
  HyperCuboid      _limits;
  std::vector<int> _nLocalBins; 

  public:
  
  UniformBinning(HyperCuboid limits, int nLocalBins);
  UniformBinning(HyperCuboid limits, std::vector<int> nLocalBins);
  
  virtual ~UniformBinning();

  int getNumLocalBins(int dimension) const;

  int getGlobalBinNumber( std::vector<int> localBinNumbers ) const;
  std::vector<int> getLocalBinNumbers( int globalBinNumber ) const;

  int getLocalBinNumber(int dim, double val) const;
  std::vector<int> getLocalBinNumbers(const HyperPoint& coords) const;

  double getLowBinEdgeLocal(int dim, int localBinNum) const;
  double getHighBinEdgeLocal(int dim, int localBinNum) const;

  HyperPoint getLowCorner(int globalBinNum) const;
  HyperPoint getHighCorner(int globalBinNum) const;

  /* Virtual functions that need to implemented from BinningBase */
  /*    These will be implemented in the derived classes */

  virtual void load(TString filename, TString option = "READ");
  virtual BinningBase* clone() const;

  virtual void save(TString filename) const;
  virtual void save() const; 

  virtual void mergeBinnings( const BinningBase& other );

  virtual int getNumBins() const;
  virtual int getBinNum(const HyperPoint& coords) const;

  virtual HyperVolume getBinHyperVolume(int binNumber) const;

  virtual HyperPoint  getAverageBinWidth() const;
  virtual HyperCuboid getLimits()          const;



};



#endif

