/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * HyperFunction takes a HyperPoint and returns
 * a double. This can be used to reweight HyperPointSets etc.
 *
 **/

 
#ifndef HYPERFUNCTION_HH
#define HYPERFUNCTION_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPointSet.h"
#include "AmpGen/HyperPlot/HyperCuboid.h"
#include "AmpGen/HyperPlot/RootPlotter2D.h"

// Root includes
#include "TH2D.h"

// std includes


class HyperFunction{  

  private:
  
  HyperCuboid _limits;

  public:

  HyperFunction(); /**< Constructor */
  HyperFunction(const HyperCuboid& limits); /**< Constructor */

  virtual double getVal(const HyperPoint& point) const = 0; /**< Virtual function that defines a HyperFunction (Map from HyperPoint -> double) */
  
  void reweightDataset(HyperPointSet& points);
  
  void setFuncLimits(const HyperCuboid& limits);

  const HyperCuboid& getFuncLimits() const{return _limits;}

  TH2D make2DFuncSlice   (TString name, int sliceDimX, int sliceDimY, const HyperPoint& slicePoint, int nbins = 100) const;
  void draw2DFuncSlice   (TString path, int sliceDimX, int sliceDimY, const HyperPoint& slicePoint, int nbins = 100) const;
  void draw2DFuncSliceSet(TString path, int sliceDimX, int sliceDimY, int sliceSetDim, int nSlices, const HyperPoint& slicePoint, int nbins = 100) const;
  void draw2DFuncSliceSet(TString path, int sliceDimX, int sliceDimY, int nSlices, const HyperPoint& slicePoint, int nbins = 100) const;
  void draw2DFuncSliceSet(TString path, int nSlices, const HyperPoint& slicePoint, int nbins = 100) const;
  
  double getDifference(const HyperFunction& other, const HyperPoint& point);
  /**< Get the difference between this HyperFunction and another HyperFunction at a point
  in n-dim space. (diff = this - other) */

  void fillCorrelations(TH2D& hist, const HyperFunction& other, const HyperPointSet& points);

  virtual ~HyperFunction(){;} /**< Destructor */

};


#endif
