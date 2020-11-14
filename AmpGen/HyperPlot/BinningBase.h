/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * The base class for any Binning object 
 *
 **/

 
#ifndef BINNINGBASE_HH
#define BINNINGBASE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HyperPoint.h"
#include "AmpGen/HyperPlot/HyperVolume.h"
#include "AmpGen/HyperPlot/HyperName.h"


// Root includes


// std includes


class BinningBase {
  
  private:

  int       _dimension;
  /**< 
    Dimensionality of the binning. Can only be set ONCE with setDimension();
  */  

  HyperName _axisNames;
  /**< 
    Possible to assign names to each dimension, which will be used for plotting etc
  */  

  TString _binningType; 
  /**< 
    Useful to know what kind of binning this actually is (will be set by derrived type)
  */  

  protected:
  
  virtual void setDimension (int dimension);  /**< set the dimensionality of the binning */
  void setBinningType(TString binningType);

  public:

  BinningBase();

  void      setNames( HyperName names );  /**< set the a name  for each dimension */
  HyperName getNames() const;             /**< get the a names for each dimension */
  const int& getDimension () const;       /**< get the dimensionality of the binning */

  double getMin(int dimension) const;
  double getMax(int dimension) const;
  
  TString getBinningType() const;
  bool isSameBinningType(const BinningBase& other) const;

  virtual bool isDiskResident() const;
  virtual TString filename() const;

  virtual ~BinningBase();
  
  //Purely virtual functions

  virtual void load(TString filename, TString option = "READ") = 0;

  virtual void save(TString filename) const = 0;
  virtual void save() const = 0; 

  virtual BinningBase* clone() const = 0;
  
  virtual int getNumBins() const = 0;
  virtual int getBinNum(const HyperPoint& coords) const = 0;

  virtual HyperVolume getBinHyperVolume(int binNumber) const = 0;

  virtual void mergeBinnings( const BinningBase& other ) = 0;

  virtual HyperPoint getAverageBinWidth() const = 0;
  virtual HyperCuboid getLimits() const = 0;

  virtual void reserveCapacity(int nElements); 
  
  virtual std::vector<int> getBinNum(const HyperPointSet& coords) const;


};



#endif

