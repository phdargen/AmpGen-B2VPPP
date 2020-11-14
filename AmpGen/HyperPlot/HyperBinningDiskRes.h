/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * HyperBinningDiskRes is a binning scheme where each bin volume is 
 * defined by a HyperVolume.
 *
 **/

/** \class HyperBinningDiskRes

Finding the correct bin number is quite a compuationally
slow process if one has to loop over every bin and
check if a HyperPoint falls within that bin volume. It's not unusual to
have millions of HyperPoints that need to be sorted into tens of thousands of bins. 
This would require billions of calulations.
To speed this process up there is the option to build a hierarchy of bins. A schematic below
shows a 1D example. 

~~~ {.cpp}

       HyperVolume Numbers 

 |-------------0-------------| 

 |------1------|------2------| 

 |--3---|---4--|---5---|--6--|

 |-7-|-8| 

           Bin Numbers

 | 0 | 1|   2  |   3   |  4  |

~~~

Imagine we have a HyperPoint that falls into Bin 0. One would first check if it 
falls into HyperVolume 0 (note the distiction here between Bin/HyperVolume Numbers as
indicated by the figure). First we check if it falls into HyperVolume 0, then HyperVolume 1 or 2, then 
HyperVolume 3 or 4, then HyperVolume 7 or 8. 

In this simple example, it took 7 operations, whereas checking each Bin would have taken 
5. Clearly as the number of bins increases, it becomes computationally much less 
expensive to follow this hierarchy approach.

The 

*/



 
#ifndef HYPERBINNINGDISKRES_HH
#define HYPERBINNINGDISKRES_HH

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
#include "AmpGen/HyperPlot/HyperBinning.h"


// Root includes
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

// std includes
#include <algorithm>
#include <sstream>

class HyperBinningDiskRes : public HyperBinning {

  private:
  
  mutable TFile* _file;
  bool   _writeable;

  mutable TTree* _tree;
  mutable HyperCuboid _cuboid;
  mutable std::vector<int>* _linkedBins;
  mutable int _volumeNumber;
  mutable int _currentEntry;

  mutable TTree* _treePrimVol;
  mutable int _primVolNum;


  protected:

  void getEntry(int volumeNumber) const;

  void loadHyperBinningTree   ();
  void loadPrimaryVolumeTree  ();

  void createHyperBinningTree   ();
  void createPrimaryVolumeTree  ();  

  protected:


  public:
  
  HyperBinningDiskRes();
  
  HyperBinningDiskRes(const HyperBinningDiskRes& other);


  virtual ~HyperBinningDiskRes();

  //Functions we are required to implement from HyperBinning

  virtual void reserveCapacity(int nElements);

  virtual void setDimension(int dim);

  virtual void addPrimaryVolumeNumber(int volumeNumber);
  virtual bool addHyperVolume(const HyperVolume& hyperVolume, std::vector<int> linkedVolumes = std::vector<int>(0, 0));
  
  virtual int getNumHyperVolumes() const;  
  virtual HyperVolume getHyperVolume(int volumeNumber) const; /**< get one of the HyperVolumes */
  virtual std::vector<int> getLinkedHyperVolumes( int volumeNumber ) const;

  virtual int getNumPrimaryVolumes  () const;  
  virtual int getPrimaryVolumeNumber(int i) const;  

  //Functions we are required to implement from BinningBase that were not implemented in HyperBinning

  virtual bool isDiskResident() const;
  virtual TString filename() const;

  virtual void load(TString filename, TString option = "READ");

  virtual BinningBase* clone() const;


};



#endif

