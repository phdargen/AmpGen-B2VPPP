/**
  <B>HyperPlot</B>,
  Author: Sam Harnew, sam.harnew@gmail.com ,
  Date: Dec 2015
 

 Binning Algorithms:

  - HyperBinningAlgorithms::SMART             (see HyperBinningMakerSmart for details on the algorithm          )
  - HyperBinningAlgorithms::MINT              (see HyperBinningMakerMint for details on the algorithm           )
  - HyperBinningAlgorithms::MINT_SMART        (see HyperBinningMakerMintSmart for details on the algorithm      )
  - HyperBinningAlgorithms::MINT_RANDOM       (see HyperBinningMakerMintRandomise for details on the algorithm  )
  - HyperBinningAlgorithms::SMART_RANDOM      (see HyperBinningMakerSmartRandomise for details on the algorithm )
  - HyperBinningAlgorithms::LIKELIHOOD        (see HyperBinningMakerLikelihood for details on the algorithm     )
  - HyperBinningAlgorithms::SMART_LIKELIHOOD  (see HyperBinningMakerSmartLikelihood for details on the algorithm)

Binning Algorithm Options:

  - AlgOption::StartDimension     (int dim                  )
  - AlgOption::BinningDimensions  (std::vector<int> dims    )
  - AlgOption::RandomSeed         (int seed                 )
  - AlgOption::MinBinWidth        (double width             )
  - AlgOption::MinBinWidth        (HyperPoint widths        ) 
  - AlgOption::MinBinContent      (double val               )
  - AlgOption::MinShadowBinContent(double val               )
  - AlgOption::UseWeights         (bool   val = true        )
  - AlgOption::UseShadowData      (const HyperPointSet& data)
  - AlgOption::Empty              (                         )
 
 **/
 
#ifndef HYPERHISTOGRAM_HH
#define HYPERHISTOGRAM_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/HistogramBase.h"
#include "AmpGen/HyperPlot/HyperFunction.h"
#include "AmpGen/HyperPlot/HyperName.h"
#include "AmpGen/HyperPlot/BinningBase.h"
#include "AmpGen/HyperPlot/HyperBinning.h"
#include "AmpGen/HyperPlot/HyperBinningDiskRes.h"
#include "AmpGen/HyperPlot/HyperBinningAlgorithms.h"

// Root includes
#include "TRandom.h"

// std includes
#include <iostream>
#include <fstream>
#include <iomanip>

class HyperHistogram : public HistogramBase, public HyperFunction {

  protected:

  BinningBase* _binning; /**< The HyperVolumeBinning used for the HyperHistogram */
  
  HyperHistogram();
  
  //BinningBase& getBinning() { return (*_binning); }  /**< get the HyperVolumeBinning */


  public:
  
  HyperHistogram(
    const HyperCuboid&   binningRange, 
    const HyperPointSet& points, 

    HyperBinningAlgorithms::Alg alg = HyperBinningAlgorithms::MINT, 

    AlgOption opt0 = AlgOption::Empty(),
    AlgOption opt1 = AlgOption::Empty(),  
    AlgOption opt2 = AlgOption::Empty(),
    AlgOption opt3 = AlgOption::Empty(),
    AlgOption opt4 = AlgOption::Empty(),
    AlgOption opt5 = AlgOption::Empty(),
    AlgOption opt6 = AlgOption::Empty(),
    AlgOption opt7 = AlgOption::Empty(),
    AlgOption opt8 = AlgOption::Empty(),
    AlgOption opt9 = AlgOption::Empty()
  );
 
  HyperHistogram(const BinningBase& binning);
  HyperHistogram(TString filename, TString option = "MEMRES READ");
  HyperHistogram(std::vector<TString> filename);
  HyperHistogram(TString targetFilename, std::vector<TString> filename);

  HyperHistogram(const HyperHistogram& other);

  HyperHistogram& operator=(const HyperHistogram& other);

  int getDimension() const;
  
  void setNames( HyperName names );
  HyperName getNames() const;

  int  fill(const HyperPoint& coords, double weight);
  int  fill(const HyperPoint& coords);
  void fill(const HyperPointSet& points);

  virtual void merge( const HistogramBase& other );

  void merge( TString filenameother );

  
  int estimateCapacity(std::vector<TString> filename, TString binningType);


  //Project the HyperHistograms down into one dimension

  void project(TH1D* histogram, const HyperCuboid& cuboid, double content, int dimension) const;
  void project(TH1D* histogram, const HyperVolume& hyperVolume, double content, int dimension) const;
  TH1D project(int dim = 0, int bins = 100, TString name = "projection") const;
  
  void drawProjection    (TString path, int dim = 0, int bins = 100) const;
  void drawAllProjections(TString path, int bins) const;

  void compareProjection    (TString path, int dim, const HyperHistogram& other, int bins = 100) const;
  void compareAllProjections(TString path, const HyperHistogram& other, int bins = 100) const;

  

  HyperHistogram slice(std::vector<int> sliceDims, const HyperPoint& slicePoint) const;
  std::vector<HyperHistogram> slice(std::vector<int> sliceDims, const HyperPointSet& slicePoints) const;

  void draw2DSlice   (TString path, int sliceDimX, int sliceDimY, const HyperPoint& slicePoint, TString options = "") const;
  void draw2DSliceSet(TString path, int sliceDimX, int sliceDimY, int sliceSetDim, int nSlices, const HyperPoint& slicePoint, TString options = "") const;
  void draw2DSliceSet(TString path, int sliceDimX, int sliceDimY, int nSlices, const HyperPoint& slicePoint, TString options = "") const;
  void draw2DSliceSet(TString path, int nSlices, const HyperPoint& slicePoint, TString options = "") const;
  void drawRandom2DSlice(TString path, TRandom* random = gRandom, TString options = "") const;


  HyperCuboid getLimits() const;


  const BinningBase& getBinning() const { return (*_binning); }  /**< get the HyperVolumeBinning */
  
  virtual double getVal(const HyperPoint& point) const;
  std::vector<double> getVal(const HyperPointSet& points) const; 


  virtual double getBinVolume(int bin) const;

  void save(TString filename);

  TString getBinningType(TString filename);

  void load     (TString filename, TString option = "MEMRES READ");
  void loadEmpty(TString filename, TString option = "MEMRES READ", TString binningType = "HyperBinning");


  void setContentsFromFunc(const HyperFunction& func);
  
  void printFull() const;
  
  void saveToTxtFile(TString filename, bool incError = true) const;
  void saveToTxtFileNoLinks(TString filename, bool incError = true) const;


  void draw(TString path, TString options = "");
  void drawDensity(TString path, TString options = "");
  
  void mergeBinsWithSameContent();
  

  virtual ~HyperHistogram();

};


  //mutable double _nIntegrationsWtrick;
  //mutable double _nIntegrationsWOtrick;
  
  //void printOptimisationStatistics();
  //interpolation stuff
  //HyperPointSet makePointsAtGaussianExtremes(const HyperPoint& mean, const HyperPoint& sigmas, double sigma) const;
  //double intgrateGaussianOverHyperCuboid(const HyperPoint& mean , const HyperPoint& sigmas, const HyperCuboid& cuboid) const;
  //double intgrateGaussianOverHyperVolume(const HyperPoint& point, const HyperPoint& sigmas, const HyperVolume& volume) const;
  //double intgrateGaussianOverBin   (const HyperPoint& point, const HyperPoint& sigmas, int bin) const;
  //double gaussianKernal                 (const HyperPoint& point, const HyperPoint& sigmas) const;
  //HyperPoint findAdaptiveSigma(const HyperPoint& point, const HyperPoint& sigmas) const;
  //std::vector<int> findNHighestContributingKernalBins(const HyperPoint& point, const HyperPoint& sigmas, int n) const;
  //double adaptiveGaussianKernal(const HyperPoint& point, double smoothing = 1.0) const;
  //void reweightDatasetWithAdaptiveGaussianKernal(HyperPointSet& points, double smoothing = 1.0) const;



#endif

