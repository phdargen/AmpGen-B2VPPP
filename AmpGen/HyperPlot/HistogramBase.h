/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * The base class for any histogram object (my version of a TH1)  
 *
 **/

 
#ifndef HISTOGRAMBASE_HH
#define HISTOGRAMBASE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/RootPlotter1D.h"
#include "AmpGen/HyperPlot/RootPlotter2D.h"
#include "AmpGen/HyperPlot/MessageService.h"

// Root includes
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"  
#include "TMath.h"

// std includes


class HistogramBase {
  
  protected:
  
  int     _nBins;                   /**< Number of bins in the histogram */
  std::vector<double> _binContents; /**< Bin contents (note that bin nBins is underflow/overflow) */
  std::vector<double> _sumW2;       /**< Sum of weights^2 for each bin (note that bin nBins is underflow/overflow) */

  double _min;         /**< Minimum bin content (for plotting) */
  double _max;         /**< Maximum bin content (for plotting) */
  double _minDensity;  /**< Minimum bin density (bin content / bin volume)  (for plotting) */
  double _maxDensity;  /**< Maximum bin density (bin content / bin volume)  (for plotting) */

  public:

  HistogramBase(int nBins);
  virtual ~HistogramBase();
  
  int checkBinNumber(int bin) const;

  void resetBinContents(int nBins);
  void clear();

  void fillBase(int binNum, double weight);
  
  void setBinContent(int bin, double val);
  void setBinError  (int bin, double val);
  
  virtual void merge( const HistogramBase& other ); 

  double getBinContent(int bin) const;
  double getBinError  (int bin) const;

  int getNBins() const{return _nBins;}
  /**< Get the number of bins in the histogram */

  void divide(const HistogramBase& other);
  void multiply(const HistogramBase& other);
  void add(const HistogramBase& other);
  void minus(const HistogramBase& other);

  void pulls(const HistogramBase& other);
  void pulls(const HistogramBase& other1, const HistogramBase& other2);
  void asymmetry(const HistogramBase& other);
  void asymmetry(const HistogramBase& other1, const HistogramBase& other2);

  void drawPullHistogram(const HistogramBase& other, TString name, int nBins = 50, double pmLimits = 3.5) const;
  double chi2(const HistogramBase& other) const;
  double pvalue(const HistogramBase& other, int ndof = -1) const;
  double chi2sig(const HistogramBase& other, int ndof = -1) const;


  double integral() const;
  double integralError() const;
  
  void randomiseWithinErrors(int seed);

  double getMin() const;
  double getMax() const;

  void setMin(double min){_min = min;}  /**< Set the minimum bin content for plotting */
  void setMax(double max){_max = max;}  /**< Set the maximum bin content for plotting */

  double getMinDensity() const;
  double getMaxDensity() const;
  void setMinDensity(double min){_minDensity = min;}  /**< Set the minimum frequency density for plotting */
  void setMaxDensity(double max){_maxDensity = max;}  /**< Set the maximum frequency density for plotting */

  void saveBase();
  void saveBase(TString filename);

  void loadBase(TString filename);
  
  void normalise(double area = 1.0);
  
  virtual double getBinVolume(int bin) const;
  double getFrequencyDensity(int bin) const;
  
  void reserveCapacity(int nElements); 

  void makeFrequencyDensity();
  
  void print();

};



#endif

