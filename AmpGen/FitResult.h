#ifndef AMPGEN_FITRESULT_H
#define AMPGEN_FITRESULT_H

#include "TMatrixD.h"
#include "TClass.h"
#include "TFile.h"
#include "TTree.h"

#include "AmpGen/FitFraction.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/Utilities.h"

namespace AmpGen
{
  class Minimiser;
  class LinearErrorPropagator; 

  class FitResult
  {
  public:
    ~FitResult(){};
    FitResult(); 
    explicit FitResult( const FitResult& other );
    explicit FitResult( const std::string& filename );
    explicit FitResult( const Minimiser& mini );
    FitResult( MinuitParameterSet& mps, const TMatrixD& covMini );

    void addObservable( const std::string& name, const double& F );
    void addChi2( const double& chi2, const double& nBins );
    void addFractions( const std::vector<FitFraction>& fractions );
    void addFraction( const std::string& name, const double& frac, const double& err );
    void setCov( const size_t& x, const size_t& y, const double& F );
    void setSystematic( const std::string& sys ){ m_sys = sys;}  
    void writeToFile( const std::string& fname );
    void writeToFileMod( const std::string& fname );
    void writeToOptionsFile( const std::string& fname, int fixParams = 0 );
    void writeToRootFile( TFile * output,  unsigned seed = 0, int verbose = 0, double nll = 0, unsigned numAmps = 0, double Ns = 0, std::vector<double> thresholds = {}, std::vector<double> numFracAboveThresholds = {} );

    void plotSpline(const std::string& name, const std::string& outDir = ".");
      
    std::string latexName(std::string name);
    void printToLatexTable( const std::string& fname );

    void clearFitFractions();

    bool readFile( const std::string& fname, bool verbose = false );

    double chi2() const;
    double LL()   const;
    double dof()  const;
    double cov(const size_t& x, const size_t& y ) const;
    double cov(const std::string& x, const std::string& y ) const; 
    double correlation( const std::string& x, const std::string& y ) const;
    
    int status() const;
    int nParam() const;
    int nBins()  const; 

    std::map<std::string, double> observables() const;
    MinuitParameterSet* mps()   const;

    std::vector<FitFraction> fitFractions()     const;
    std::vector<MinuitParameter*> parameters()  const;
    std::vector<MinuitParameter*> floating(const bool& extended = false) const;
    TMatrixD cov() const;
    
    void print() const;

    TMatrixD getReducedCovariance( const bool& extended = false ) const;
    LinearErrorPropagator getErrorPropagator( const bool& extended = false ) const;

  private:
    MinuitParameterSet*                 m_mps    = {nullptr};
    double                              m_chi2   = {0};
    double                              m_LL     = {-999};
    double                              m_nBins  = {0};
    double                              m_nParam = {0};
    int                                 m_status = {-1};
    bool                                m_fitted = {false};
    std::string                         m_sys    = {""};  
    std::map<std::string, double>       m_observables;
    std::vector<FitFraction>            m_fitFractions;
    TMatrixD                            m_covarianceMatrix;
    std::map<std::string, unsigned int> m_covMapping;

    void addToParameters( const std::string& line );
    void addToObservables( const std::string& line );
    void setFitQuality( const std::string& line );
  };
} // namespace AmpGen

#endif
