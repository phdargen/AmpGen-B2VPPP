#include <TMatrixDfwd.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixTUtils.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <memory.h>
#include <iomanip>
#include <istream>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FitFraction.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Utilities.h"

using namespace AmpGen;

FitResult::FitResult() = default; 

FitResult::FitResult( const FitResult& other )
  : m_mps( other.mps()  )
  , m_chi2( other.chi2() )
  , m_LL( other.LL() )
  , m_nBins( other.nBins() )
  , m_nParam( other.nParam() )
  , m_status( other.status() )
  , m_observables( other.observables() )
  , m_fitFractions( other.fitFractions() )
  , m_covarianceMatrix(other.cov())
  , m_covMapping(other.m_covMapping )
{
}

FitResult::FitResult( const std::string& filename ) :
  m_mps( new MinuitParameterSet() ) 
{
  m_fitted = readFile( filename );
}

FitResult::FitResult( const Minimiser& mini )
  : m_mps  ( mini.parSet() )
  , m_LL   ( mini.FCN() )
  , m_nParam( 0 )
  , m_status( mini.status() )
  , m_covarianceMatrix( mini.covMatrixFull() ) {
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->at(i)->isFree() ) m_nParam++;
    m_covMapping[ m_mps->at(i)->name() ] = i;
  }
}

FitResult::FitResult( MinuitParameterSet& mps, const TMatrixD& covMini ) : m_mps(&mps)
{
  if ( int( mps.size() ) != covMini.GetNcols() ) {
    ERROR( "Minuit parameter set size does not match covariance matrix size!" );
  }
  m_covarianceMatrix.ResizeTo( covMini.GetNcols(), covMini.GetNrows() );
  m_covarianceMatrix = covMini;
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    if ( m_mps->at(i)->isFree()) m_nParam++;
  }
}

bool FitResult::readFile( const std::string& fname, bool verbose )
{
  std::ifstream checkIsClosed( fname );
  if ( !checkIsClosed.is_open() 
      || checkIsClosed.peek() == std::ifstream::traits_type::eof() ) return false;
  checkIsClosed.close();
  auto lines = vectorFromFile(fname);
  if( *lines.rbegin() != "End Log" ){
    ERROR("File not properly closed: " << *lines.rbegin() );
    return false; 
  }
  INFO("Reading fit results from file " << fname);   
  std::vector<std::string> parameterLines;
  for( auto& line : lines ){
    const std::string name = split( line, ' ' )[0];
    if ( name == "Parameter" ) parameterLines.push_back( line );
    else if ( name == "FitQuality" )  this->setFitQuality( line );
    else if ( name == "FitFraction" ) this->m_fitFractions.emplace_back(line);
    else if ( name == "Observable" )  this->addToObservables( line );
    else if ( name == "Systematic" )  this->setSystematic( line );
  }
  size_t nParameters = parameterLines.size();
  m_covarianceMatrix.ResizeTo( parameterLines.size(), parameterLines.size() );
  for (size_t i = 0; i < nParameters; ++i ) {
    auto tokens             = split( parameterLines[i], ' ' );
    m_covMapping[tokens[1]] = i;
    if(verbose)INFO(tokens[1]);
    m_mps->add( new MinuitParameter( tokens[1], parse<Flag>(tokens[2]), stod( tokens[3] ), stod( tokens[4] ), 0, 0 ) );
    for (size_t j = 0; j < nParameters; ++j ) m_covarianceMatrix( i, j ) = stod( tokens[5 + j] );
  }
      
  return true;
}

void FitResult::setFitQuality( const std::string& line )
{
  auto tokens = split( line, ' ' );
  bool status=true;
  if( tokens.size() != 6 ){
    WARNING("Cannot pass FitQuality line: " << line );
    return;
  }
  m_chi2      = lexical_cast<double>( tokens[1] , status );
  m_nBins     = lexical_cast<double>( tokens[2] , status );
  m_nParam    = lexical_cast<double>( tokens[3] , status );
  m_LL        = lexical_cast<double>( tokens[4] , status );
  m_status    = lexical_cast<int>   ( tokens[5] , status );
}

void FitResult::addToObservables( const std::string& line )
{
  auto tokens              = split( line, ' ' );
  m_observables[tokens[1]] = stod( tokens[2] );
}

void FitResult::addObservable( const std::string& name, const double& F ) { m_observables[name] = F; }

void FitResult::writeToFileMod( const std::string& fname )
{
  std::ofstream outlog;
  outlog.open( fname );
  /*
  for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
    auto param = m_mps->at(i);
    outlog << "Parameter"
      << " " << param->name() << " " << to_string<Flag>(param->flag()) << " " << param->mean() << " "
      << ( param->isFree() ? m_mps->at(i)->err() : 0 ) << " ";
    for (size_t j = 0; j < (size_t)m_covarianceMatrix.GetNcols(); ++j ) outlog << m_covarianceMatrix[i][j] << " ";
    outlog << std::endl;
  }
  */
  outlog << std::setprecision( 8 );
  outlog << "FitQuality " << m_chi2 << " " << m_nBins << " " << m_nParam << " " << m_LL    << " " << m_status << "\n";
  for ( auto& f : m_fitFractions )  outlog << "FitFraction " << f.name() << " " << f.val() << " +/- " << f.err()  << "\n";
  //for ( auto& o : m_observables )   outlog << "Observable "  << o.first  << " " << o.second << "\n";
  for (size_t i = 0; i < (size_t)m_mps->size(); ++i ) {
        auto param = m_mps->at(i);
        if(!param->isFree() || param->name().find( "::Spline" ) != std::string::npos )continue;
        if(param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos) continue;
        outlog << "Parameter "  << param->name()  << " " << param->mean() << " +/- " << param->err() <<  "\n";
    }

  outlog << "End Log\n";
  outlog.close();
}

void FitResult::writeToFile( const std::string& fname )
{
    std::ofstream outlog;
    outlog.open( fname );
    for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
        auto param = m_mps->at(i);
        outlog << "Parameter"
        << " " << param->name() << " " << to_string<Flag>(param->flag()) << " " << param->mean() << " "
        << ( param->isFree() ? m_mps->at(i)->err() : 0 ) << " ";
        for (size_t j = 0; j < (size_t)m_covarianceMatrix.GetNcols(); ++j ) outlog << m_covarianceMatrix[i][j] << " ";
        outlog << std::endl;
    }
    outlog << std::setprecision( 8 );
    outlog << "FitQuality " << m_chi2 << " " << m_nBins << " " << m_nParam << " " << m_LL    << " " << m_status << "\n";
    for ( auto& f : m_fitFractions )  outlog << "FitFraction " << f.name() << " " << f.val() << " " << f.err()  << "\n";
    for ( auto& o : m_observables )   outlog << "Observable "  << o.first  << " " << o.second << "\n";
    outlog << "Systematic " << m_sys << "\n";
    outlog << "End Log\n";
    outlog.close();
}

std::string FitResult::latexName(std::string name){
    
    TString n(name);

    n.ReplaceAll("{","\\{");
    n.ReplaceAll("}","\\}");
    
    n.ReplaceAll("rhoOmega10","rho(770)0");
    n.ReplaceAll("rhoOmega20","rho(770)0");
    n.ReplaceAll("rhoOmega30","rho(770)0");
    n.ReplaceAll("rhoOmega40","rho(770)0");
    n.ReplaceAll("rhoOmega00","rho(770)0");
    
    n.ReplaceAll("*","^{*}");
    n.ReplaceAll("+","^{+}");
    n.ReplaceAll("-","^{-}");
    n.ReplaceAll(")0",")^{0}");
    n.ReplaceAll(",","");
    n.ReplaceAll(";","");

    n.ReplaceAll("Sum_","\\text{Sum } ");

    n.ReplaceAll("J/psi0","J/\\psi");
    n.ReplaceAll("psi(2S)0","\\psi(2S)");
    n.ReplaceAll("rho","\\rho");
    n.ReplaceAll("pi","\\pi");
    
    n.ReplaceAll("[GounarisSakurai.Omega]","");
    n.ReplaceAll("[GounarisSakurai]","");
    n.ReplaceAll("[LASS]","");
    n.ReplaceAll("[GLASS]","");
    n.ReplaceAll("[Bugg]","");
    n.ReplaceAll("[Kappa]","");
    n.ReplaceAll("[Dabba]","");
    n.ReplaceAll("[GSpline.EFF]","");
    n.ReplaceAll("[GSpline]","");
    n.ReplaceAll("[SBW]","");
    n.ReplaceAll("[ExpNR]","");
    n.ReplaceAll("[Flatte]","");
    n.ReplaceAll("[DGSpline]","");

    n.ReplaceAll("GounarisSakurai.Omega","");
    n.ReplaceAll("GounarisSakurai","");
    n.ReplaceAll("LASS","");
    n.ReplaceAll("GLASS","");
    n.ReplaceAll("Bugg","");
    n.ReplaceAll("Kappa","");
    n.ReplaceAll("Dabba","");
    n.ReplaceAll("GSpline.EFF","");
    n.ReplaceAll("SBW","");
    n.ReplaceAll("ExpNR","");

    n.ReplaceAll("(0)","_{0}");
    n.ReplaceAll("(1)","_{1}");
    n.ReplaceAll("(2)","_{2}");

    n.ReplaceAll("NonResP20","NonResP0");    
    n.ReplaceAll("NonResP30","NonResP0");    
    n.ReplaceAll("NonResP40","NonResP0");    
    n.ReplaceAll("NonResP50","NonResP0");    
    n.ReplaceAll("NonResS20","NonResS0");    
    n.ReplaceAll("NonResS30","NonResS0");    
    n.ReplaceAll("NonResS40","NonResS0");    
    n.ReplaceAll("NonResS50","NonResS0");    
    n.ReplaceAll("NonResV20","NonResV0");    
    n.ReplaceAll("NonResV30","NonResV0");    
    n.ReplaceAll("NonResV40","NonResV0");    
    n.ReplaceAll("NonResV50","NonResV0");    
    n.ReplaceAll("NonResA20","NonResA0");    
    n.ReplaceAll("NonResA30","NonResA0");    
    n.ReplaceAll("NonResA40","NonResA0");    
    n.ReplaceAll("NonResA50","NonResA0");    
    n.ReplaceAll("NonResT20","NonResT0");    
    n.ReplaceAll("NonResT30","NonResT0");    
    n.ReplaceAll("NonResT40","NonResT0");    
    n.ReplaceAll("NonResT50","NonResT0");    
    n.ReplaceAll("NonResPT20","NonResPT0");    
    n.ReplaceAll("NonResPT30","NonResPT0");    
    n.ReplaceAll("NonResPT40","NonResPT0");    
    n.ReplaceAll("NonResPT50","NonResPT0");   

    n.ReplaceAll("KPi10[FOCUS.K\\pi]\\{K^{+}\\pi^{-}\\}","[K^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("KPi10","[K^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("KPi20","[K^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("KPi30","[K^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("KPi40","[K^{+}\\pi^{-}]_{\\text{S}}");

    n.ReplaceAll("PiPi00","[\\pi^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("PiPi10","[\\pi^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("PiPi20","[\\pi^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("PiPi30","[\\pi^{+}\\pi^{-}]_{\\text{S}}");
    n.ReplaceAll("PiPi40","[\\pi^{+}\\pi^{-}]_{\\text{S}}");

    n.ReplaceAll("\\{NonResS0\\{psi(2S)^{0}\\pi^{+}\\}","\\{[psi(2S)^{0}\\pi^{+}\\}_{\\text{S}}");    
    n.ReplaceAll("\\{NonResP0\\{psi(2S)^{0}\\pi^{+}\\}","\\{[psi(2S)^{0}\\pi^{+}\\}_{\\text{P}}");    
    n.ReplaceAll("\\{NonResA0\\{psi(2S)^{0}\\pi^{+}\\}","\\{[psi(2S)^{0}\\pi^{+}\\}_{\\text{A}}");    
    n.ReplaceAll("\\{NonResV0\\{psi(2S)^{0}\\pi^{+}\\}","\\{[psi(2S)^{0}\\pi^{+}\\}_{\\text{V}}");    
    n.ReplaceAll("\\{NonResT0\\{psi(2S)^{0}\\pi^{+}\\}","\\{[psi(2S)^{0}\\pi^{+}\\}_{\\text{T}}");    
    n.ReplaceAll("\\{NonResPT0\\{psi(2S)^{0}\\pi^{+}\\}","\\{[psi(2S)^{0}\\pi^{+}\\}_{\\text{PT}}");  
    
    n.ReplaceAll("\\{NonResS0\\{psi(2S)^{0}\\pi^{-}\\}","\\{[psi(2S)^{0}\\pi^{-}\\}_{\\text{S}}");    
    n.ReplaceAll("\\{NonResP0\\{psi(2S)^{0}\\pi^{-}\\}","\\{[psi(2S)^{0}\\pi^{-}\\}_{\\text{P}}");    
    n.ReplaceAll("\\{NonResA0\\{psi(2S)^{0}\\pi^{-}\\}","\\{[psi(2S)^{0}\\pi^{-}\\}_{\\text{A}}");    
    n.ReplaceAll("\\{NonResV0\\{psi(2S)^{0}\\pi^{-}\\}","\\{[psi(2S)^{0}\\pi^{-}\\}_{\\text{V}}");    
    n.ReplaceAll("\\{NonResT0\\{psi(2S)^{0}\\pi^{-}\\}","\\{[psi(2S)^{0}\\pi^{-}\\}_{\\text{T}}");    
    n.ReplaceAll("\\{NonResPT0\\{psi(2S)^{0}\\pi^{-}\\}","\\{[psi(2S)^{0}\\pi^{-}\\}_{\\text{PT}}");    
    
    n.ReplaceAll("\\{NonResS0\\{psi(2S)^{0}K^{+}\\}","\\{[psi(2S)^{0}K^{+}\\}_{\\text{S}}");    
    n.ReplaceAll("\\{NonResP0\\{psi(2S)^{0}K^{+}\\}","\\{[psi(2S)^{0}K^{+}\\}_{\\text{P}}");    
    n.ReplaceAll("\\{NonResA0\\{psi(2S)^{0}K^{+}\\}","\\{[psi(2S)^{0}K^{+}\\}_{\\text{A}}");    
    n.ReplaceAll("\\{NonResV0\\{psi(2S)^{0}K^{+}\\}","\\{[psi(2S)^{0}K^{+}\\}_{\\text{V}}");    
    n.ReplaceAll("\\{NonResT0\\{psi(2S)^{0}K^{+}\\}","\\{[psi(2S)^{0}K^{+}\\}_{\\text{T}}");    
    n.ReplaceAll("\\{NonResPT0\\{psi(2S)^{0}K^{+}\\}","\\{[psi(2S)^{0}K^{+}\\}_{\\text{PT}}");   
    
    n.ReplaceAll("\\rho(770)^{0}\\{\\pi^{+}\\pi^{-}\\}","\\rho(770)^{0}");    
    n.ReplaceAll("\\rho(1450)^{0}\\{\\pi^{+}\\pi^{-}\\}","\\rho(1450)^{0}");    
    n.ReplaceAll("K^{*}(892)^{0}\\{K^{+}\\pi^{-}\\}","K^{*}(892)^{0}");    
    n.ReplaceAll("K^{*}(1410)^{0}\\{K^{+}\\pi^{-}\\}","K^{*}(1410)^{0}");    

    n.ReplaceAll("\\{NonResP0\\{psi(2S)^{0}\\rho(770)^{0}\\}","\\{[psi(2S)^{0}\\rho(770)^{0}\\}_{\\text{P}}");    
    n.ReplaceAll("\\{NonResS0\\{psi(2S)^{0}\\rho(770)^{0}\\}","\\{[psi(2S)^{0}\\rho(770)^{0}\\}_{\\text{S}}");    
    n.ReplaceAll("\\{NonResV0\\{psi(2S)^{0}\\rho(770)^{0}\\}","\\{[psi(2S)^{0}\\rho(770)^{0}\\}_{\\text{V}}");    
    n.ReplaceAll("\\{NonResA0\\{psi(2S)^{0}\\rho(770)^{0}\\}","\\{[psi(2S)^{0}\\rho(770)^{0}\\}_{\\text{A}}");    
    n.ReplaceAll("\\{NonResT0\\{psi(2S)^{0}\\rho(770)^{0}\\}","\\{[psi(2S)^{0}\\rho(770)^{0}\\}_{\\text{T}}");    
    n.ReplaceAll("\\{NonResPT0\\{psi(2S)^{0}\\rho(770)^{0}\\}","\\{[psi(2S)^{0}\\rho(770)^{0}\\}_{\\text{PT}}");    
   
    n.ReplaceAll("\\{NonResP0\\{psi(2S)^{0}K^{*}(892)^{0}\\}","\\{[psi(2S)^{0}K^{*}(892)^{0}\\}_{\\text{P}}");    
    n.ReplaceAll("\\{NonResS0\\{psi(2S)^{0}K^{*}(892)^{0}\\}","\\{[psi(2S)^{0}K^{*}(892)^{0}\\}_{\\text{S}}");    
    n.ReplaceAll("\\{NonResV0\\{psi(2S)^{0}K^{*}(892)^{0}\\}","\\{[psi(2S)^{0}K^{*}(892)^{0}\\}_{\\text{V}}");    
    n.ReplaceAll("\\{NonResA0\\{psi(2S)^{0}K^{*}(892)^{0}\\}","\\{[psi(2S)^{0}K^{*}(892)^{0}\\}_{\\text{A}}");    
    n.ReplaceAll("\\{NonResT0\\{psi(2S)^{0}K^{*}(892)^{0}\\}","\\{[psi(2S)^{0}K^{*}(892)^{0}\\}_{\\text{T}}");    
    n.ReplaceAll("\\{NonResPT0\\{psi(2S)^{0}K^{*}(892)^{0}\\}","\\{[psi(2S)^{0}K^{*}(892)^{0}\\}_{\\text{PT}}");   
    
    n.ReplaceAll(")S0",")_{\\text{S}}^{0}");
    n.ReplaceAll(")P0",")_{\\text{P}}^{0}");
    n.ReplaceAll(")V0",")_{\\text{V}}^{0}");
    n.ReplaceAll(")A0",")_{\\text{A}}^{0}");
    n.ReplaceAll(")T0",")_{\\text{T}}^{0}");
    n.ReplaceAll(")PT0",")_{\\text{PT}}^{0}");
    
    n.ReplaceAll(")S+",")_{\\text{S}}^{+}");
    n.ReplaceAll(")P+",")_{\\text{P}}^{+}");
    n.ReplaceAll(")V+",")_{\\text{V}}^{+}");
    n.ReplaceAll(")A+",")_{\\text{A}}^{+}");
    n.ReplaceAll(")T+",")_{\\text{T}}^{+}");
    n.ReplaceAll(")PT+",")_{\\text{PT}}^{+}");
    
    n.ReplaceAll(")S-",")_{\\text{S}}^{-}");
    n.ReplaceAll(")P-",")_{\\text{P}}^{-}");
    n.ReplaceAll(")V-",")_{\\text{V}}^{-}");
    n.ReplaceAll(")A-",")_{\\text{A}}^{-}");
    n.ReplaceAll(")T-",")_{\\text{T}}^{-}");
    n.ReplaceAll(")PT-",")_{\\text{PT}}^{-}");
    

    n.ReplaceAll("_Re","\\,\\text{Amp}");
    n.ReplaceAll("_Im","\\,\\text{Phase}");

    n.ReplaceAll("psi(2S)^{0}","\\psi(\\text{2S})\\,");
    n.ReplaceAll("psi(4","\\psi(4");
    n.ReplaceAll("\\{","\\rightarrow [");
    n.ReplaceAll("\\}","]");
    
    n.ReplaceAll("[","\\left[");
    n.ReplaceAll("]","\\right]");

    n.ReplaceAll("_mass"," \\text{ mass [MeV]}");
    n.ReplaceAll("_width"," \\text{ width [MeV]}");
    
    return (std::string)n;
}

void FitResult::printToLatexTable( const std::string& fname )
{
    std::ofstream SummaryFile;
    SummaryFile.open(fname.c_str(),std::ofstream::trunc);
    
    SummaryFile << "\\begin{tabular}{l r} " << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\multicolumn{1}{c}{Decay Channel} & \\multicolumn{1}{c}{$F_{i}[\\%]$} " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << std::fixed << std::setprecision(2);
    for ( auto& f : m_fitFractions ){
        SummaryFile << "$ " << latexName(f.name()) << " $ & $ " << f.val() * 100. << " \\pm " << f.err()* 100. << "$ " << " \\\\ ";
        if(f.name().find( "Sum" ) != std::string::npos) SummaryFile << "\\hline";
        SummaryFile << "\n";
    }
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n \n";

    SummaryFile << "\\begin{tabular}{l r} " << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\multicolumn{1}{c}{Parameter} & \\multicolumn{1}{c}{Value} " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";

    for (size_t i = 0; i < (size_t)m_mps->size(); ++i ) {
        auto param = m_mps->at(i);
        double scale = 1;
        if(!param->isFree() || param->name().find( "::Spline" ) != std::string::npos )continue;
        if(param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos) continue;
        if(param->name().find( "_mass" ) != std::string::npos || param->name().find( "_width" ) != std::string::npos){
            scale = 1000;
            SummaryFile << std::fixed << std::setprecision(2);            
        }
        else{
            scale = 1;
            SummaryFile << std::fixed << std::setprecision(4);            
        }
        SummaryFile << "$ " << latexName(param->name()) << " $ & $ " << param->mean() * scale << " \\pm " << param->err()* scale << "$ " << " \\\\ " << "\n";
    }
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";
}

void FitResult::writeToOptionsFile( const std::string& fname, int fixParams )
{
    std::ofstream outlog;
    outlog << std::setprecision( 4 );
    outlog.open( fname );
    double scale_error = 1.;
    for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
        auto param = m_mps->at(i);
        if(param->isHidden() || param->name().find( "::Spline" ) != std::string::npos )continue;
        //if(!param->isFree())continue;

        
        if(param->name().find( "_Re" ) != std::string::npos){
            int flag = (int)param->flag() ;
            if((param->name().find( "FOCUS" ) != std::string::npos || param->name().find( "kMatrix" ) != std::string::npos) && fixParams>0 ) flag = 2;
                
            outlog << replaceAll(param->name(),"_Re", "") << "  " << flag << " "
                   << param->mean() << " " << ( param->isFree() ? m_mps->at(i)->err()/scale_error : 0 ) << " "
                   << param->minInit() << " ";
            
            if(abs(param->mean()-param->maxInit())<10*param->stepInit()) outlog << param->maxInit()*2 << "    " ;       
            else outlog << param->maxInit() << "    " ;               
            
            auto param_im = m_mps->at(i+1);
            if(param->name().find( replaceAll(param_im->name(),"_Im", "")) == std::string::npos)
                ERROR("Found no matching imaginary part for " << param->name() << "; have " << param_im->name());

            outlog << flag<< " " << param_im->mean() << " " << (param_im->isFree() ? m_mps->at(i+1)->err()/scale_error : 0 ) << " "
            << param_im->minInit() << " " << param_im->maxInit() << " " ;               
            i++;            
        }
        else if(param->name().find( "_Im" ) != std::string::npos){
            auto param_re = m_mps->at(i+1);
            if(param->name().find( replaceAll(param_re->name(),"_Re", "")) == std::string::npos)ERROR("Found no matching real part for " << param->name() << "; have " << param_re->name());
            
            int flag = (int)param->flag() ;
            if((param->name().find( "FOCUS" ) != std::string::npos || param->name().find( "kMatrix" ) != std::string::npos) && fixParams>0 ) flag = 2;
            outlog << replaceAll(param_re->name(),"_Re", "") << "  " << flag << " "
                   << param_re->mean() << " " << ( param_re->isFree() ? m_mps->at(i+1)->err()/scale_error : 0 ) << " "
                   << param_re->minInit() << " ";
            if(abs(param_re->mean()-param_re->maxInit())<10*param_re->stepInit()) outlog << param_re->maxInit()*2 << "    " ;
            else outlog << param_re->maxInit() << "    " ;
            
            outlog << flag << " " << param->mean() << " " << (param->isFree() ? m_mps->at(i)->err()/scale_error : 0 ) << " "
            << param->minInit() << " " << param->maxInit() << " " ;
            
            
            i++;
        }
        else if(param->name().find( "_mass" ) != std::string::npos || param->name().find( "_width" ) != std::string::npos){
            outlog << param->name() << "  " << ( fixParams==0 ? (int)param->flag() : 2) << " "
            << param->mean() << " " << ( param->isFree() ? m_mps->at(i)->stepInit() : 0) << " "
            << param->minInit() << " " << param->maxInit() << "    " ;
        }
        else{
            outlog << param->name() << "  " << (int)param->flag() << " " 
            << param->mean() << " " << ( param->isFree() ? m_mps->at(i)->stepInit() : 0) << " "
            << param->minInit() << " " << param->maxInit() << "    " ;               
        }
        outlog << std::endl;
    }
    outlog.close();
}

void FitResult::writeToRootFile(TFile * output, unsigned seed, int verbose, double nll, unsigned nAmps, double Ns, std::vector<double> thresholds, std::vector<double> numFracAboveThresholds){
        
        double chi2;
        unsigned status, nPar;
        double res[m_fitFractions.size()*2];
        double par[m_mps->size()*2];

        double r[thresholds.size()];
        double AIC[thresholds.size()];
        double BIC[thresholds.size()];
        double AIC2[thresholds.size()];
        double BIC2[thresholds.size()];

        output->cd();
        TTree* outputTree = new TTree("Result","Result");
        outputTree->Branch("status", &status, "status/I");
        outputTree->Branch("nll", &nll, "nll/D");
        outputTree->Branch("chi2", &chi2, "chi2/D");
        outputTree->Branch("seed", &seed, "seed/I");
        outputTree->Branch("nPar", &nPar, "nPar/I");
        outputTree->Branch("nAmps", &nAmps, "nAmps/I");
        outputTree->Branch("nSig", &Ns);

        for (unsigned int i = 0; i < m_fitFractions.size(); i++ ){
            res[i] = m_fitFractions[i].val();
            res[i+m_fitFractions.size()] = m_fitFractions[i].err();
            outputTree->Branch(((std::string)AmpGen::programatic_name(m_fitFractions[i].name())).c_str(), &res[i]);
            if(verbose>=1)outputTree->Branch(((std::string)AmpGen::programatic_name(m_fitFractions[i].name())+"_err").c_str(), &res[i+m_fitFractions.size()]);
        }
    
        for (size_t i = 0; i < (size_t)m_mps->size(); i++ ) {
            auto param = m_mps->at(i);
            if(!param->isFree() || param->name().find( "::Spline" ) != std::string::npos )continue;
            if(verbose<2)if(param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos) continue;
            
            par[i] = param->mean();
            par[i+m_mps->size()] = param->err();
            outputTree->Branch(((std::string)AmpGen::programatic_name(param->name())).c_str(), &par[i]);
            if(verbose>=1)outputTree->Branch(((std::string)AmpGen::programatic_name(param->name())+"_err").c_str(), &par[i+m_fitFractions.size()]);
        }
    
        for (unsigned int i = 0; i < thresholds.size(); i++ ){
            AIC[i] = nll + 2 * numFracAboveThresholds[i];
            BIC[i] = nll + numFracAboveThresholds[i] * log(Ns);
            AIC2[i] = nll + 4 * numFracAboveThresholds[i];
            BIC2[i] = nll + 2 * numFracAboveThresholds[i] * log(Ns);
            r[i] = numFracAboveThresholds[i];
            
            outputTree->Branch(("AIC_" + std::to_string((int)(thresholds[i]*10))).c_str(), &AIC[i]);
            outputTree->Branch(("AIC2_" + std::to_string((int)(thresholds[i]*10))).c_str(), &AIC2[i]);
            outputTree->Branch(("BIC_" + std::to_string((int)(thresholds[i]*10))).c_str(), &BIC[i]);
            outputTree->Branch(("BIC2_" + std::to_string((int)(thresholds[i]*10))).c_str(), &BIC2[i]);
            outputTree->Branch(("r_" + std::to_string((int)(thresholds[i]*10))).c_str(), &r[i]);            
        }    
            
        chi2 = m_chi2 / dof();
        status = m_status;
        nPar = nParam();
    
        outputTree->Fill();
        outputTree->Write();
}

std::complex<double> BW_val(const double& s, const double& m, const double& gamma){
    //std::complex<double> BW = -complex<double>(0,1) * m * gamma/(m*m - s -  complex<double>(0,1) * m * gamma);
    std::complex<double> BW = -std::complex<double>(0,1) * m * gamma/(m*m - s -  std::complex<double>(0,1) * sqrt(s) * gamma);
    return BW;
}

void FitResult::plotSpline( const std::string& name, const std::string& outDir ) {
        
    double min, max, nBins( 0 );
    auto spline_params = NamedParameter<double>( name + "::Spline").getVector();
    if( spline_params.size() == 3 ){
        nBins = int( spline_params[0] );
        min   =         spline_params[1] ; 
        max   =         spline_params[2];
    }
    else if( spline_params.size() == 4 ){
        nBins = size_t( spline_params[0] );
        min   =       spline_params[1] * spline_params[1] ; 
        max   =       spline_params[2] * spline_params[2];
    }
    else {
        nBins = NamedParameter<int>( name + "::Spline::N"  , 0. );
        min   = NamedParameter<double>( name + "::Spline::Min", 0. );
        max   = NamedParameter<double>( name + "::Spline::Max", 0. );
    }
    
    TGraphErrors* g_amp = new TGraphErrors();
    TGraphErrors* g_phase = new TGraphErrors();
    TGraphErrors* g_argand = new TGraphErrors();
    
    TGraphErrors* g_amp_bw = new TGraphErrors();
    TGraphErrors* g_phase_bw = new TGraphErrors();
    TGraphErrors* g_argand_bw = new TGraphErrors();

    double amp,phase,amp_err,phase_err;
    
    std::ofstream outlog;
    outlog << std::setprecision( 4 );
    outlog.open( (outDir+"/"+name+"_spline.txt").c_str() );
    outlog << "Head " + name  << std::endl;
    outlog << name + "::Spline";
    for(auto& p : spline_params)outlog << " " << p;
    outlog << std::endl;
    
    auto mass = m_mps->find(name+"_mass");
    auto gamma = m_mps->find(name+"_width");
    
    outlog << mass->name() << "  " << (int)mass->flag() << " "
    << mass->mean() << " " << mass->stepInit()  << " "
    << mass->minInit() << " " << mass->maxInit() << "    " ;
    outlog << std::endl;
    
    outlog << gamma->name() << "  " << (int)gamma->flag() << " "
    << gamma->mean() << " " << gamma->stepInit()  << " "
    << gamma->minInit() << " " << gamma->maxInit() << "    " ;
    outlog << std::endl;
    
    int nBins_bw = 100;
    double min_bw = min; //pow(m - 2 * gamma,2);
    double max_bw = max;//pow(m + 2 * gamma,2);
    double st_bw = (max_bw-min_bw)/double(nBins_bw-1);
    for ( int c = 0; c < nBins_bw; ++c ) {
            double s = min + double( c ) * st_bw;
        
            std::complex<double> BW = BW_val(s,mass->mean(),gamma->mean());

            amp = std::abs(BW);
            amp_err = 0;
            
            phase = std::arg(BW)*180./3.141;
            phase_err = 0;
            
            g_amp_bw->SetPoint( c, sqrt(s), amp );
            g_phase_bw->SetPoint( c, sqrt(s), phase );
        
            g_amp_bw->SetPointError( c, 0, amp_err );
            g_phase_bw->SetPointError( c, 0, phase_err );
        
            g_argand_bw->SetPoint( c, amp * cos(phase/180.*M_PI), amp * sin(phase/180.*M_PI)  );
            g_argand_bw->SetPointError( c, sqrt( pow(amp_err * cos(phase/180.*M_PI),2) + pow(amp*sin(phase/180.*M_PI)*phase_err/180.*M_PI, 2) ) , sqrt( pow(amp_err * sin(phase/180.*M_PI),2) + pow(amp*cos(phase/180.*M_PI)*phase_err/180.*M_PI, 2) )   );
    }

    double st = (max-min)/double(nBins-1);

    int c_norm = NamedParameter<int>( "Spline::NormBinBW"  , nBins/2 - 1 );
    double s_norm = min + ( c_norm ) * st;

    std::complex<double> BW_norm = BW_val(s_norm,mass->mean(),gamma->mean());
    double amp_bw_norm = std::abs(BW_norm);
    double phase_bw_norm = std::arg(BW_norm)*180./3.141;

    double amp_norm, phase_norm;
    for (size_t i = 0; i < (size_t)m_mps->size(); ++i ) {
        auto param = m_mps->at(i);
        if(param->name().find( name + "::Spline::") != std::string::npos){
            if(param->name() == name + "::Spline::Re::" + std::to_string(c_norm)){
                amp_norm = param->mean() ;
            }
            else if(param->name() == name + "::Spline::Im::" + std::to_string(c_norm)){
                phase_norm = param->mean() ;
            }
        }
    }

    for ( int c = 0; c < nBins; ++c ) {
            double s = min + double( c ) * st;
            
            for (size_t i = 0; i < (size_t)m_mps->size(); ++i ) {
                auto param = m_mps->at(i);
                if(param->name().find( name + "::Spline::") != std::string::npos){      
                    if(param->name() == name + "::Spline::Re::" + std::to_string(c)){
                        amp = param->mean() / amp_norm * amp_bw_norm;
                        amp_err = param->err() / amp_norm * amp_bw_norm;
                        outlog << param->name() << "  " << (int)param->flag() << " "
                        << param->mean() << " " << param->err()  << " "
                        << param->minInit() << " " << param->maxInit() << "    " ;
                        outlog << std::endl;
                    }
                    else if(param->name() == name + "::Spline::Im::" + std::to_string(c)){
                        phase = param->mean() - phase_norm + phase_bw_norm;
                        phase_err = param->err();
                        outlog << param->name() << "  " << (int)param->flag() << " "
                        << param->mean() << " " << param->err()  << " "
                        << param->minInit() << " " << param->maxInit() << "    " ;
                        outlog << std::endl;
                    }
                }
            }
            g_amp->SetPoint( c, sqrt(s), amp );
            g_phase->SetPoint( c, sqrt(s), phase );
        
            g_amp->SetPointError( c, 0, amp_err );
            g_phase->SetPointError( c, 0, phase_err );
        
            g_argand->SetPoint( c, amp * cos(phase/180.*M_PI), amp * sin(phase/180.*M_PI)  );
            g_argand->SetPointError( c, sqrt( pow(amp_err * cos(phase/180.*M_PI),2) + pow(amp*sin(phase/180.*M_PI)*phase_err/180.*M_PI, 2) ) , sqrt( pow(amp_err * sin(phase/180.*M_PI),2) + pow(amp*cos(phase/180.*M_PI)*phase_err/180.*M_PI, 2) )   );
    }
    
    outlog.close();
    
    TCanvas* c = new TCanvas();
    g_amp->SetMarkerColor(4);
    g_amp->SetMarkerStyle(20);
    g_amp->Draw("APC");
    c->Print((outDir+"/"+name+"_amp.pdf").c_str());
    g_phase->SetMarkerColor(4);
    g_phase->SetMarkerStyle(20);
    g_phase->Draw("APC");
    c->Print((outDir+"/"+name+"_phase.pdf").c_str());
    g_argand->SetMarkerColor(4);
    g_argand->SetMarkerStyle(20);
    g_argand->Draw("APC");
    c->Print((outDir+"/"+name+"_argand.pdf").c_str());
    
    g_amp->SetMarkerColor(4);
    g_amp->SetMarkerStyle(20);
    g_amp->Draw("APL");
    g_amp_bw->SetLineColor(kRed);
    g_amp_bw->SetLineWidth(3);
    g_amp_bw->Draw("CSAME");
    c->Print((outDir+"/"+name+"_amp_bw.pdf").c_str());
    
    g_phase->SetMarkerColor(4);
    g_phase->SetMarkerStyle(20);
    g_phase->Draw("APL");
    g_phase_bw->SetLineColor(kRed);
    g_phase_bw->SetLineWidth(3);
    g_phase_bw->Draw("CSAME");
    c->Print((outDir+"/"+name+"_phase_bw.pdf").c_str());
    
    g_argand->SetMarkerColor(4);
    g_argand->SetMarkerStyle(20);
    g_argand->Draw("APL");
    g_argand_bw->SetLineColor(kRed);
    g_argand_bw->SetLineWidth(3);
    g_argand_bw->Draw("CSAME");
    c->Print((outDir+"/"+name+"_argand_bw.pdf").c_str());

}


void FitResult::print() const
{
  INFO( "Chi2 per bin = " << m_chi2 / m_nBins );
  INFO( "Chi2 per dof = " << m_chi2 / dof() );
  INFO( "-2LL         = " << m_LL );
  INFO( "Fit Status   = " << m_status );
}

std::vector<MinuitParameter*> FitResult::parameters() const
{
  std::vector<MinuitParameter*> params( m_mps->size() );
  std::copy( m_mps->begin(), m_mps->end(), params.begin() );
  return params;
}

std::vector<MinuitParameter*> FitResult::floating( const bool& extended ) const
{
  std::vector<MinuitParameter*> floating;
  for ( auto& param : *m_mps ) {
    if ( ( param->isFree() || extended ) && param->err() > 1e-6 ) floating.push_back( param );
  }
  return floating;
}

TMatrixD FitResult::getReducedCovariance( const bool& extended ) const
{
  std::vector<unsigned int> floating_indices;
  for (size_t i = 0; i < m_mps->size(); ++i ) {
    auto param = m_mps->at(i);
    if ( ( param->isFree() || extended ) && param->err() > 1e-6 ) floating_indices.push_back( i );
  }
  TMatrixD reducedCov( floating_indices.size(), floating_indices.size() );
  if ( int( floating_indices.size() ) > m_covarianceMatrix.GetNrows() ) {
    ERROR( "Looking for more floating indices than real indices: " << m_covarianceMatrix.GetNrows() << " "
        << floating_indices.size() );
    return reducedCov;
  }
  for ( unsigned int i = 0; i < floating_indices.size(); ++i ) {
    for ( unsigned int j = 0; j < floating_indices.size(); ++j ) {
      reducedCov( i, j ) = m_covarianceMatrix( floating_indices[i], floating_indices[j] );
    }
  }
  return reducedCov;
}

LinearErrorPropagator FitResult::getErrorPropagator( const bool& extended ) const
{
  return LinearErrorPropagator( getReducedCovariance( extended ), floating( extended ) );
}

void FitResult::addChi2( const double& chi2, const double& nBins )
{
  m_chi2  = chi2;
  m_nBins = nBins;
}

void FitResult::addFractions( const std::vector<FitFraction>& fractions ) 
{ 
  m_fitFractions = fractions; 
}

double FitResult::chi2() const { return m_chi2; }
double FitResult::LL() const { return m_LL; }
int FitResult::status() const { return m_status; }
int FitResult::nParam() const { return m_nParam; }
int FitResult::nBins() const { return m_nBins; }

std::map<std::string, double> FitResult::observables() const { return m_observables; }
MinuitParameterSet* FitResult::mps() const { return m_mps; }

double FitResult::dof() const { return m_nBins - m_nParam - 1; }
std::vector<FitFraction> FitResult::fitFractions() const { return m_fitFractions; }
TMatrixD FitResult::cov() const { return m_covarianceMatrix; }
double FitResult::cov( const size_t& x, const size_t& y ) const { return m_covarianceMatrix( x, y ); }
double FitResult::cov( const std::string& x, const std::string& y ) const 
{
  auto tx = m_covMapping.find(x);
  auto ty = m_covMapping.find(y);
  if( tx == m_covMapping.end() || ty == m_covMapping.end() ){
    ERROR("Parameter not found: " << x << ", " << y );
    for (const auto& [key, value] : m_covMapping)
          std::cout << '[' << key << "] = " << value << "; ";
    return -1;
  }
  return cov(tx->second, ty->second);
}

void FitResult::addFraction( const std::string& name, const double& frac, const double& err )
{
  m_fitFractions.emplace_back( name, frac, err );
}
void FitResult::clearFitFractions() { m_fitFractions.clear(); }
void FitResult::setCov( const size_t& x, const size_t& y, const double& F ) { m_covarianceMatrix( x, y ) = F; }
double FitResult::correlation( const std::string& x, const std::string& y ) const
{
  auto tx = m_covMapping.find(x);
  auto ty = m_covMapping.find(y);
  if( tx == m_covMapping.end() || ty == m_covMapping.end() ){
    ERROR("Parameter not found: " << x << ", " << y );
    return -1;
  }
  return cov(tx->second, ty->second)/sqrt(cov(tx->second, tx->second)*cov(ty->second, ty->second));
}
