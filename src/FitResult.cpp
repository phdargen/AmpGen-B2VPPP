#include <TMatrixDfwd.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixTUtils.h>
#include <TGraph.h>
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

bool FitResult::readFile( const std::string& fname )
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
  }
  size_t nParameters = parameterLines.size();
  m_covarianceMatrix.ResizeTo( parameterLines.size(), parameterLines.size() );
  for (size_t i = 0; i < nParameters; ++i ) {
    auto tokens             = split( parameterLines[i], ' ' );
    m_covMapping[tokens[1]] = i;
    INFO(tokens[1]);
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
    outlog << "End Log\n";
    outlog.close();
}

std::string FitResult::latexName(std::string name){
    
    TString n(name);

    n.ReplaceAll("{","\\{");
    n.ReplaceAll("}","\\}");
    
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
    n.ReplaceAll("[NONRESEXP]","");
    n.ReplaceAll("[Flatte]","");

    n.ReplaceAll("GounarisSakurai.Omega","");
    n.ReplaceAll("GounarisSakurai","");
    n.ReplaceAll("LASS","");
    n.ReplaceAll("GLASS","");
    n.ReplaceAll("Bugg","");
    n.ReplaceAll("Kappa","");
    n.ReplaceAll("Dabba","");
    n.ReplaceAll("GSpline.EFF","");
    n.ReplaceAll("SBW","");
    n.ReplaceAll("NONRESEXP","");

    n.ReplaceAll("(0)","_{0}");
    n.ReplaceAll("(1)","_{1}");
    n.ReplaceAll("(2)","_{2}");

    n.ReplaceAll("(S)","_{S}");
    n.ReplaceAll("(P)","_{P}");
    n.ReplaceAll("(V)","_{V}");
    n.ReplaceAll("(A)","_{A}");
    n.ReplaceAll("(T)","_{T}");
    n.ReplaceAll("(PT)","_{PT}");

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

void FitResult::writeToOptionsFile( const std::string& fname )
{
    std::ofstream outlog;
    outlog << std::setprecision( 4 );
    outlog.open( fname );
    for (size_t i = 0; i < (size_t)m_covarianceMatrix.GetNrows(); ++i ) {
        auto param = m_mps->at(i);
        if(param->isHidden() || param->name().find( "::Spline" ) != std::string::npos )continue;
        if(param->name().find( "_Re" ) != std::string::npos){
            outlog << replaceAll(param->name(),"_Re", "") << "  " << (int)param->flag() << " " 
                   << param->mean() << " " << ( param->isFree() ? m_mps->at(i)->err()/5 : 0 ) << " "
                   << param->minInit() << " ";
            if(abs(param->mean()-param->maxInit())<10*param->stepInit()) outlog << param->maxInit()*2 << "    " ;       
            else outlog << param->maxInit() << "    " ;               
            auto param_im = m_mps->at(i+1);
            outlog << (int)param_im->flag() << " " << param_im->mean() << " " << (param_im->isFree() ? m_mps->at(i)->err()/5 : 0 ) << " "
            << param_im->minInit() << " " << param_im->maxInit() << " " ;               
            i++;            
        }
        else{
            outlog << param->name() << "  " << (int)param->flag() << " " 
            << param->mean() << " " << ( param->isFree() ? m_mps->at(i)->err()/5 : 0 ) << " "
            << param->minInit() << " " << param->maxInit() << "    " ;               
        }
        outlog << std::endl;
    }
    outlog.close();
}

void FitResult::writeToRootFile(TFile * output, unsigned seed, int verbose, unsigned nAmps, double Ns, std::vector<double> thresholds, std::vector<double> numFracAboveThresholds){
        
        double nll,chi2;
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
    
        for (size_t i = 0; i < (size_t)m_mps->size(); ++i ) {
            auto param = m_mps->at(i);
            if(!param->isFree() || param->name().find( "::Spline" ) != std::string::npos )continue;
            if(verbose<2)if(param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos) continue;
            
            par[i] = param->mean();
            par[i+m_mps->size()] = param->err();
            outputTree->Branch(((std::string)AmpGen::programatic_name(param->name())).c_str(), &par[i]);
            if(verbose>=1)outputTree->Branch(((std::string)AmpGen::programatic_name(param->name())+"_err").c_str(), &par[i+m_fitFractions.size()]);
        }
    
        for (unsigned int i = 0; i < thresholds.size(); i++ ){
            AIC[i] = m_LL + 2 * numFracAboveThresholds[i];
            BIC[i] = m_LL + numFracAboveThresholds[i] * log(Ns);
            AIC2[i] = m_LL + 4 * numFracAboveThresholds[i];
            BIC2[i] = m_LL + 2 * numFracAboveThresholds[i] * log(Ns);
            r[i] = numFracAboveThresholds[i];
            
            outputTree->Branch(("AIC_" + std::to_string((int)(thresholds[i]*10))).c_str(), &AIC[i]);
            outputTree->Branch(("AIC2_" + std::to_string((int)(thresholds[i]*10))).c_str(), &AIC2[i]);
            outputTree->Branch(("BIC_" + std::to_string((int)(thresholds[i]*10))).c_str(), &BIC[i]);
            outputTree->Branch(("BIC2_" + std::to_string((int)(thresholds[i]*10))).c_str(), &BIC2[i]);
            outputTree->Branch(("r_" + std::to_string((int)(thresholds[i]*10))).c_str(), &r[i]);            
        }    
            
        nll = m_LL;
        chi2 = m_chi2 / dof();
        status = m_status;
        nPar = nParam();
    
        outputTree->Fill();
        outputTree->Write();
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
    
    TGraph* g_amp = new TGraph();
    TGraph* g_phase = new TGraph();
    TGraph* g_argand = new TGraph();

    double st = (max-min)/double(nBins-1);
    double amp,phase;

    for ( int c = 0; c < nBins; ++c ) {
            double s = min + double( c ) * st;
            
            for (size_t i = 0; i < (size_t)m_mps->size(); ++i ) {
                auto param = m_mps->at(i);
                if(param->name().find( name + "::Spline::") != std::string::npos){      
                            if(param->name().find( "Re::" + std::to_string(c) ) != std::string::npos) amp = param->mean();
                            else if(param->name().find( "Im::" + std::to_string(c) ) != std::string::npos) phase = param->mean();
                }
            }
            g_amp->SetPoint( g_amp->GetN(), sqrt(s), amp );
            g_phase->SetPoint( g_phase->GetN(), sqrt(s), phase );
            g_argand->SetPoint( g_argand->GetN(), amp, phase );        
    }
    
    TCanvas* c = new TCanvas();
    g_amp->Draw("AC*");
    c->Print((outDir+"/"+name+"_amp.eps").c_str());
    g_phase->Draw("AC*");
    c->Print((outDir+"/"+name+"_phase.eps").c_str());
    g_argand->Draw("AC*");
    c->Print((outDir+"/"+name+"_argand.eps").c_str());
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
double FitResult::cov( const std::string& x, const std::string& y ) const { return m_covarianceMatrix( m_covMapping.at(x), m_covMapping.at(y) ); }

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
