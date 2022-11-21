#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/LHCbStyle.h"
#include <TH1.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/FitResult.h"
#include "TMatrixD.h"

using namespace std;
using namespace AmpGen;


void regulateParameters( FitResult& fr ){  
  auto params = fr.parameters();
  std::vector< std::pair<size_t,size_t > > amplitude_parameters; 
  std::vector<double> scale_factors( params.size(), 1 );
  for( size_t i = 0 ; i < params.size(); ++i ){
    if( ! params[i]->isFree() ) continue;
    std::string name = params[i]->name();
    if( name.find("_Re") == std::string::npos ) continue;  
    auto im_name = replaceAll(name,"_Re","_Im");
    size_t j = 0 ; 
    for( ; j < params.size() ; ++j){
      if( params[j]->name() == im_name ) break;
    }
    auto re=params[i];
    auto im=params[j];
    double mean_re = re->mean(); 
    double mean_im = im->mean(); 
    if(  mean_re < 0 ){
      mean_re = -mean_re;
      mean_im = mean_im + M_PI;
      scale_factors[ i ] = -1 ;
    }
    mean_im = atan2( sin(mean_im), cos(mean_im) ); 
    re->setCurrentFitVal( mean_re );
    im->setCurrentFitVal( mean_im );
  }  
  for( size_t x = 0 ; x < params.size(); ++x){
    for( size_t y = 0 ; y < params.size() ;++y){
      double sf = scale_factors[x] * scale_factors[y];
      fr.setCov(x,y, fr.cov(x,y) *sf );
    } 
  }
}

double form_ring( const double& x, const double& x0) {
  double minimum_distance=(x-x0)*(x-x0);
  double distance = minimum_distance ;
  double xc=x;
  bool higherOrLower=false;
  do {
    minimum_distance = distance ;
    double minus_one =  (xc-2*M_PI -x0)*(xc-2*M_PI-x0);
    double plus_one =   (xc+2*M_PI -x0)*(xc+2*M_PI-x0);
    higherOrLower = minus_one > plus_one;
    distance =  higherOrLower ? plus_one : minus_one;
    xc        = higherOrLower ? xc + 2*M_PI : xc - 2 * M_PI;

  } while( distance < minimum_distance ) ;
  xc = higherOrLower ? xc - 2*M_PI : xc + 2*M_PI;

  return xc;
}

void EnsureRing( const FitResult& fit, const FitResult& starterFit ){
  auto params = starterFit.floating();
  for( auto param : params ){
    if( param->name().find("_Im") == std::string::npos ) continue; 
    MinuitParameter* other_param = fit.mps()->find(param->name() ) ; 
    if( other_param == nullptr ){
         ERROR(param->name() << " not found in fit! [map size = " << fit.mps()->size() << "]");
         return;
    }
    other_param->setCurrentFitVal( form_ring( other_param->mean(), param->mean() ) );
  }
}

template< typename InputIterator, typename ReturnType = double , typename InputType = typename std::remove_pointer<decltype(std::declval<InputIterator>().p)>::type > 
    std::pair<ReturnType,ReturnType> auto_stat( InputIterator first, InputIterator last,  const std::function<ReturnType(const InputType&)>& func, const std::function<double(const InputType&)>& wf = nullptr ){
    ReturnType x=0;
    ReturnType xx=0;
    double N=0;
    while( first != last ){
        ReturnType f = func(*first);
        double weight = wf == nullptr ? 1 : wf(*first);
        x+=f * weight;
        xx+=f*f * weight ;
        N+=weight;
        ++first;
    }
    x/=N;
    xx = xx/N -x*x;
    return std::pair<ReturnType,ReturnType>(x,sqrt(xx));
}

TMatrixD diffMethod( const std::vector<FitResult>& fits, const FitResult& starterFit){
  auto params = starterFit.floating();
  TMatrixD covMatrix( params.size(), params.size() );
  
  //std::vector<std::pair<double,double>> mean_test;
  for( unsigned int i = 0 ; i < params.size(); ++i ){
    std::string name_i = params[i]->name();
      
    std::function<double(const FitResult& f)> getMean = [&](const FitResult& f){ return f.mps()->find(name_i)->mean() ; };
    std::function<double(const FitResult& f)> getErr = [&](const FitResult& f){ return f.mps()->find(name_i)->err() ; };
    auto meanAndSigma = auto_stat( fits.begin(), fits.end(), getMean ); 
      
    covMatrix(i,i) = pow( meanAndSigma.first - params[i]->mean() ,2 );
    if( params[i]->isFree() ){
      INFO( params[i]->name() << "   " << params[i]->mean() << " +/- " << params[i]->err() << "; Mean = "  << meanAndSigma.first  <<  " sigma_sys/sigma_stat = " << sqrt( covMatrix(i,i) ) / params[i]->err()  );
    }
    //mean_test.push_back( meanAndSigma );
  }
  return covMatrix; 
}

TMatrixD covMethod( const std::vector<FitResult>& fits, const FitResult& starterFit ){
      
  auto params = starterFit.floating();
  TMatrixD SIGMA( params.size(), params.size() );
  TVectorD param_vector( params.size() );
    
  for( size_t iFit = 0; iFit<fits.size();++iFit ){
    auto& fit = fits[iFit];
    fit.print();  
    TMatrixD rc( params.size(), params.size() );
    for( unsigned int i = 0 ; i < params.size(); ++i ){
      for( unsigned int j = 0 ; j < params.size(); ++j ){
          rc(i,j) =  fit.cov( params[i]->name(), params[j]->name() );
      }
    }
    auto irc = rc.Invert();
    SIGMA += irc; 
    TVectorD xvec( params.size() );
    for( unsigned int i = 0 ; i < params.size(); ++i){
      xvec[i] = fit.mps()->find( params[i]->name() )->mean();
    }
    param_vector += irc * xvec ;
  }
      
  auto M_matrix = SIGMA.Invert();
  TVectorD combined_vector = M_matrix * param_vector;
 
  //auto allParameters = starterFit.parameters();
  auto allParameters = starterFit.floating();
  auto size = allParameters.size();
  
  std::vector<int> floatingToFull( params.size() , -1 );
  std::vector<double> diagonalErrors( params.size() ,0);
  for( size_t i = 0 ; i < params.size(); ++i ){
    double param  = params[i]->mean();
    double effSys = combined_vector[i];
    diagonalErrors[i] = (param-effSys)*(param-effSys);
    INFO( params[i]->name() << " " << param << " " << effSys << " " << params[i]->err() );
    for( size_t j = 0 ; j < allParameters.size(); ++j ){
      if( allParameters[j]->name() == params[i]->name() ){
        floatingToFull[i] = j;
        break;
      }
    }
    if( floatingToFull[i] == -1 ){
      ERROR("Parameter: " << params[i]->name() << " not found in full covariance matrix");
    }
    auto paramName = params[i]->name();
    //auto hist = auto_histogram( fits, [&paramName](auto& fit){ return fit.parameter(paramName )->mean() ; } , paramName+"_regulated" ,50);
    //hist->Write();
    //delete hist;
  }
  auto reducedCovariance = starterFit.getReducedCovariance();
  TMatrixD covMatrix( size, size );
      
  for( size_t i = 0 ; i < params.size(); ++i){
    for( size_t j = 0 ; j < params.size(); ++j){
      covMatrix(floatingToFull[i],floatingToFull[j]) =  reducedCovariance(i,j) * sqrt( diagonalErrors[i] * diagonalErrors[j] ) / sqrt( reducedCovariance(i,i) * reducedCovariance(j,j) ); 
    }
    auto param = params[i];
    double cov = sqrt( covMatrix( floatingToFull[i] , floatingToFull[i] ) );
    INFO( std::setw(60) << std::left << param->name() << " " << std::setprecision(4) 
        << std::setw(10) <<  std::right << round(  param->mean(),4) << " ± " 
        << std::setw(8) <<  round( param->err() ,4) << " ± " 
        << std::setw(8) <<  round(  cov ,4) << " ( " 
        << std::setw(8) << round( cov /  param->err() ,4) << " %)" );
  }
  return covMatrix; 
}

TMatrixD combMethod( const std::vector<FitResult>& fits ){
  const FitResult& f0 = fits[0];
  auto mps = f0.floating();
  std::vector<double> means( mps.size(), 0 ); 
  std::vector<std::string> names (mps.size(), "");
  INFO("nParameters = " << mps.size() );
    
  std::vector< std::vector<MinuitParameter* > > params(  mps.size(), std::vector<MinuitParameter*>(fits.size(),0) );
  for( size_t x = 0 ; x < mps.size(); ++x){
    names[x] = mps[x]->name();
    for( size_t iFit = 0; iFit < fits.size() ;++iFit ){
      params[x][iFit] = fits[iFit].mps()->find( names[x] );
      means[x] += params[x][iFit]->mean();
    }
    means[x] /= double( fits.size() );
  }
  INFO("Calculated means, calculating covariance matrix");
  TMatrixD covMatrix( means.size(), means.size() );
  for( unsigned int i = 0 ; i < means.size(); ++i){ 
    auto   name_i = names[i];
    double val_i  = means[i];
    for( size_t j = 0 ; j < means.size(); ++j ){
      auto   name_j = names[j]; 
      double val_j  = means[j];
      for( size_t iFit = 0 ; iFit < fits.size(); ++iFit ){
        covMatrix(i,j) += 
          ( params[i][iFit]->mean() - val_i )*( params[j][iFit]->mean() - val_j );
      }
    }
  }
  for( unsigned int i = 0 ; i < means.size(); ++i ){
    for( unsigned int j = 0 ; j < means.size(); ++j ){
      covMatrix(i,j) = covMatrix(i,j) / double( fits.size() );
    }
  }
  return covMatrix; 
}


/*
void doAnalysis(const std::vector<Fit>& fits, const Fit& starterFit ){
  INFO("Doing systematic analysis");
  Fit& mutable_fit = const_cast<Fit&>(starterFit) ; /// evil hack /// 
  auto params = starterFit.fitResult().floating();
  regulateParameters( mutable_fit.fitResultMutable() );
  for( auto& fit : fits ){
    Fit& this_mutable_fit = const_cast<Fit&>(fit) ; /// evil hack /// 
    regulateParameters( this_mutable_fit.fitResultMutable()  );
    EnsureRing( fit, starterFit );   
  }
  std::string mode = NamedParameter<std::string>("Mode","Combine"); // sigma or diff ///
  INFO("Mode = " << mode );
  const FitResult& fr = mode=="Combine"? fits[0].fitResult() : starterFit.fitResult();
  size_t     size     = fr.mps()->size();
  TMatrixD covMatrix( size, size );
  if( mode == "Combine" ){
    INFO("Using mode combine");
    covMatrix = combMethod( fits, starterFit );
  }
  else if( mode == "Diff"    ) covMatrix = diffMethod( fits, starterFit );
  else if( mode == "Cov"     ){
    covMatrix = CovMethod(  fits, starterFit );
  }
  else {
    ERROR("Systematic derivation : " << mode << " not recognised"); 
  }
  AmpGen::MinuitParameterSet output_params; 
  for( unsigned int i = 0 ; i < size ; ++i){
    auto iparam = fr.parameters()[i];
    std::string paramName = iparam->name();
    std::function<double(const Fit& fit)> val = [&paramName](const Fit& fit){ return fit.parameter(paramName)->mean(); };
    std::function<double(const Fit& fit)> err = [&paramName](const Fit& fit){ return fit.parameter(paramName)->err(); };
    std::pair<double,double> mu = auto_stat( fits.begin(), fits.end(),val);
    std::pair<double,double> sig = auto_stat( fits.begin(), fits.end(),err);
    double error = sqrt( covMatrix(i,i) );
    AmpGen::MinuitParameter* param = 
      new MinuitParameter(paramName, iparam->flag(), mu.first, error , iparam->minInit(), iparam->maxInit() );
    output_params.add( param);
    
    INFO( std::setw(60) << std::left << param->name() << " " << std::setprecision(4) 
        << std::setw(10) <<  std::right << round(  param->mean(),4) << " ± " 
        << std::setw(8) <<  round( sig.first ,4) << " ± " 
        << std::setw(8) <<  round(  sqrt( covMatrix(i,i)  ),4) << " ( " 
        << std::setw(8) << round( sqrt(covMatrix(i,i) )  / sig.first,4) << " %)" );
    
    //auto hist = auto_histogram( fits, [&paramName](auto& fit){ return fit.parameter(paramName )->mean() ; } , paramName ,50);
    //hist->Write();
    //delete hist;
  }
  AmpGen::FitResult output( output_params, covMatrix );
  if( NamedParameter<bool>("CalculateFractions",0)==1 ){
    for( auto& ff : starterFit.processes() ){
      auto name = ff.first; 
      std::function<double(const Fit& fit)> frac = [&name](const Fit& fit){ return fit.fraction(name)->val(); };
      std::pair<double,double> mus = auto_stat( fits.begin(), fits.end(),frac);
      output.addFraction( name, mus.first, mus.second ); 
    }
  }
  std::function<double(const Fit&)> LL_func = [](auto& fit){ return fit.likelihood() ; };
  auto    chi2 = [](auto& fit){ return fit.fitResult().chi2(); };
  //double avg_ll = auto_stat( fits.begin(), fits.end(), LL_func ).first;
  //auto dLL_func = [&avg_ll](auto& fit){ return fit.likelihood() - avg_ll; };
  
//  auto_histogram( fits, LL_func, "delta_LL",50 )->Write();
//  auto_histogram( fits, chi2, "chi2",50)->Write();
//  std::string outputFile = AmpGen::NamedParameter<std::string>("Project","").getVal()+"/results.dat";
  output.writeToFile(outputFile);
}
*/


void analyzeResults(){

  INFO("Doing systematic analysis");    
  std::string baseLineFit = NamedParameter<std::string>("baseLineFit", "log.txt");  
    
  std::string inDir = NamedParameter<std::string>("inDir", "sys/");
  auto sysFiles = NamedParameter<std::string>("sysFiles", std::vector<std::string>()).getVector();

  FitResult starterFit(baseLineFit);  
  starterFit.print();
    
  vector<FitResult> fits;  
    
  for(auto&f : sysFiles){
      FitResult fr(inDir+f);
      EnsureRing( fr, starterFit );    
      fits.push_back(fr);
  }

  TMatrixD cov = diffMethod( fits, starterFit);
  cov.Print();
    
  cov = combMethod( fits);
  cov.Print();
    
  TMatrixD cov3 = covMethod( fits, starterFit);
  cov3.Print();

  //starterFit.setCov((int)starterFit.parameters().size(),(int)starterFit.parameters().size() ,cov);
  //TMatrixD cov_reduced =  starterFit.getReducedCovariance();
  //cov_reduced.Print();
    
}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  analyzeResults();

  return 0;
}
