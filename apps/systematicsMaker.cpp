#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
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
  vector< pair<size_t,size_t > > amplitude_parameters; 
  vector<double> scale_factors( params.size(), 1 );
  for( size_t i = 0 ; i < params.size(); ++i ){
    if( ! params[i]->isFree() ) continue;
    string name = params[i]->name();
    if( name.find("_Re") == string::npos ) continue;  
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
    if( param->name().find("_Im") == string::npos ) continue; 
    MinuitParameter* other_param = fit.mps()->find(param->name() ) ; 
    if( other_param == nullptr ){
         //ERROR(param->name() << " not found in fit! [map size = " << fit.mps()->size() << "]");
         return;
    }
    other_param->setCurrentFitVal( form_ring( other_param->mean(), param->mean() ) );
  }
}

template< typename InputIterator, typename ReturnType = double , typename InputType = typename remove_pointer<decltype(declval<InputIterator>().p)>::type > 
    pair<ReturnType,ReturnType> auto_stat( InputIterator first, InputIterator last,  const function<ReturnType(const InputType&)>& func, const function<double(const InputType&)>& wf = nullptr ){
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
    return pair<ReturnType,ReturnType>(x,sqrt(xx));
}

void sanityChecks(MinuitParameterSet& mps){
   for(int i=0;i<mps.size();i++){
       if((mps[i]->name().find( "_mass" ) != string::npos || mps[i]->name().find( "_width" ) != string::npos || mps[i]->name().find( "_alpha" ) != string::npos || mps[i]->name().find( "_beta" ) != string::npos ) ){

               INFO("Found " << mps[i]->name() << ", remove it from MinuitParameterSet");
               mps.unregister(mps[i]);  
               i=0;                
           }
    }
}  



TMatrixD diffMethod( const vector<FitResult>& fits, const FitResult& starterFit){
  auto params = starterFit.floating();
  TMatrixD covMatrix( params.size(), params.size() );
  
  //vector<pair<double,double>> mean_test;
  for( unsigned int i = 0 ; i < params.size(); ++i ){
    string name_i = params[i]->name();
      
    function<double(const FitResult& f)> getMean = [&](const FitResult& f){ return f.mps()->find(name_i)->mean() ; };
    function<double(const FitResult& f)> getErr = [&](const FitResult& f){ return f.mps()->find(name_i)->err() ; };
    auto meanAndSigma = auto_stat( fits.begin(), fits.end(), getMean ); 
      
    covMatrix(i,i) = pow( meanAndSigma.first - params[i]->mean() ,2 );
    if( params[i]->isFree() ){
      INFO( params[i]->name() << "   " << params[i]->mean() << " +/- " << params[i]->err() << "; Mean = "  << meanAndSigma.first  <<  " sigma_sys/sigma_stat = " << sqrt( covMatrix(i,i) ) / params[i]->err()  );
    }
    //mean_test.push_back( meanAndSigma );
  }
  return covMatrix; 
}


vector<TMatrixD> sampleVarMethod( const vector<FitResult>& fits, const FitResult& starterFit, const string& mode = "Bias"){

    // Fit parameters    
    auto params = starterFit.floating();
    TMatrixD covMatrix( params.size(), params.size() );

    for( unsigned int i = 0 ; i < params.size(); ++i ){
      string name_i = params[i]->name();
      double mean = 0;
      double maxDiff = 0;  
      for(auto& fit : fits){
          auto params_fit = fit.floating();
          for(auto& p: params_fit){
              if(name_i == p->name() ){
                  mean += p->mean();    
                  maxDiff = abs(p->mean() - params[i]->mean()) > maxDiff ? abs(p->mean() - params[i]->mean()) : maxDiff;
              }
          }
      }
      if(mean==0){
          ERROR("Parameter " << name_i << " not found ");
          covMatrix(i,i) = 0;
          continue;
      }  
        
      mean/=double(fits.size());
      
      double var = 0;
      for(auto& fit : fits){
          auto params_fit = fit.floating();
          for(auto& p: params_fit){
              if(name_i == p->name() ) var += pow(p->mean() - mean,2);            
          }
      }        
      if(mode=="Var")var/=double(fits.size()-1.);
  
      if(mode=="Bias")covMatrix(i,i) = pow(mean - params[i]->mean(),2 );
      if(mode=="Var") covMatrix(i,i) = var;
      if(mode=="MaxDiff") covMatrix(i,i) = pow(maxDiff,2);
    }  

    // Fit fractions
    auto fracs = starterFit.fitFractions();
    TMatrixD covMatrix_frac( fracs.size(), fracs.size() );

    for( unsigned int i = 0 ; i < fracs.size(); ++i ){
        string name_i = fracs[i].name();          
        double mean = 0;
        double maxDiff = 0;  
        for(auto& fit : fits){
            auto fracs_fit = fit.fitFractions();
            for(auto& f: fracs_fit){
                if(name_i == f.name() ){
                    mean += f.val();     
                    maxDiff = abs(f.val() - fracs[i].val()) > maxDiff ? abs(f.val() - fracs[i].val()) : maxDiff;
                }
            }
        }
        if(mean==0){
            ERROR("Fraction " << name_i << " not found ");
            covMatrix_frac(i,i) = 0;
            continue;
        }  
        mean/=double(fits.size());
        
        double var = 0;
        for(auto& fit : fits){
            auto fracs_fit = fit.fitFractions();
            for(auto& f: fracs_fit){
                if(name_i == f.name() )var += pow(f.val() - mean,2);           
            }
        }
        if(mode=="Var")var/=double(fits.size()-1.);

        if(mode=="Bias")covMatrix_frac(i,i) = pow(mean - fracs[i].val(),2);
        if(mode=="Var") covMatrix_frac(i,i) = var;
        if(mode=="MaxDiff") covMatrix_frac(i,i) = pow(maxDiff,2);
    }
    
    return vector<TMatrixD>({covMatrix,covMatrix_frac}); 
}

TMatrixD covMethod( const vector<FitResult>& fits, const FitResult& starterFit ){
      
  auto params = starterFit.floating();
  TMatrixD SIGMA( params.size(), params.size() );
  TVectorD param_vector( params.size() );
    
  for( size_t iFit = 0; iFit<fits.size();++iFit ){
    auto& fit = fits[iFit];
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
  
  vector<int> floatingToFull( params.size() , -1 );
  vector<double> diagonalErrors( params.size() ,0);
  for( size_t i = 0 ; i < params.size(); ++i ){
    double param  = params[i]->mean();
    double effSys = combined_vector[i];
    diagonalErrors[i] = (param-effSys)*(param-effSys);
    //INFO( params[i]->name() << " " << param << " " << effSys << " " << params[i]->err() );
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
    /*  
    INFO( setw(60) << left << param->name() << " " << setprecision(4) 
        << setw(10) <<  right << round(  param->mean(),4) << " ± " 
        << setw(8) <<  round( param->err() ,4) << " ± " 
        << setw(8) <<  round(  cov ,4) << " ( " 
        << setw(8) << round( cov /  param->err() ,4) << " %)" );
     */
  }
  return covMatrix; 
}

TMatrixD combMethod( const vector<FitResult>& fits ){
  const FitResult& f0 = fits[0];
  sanityChecks(*f0.mps());  
  auto mps = f0.floating();

  vector<double> means( mps.size(), 0 ); 
  vector<string> names (mps.size(), "");
  INFO("nParameters = " << mps.size() );
    
  vector< vector<MinuitParameter* > > params(  mps.size(), vector<MinuitParameter*>(fits.size(),0) );
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
void doAnalysis(const vector<Fit>& fits, const Fit& starterFit ){
  INFO("Doing systematic analysis");
  Fit& mutable_fit = const_cast<Fit&>(starterFit) ; /// evil hack /// 
  auto params = starterFit.fitResult().floating();
  regulateParameters( mutable_fit.fitResultMutable() );
  for( auto& fit : fits ){
    Fit& this_mutable_fit = const_cast<Fit&>(fit) ; /// evil hack /// 
    regulateParameters( this_mutable_fit.fitResultMutable()  );
    EnsureRing( fit, starterFit );   
  }
  string mode = NamedParameter<string>("Mode","Combine"); // sigma or diff ///
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
    string paramName = iparam->name();
    function<double(const Fit& fit)> val = [&paramName](const Fit& fit){ return fit.parameter(paramName)->mean(); };
    function<double(const Fit& fit)> err = [&paramName](const Fit& fit){ return fit.parameter(paramName)->err(); };
    pair<double,double> mu = auto_stat( fits.begin(), fits.end(),val);
    pair<double,double> sig = auto_stat( fits.begin(), fits.end(),err);
    double error = sqrt( covMatrix(i,i) );
    AmpGen::MinuitParameter* param = 
      new MinuitParameter(paramName, iparam->flag(), mu.first, error , iparam->minInit(), iparam->maxInit() );
    output_params.add( param);
    
    INFO( setw(60) << left << param->name() << " " << setprecision(4) 
        << setw(10) <<  right << round(  param->mean(),4) << " ± " 
        << setw(8) <<  round( sig.first ,4) << " ± " 
        << setw(8) <<  round(  sqrt( covMatrix(i,i)  ),4) << " ( " 
        << setw(8) << round( sqrt(covMatrix(i,i) )  / sig.first,4) << " %)" );
    
    //auto hist = auto_histogram( fits, [&paramName](auto& fit){ return fit.parameter(paramName )->mean() ; } , paramName ,50);
    //hist->Write();
    //delete hist;
  }
  AmpGen::FitResult output( output_params, covMatrix );
  if( NamedParameter<bool>("CalculateFractions",0)==1 ){
    for( auto& ff : starterFit.processes() ){
      auto name = ff.first; 
      function<double(const Fit& fit)> frac = [&name](const Fit& fit){ return fit.fraction(name)->val(); };
      pair<double,double> mus = auto_stat( fits.begin(), fits.end(),frac);
      output.addFraction( name, mus.first, mus.second ); 
    }
  }
  function<double(const Fit&)> LL_func = [](auto& fit){ return fit.likelihood() ; };
  auto    chi2 = [](auto& fit){ return fit.fitResult().chi2(); };
  //double avg_ll = auto_stat( fits.begin(), fits.end(), LL_func ).first;
  //auto dLL_func = [&avg_ll](auto& fit){ return fit.likelihood() - avg_ll; };
  
//  auto_histogram( fits, LL_func, "delta_LL",50 )->Write();
//  auto_histogram( fits, chi2, "chi2",50)->Write();
//  string outputFile = AmpGen::NamedParameter<string>("Project","").getVal()+"/results.dat";
  output.writeToFile(outputFile);
}
*/


void analyzeResults(){

    INFO("Doing systematic analysis");    
    string baseLineFit = NamedParameter<string>("baseLineFit", "log.txt");  
    string outDir = NamedParameter<string>("outDir", "sys/out/");  

    FitResult starterFit(baseLineFit);  
    starterFit.print();
    auto params = starterFit.floating();
    auto fracs = starterFit.fitFractions();

    vector<TMatrixD> covs;
    vector<TMatrixD> covs_frac;
        
    TMatrixD covTot(params.size(),params.size());
    TMatrixD covTot_frac(fracs.size(),fracs.size());

    auto includeSys = NamedParameter<string>("includeSys", vector<string>()).getVector();
    vector<string> sysNames;

    for(auto& sys : includeSys){

        auto sysName = NamedParameter<string>("sysName"+sys, sys);
        sysNames.push_back(sysName);
        auto sysFiles = NamedParameter<string>("sysFiles"+sys, vector<string>()).getVector(); // basename, first file number, last file number
        auto sysMethod = NamedParameter<string>("sysMethod"+sys, "Var");

        vector<FitResult> fits;  
        for(unsigned int i = stoi(sysFiles[1]) ; i <= stoi(sysFiles[2]) ; ++i){
          if(! std::ifstream((sysFiles[0] + to_string(i)+".txt").c_str()).good()){
              ERROR("File " << sysFiles[0] + to_string(i)+".txt" << " not found");
              continue;
          }
          FitResult fr(sysFiles[0] + to_string(i)+".txt");
          EnsureRing( fr, starterFit );    
          fits.push_back(fr);
        }

        vector<TMatrixD> cov({TMatrixD(params.size(),params.size()),TMatrixD(fracs.size(),fracs.size())});
        if(sysFiles.size()>0){
            cov = sampleVarMethod(fits,starterFit,sysMethod);
            cout <<  sys << " systematic with method " << sysMethod << endl;
            cout <<  "Fit parameters " << endl;    
            for( unsigned int i = 0 ; i < params.size() ; ++i){
                  INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
                      << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
                      << setw(8) <<  round( params[i]->err() ,4) << " ± " 
                      << setw(8) <<  round(  sqrt( cov[0](i,i)  ),4) << " ( " 
                      << setw(8) << round( sqrt(cov[0](i,i) )  / params[i]->err() * 100.,4) << " %)" );  
            }
            cout <<  "Fit fractions " << endl;    
            for( unsigned int i = 0 ; i < fracs.size() ; ++i){
                  if(fracs[i].err()<0.0001)continue;
                  INFO( setw(60) << left << fracs[i].name() << " " << setprecision(4) 
                      << setw(10) <<  right << round(  fracs[i].val() * 100,4) << " ± " 
                      << setw(8) <<  round( fracs[i].err() * 100,4) << " ± " 
                      << setw(8) <<  round(  sqrt( cov[1](i,i)  ) * 100,4) << " ( " 
                      << setw(8) << round( sqrt(cov[1](i,i) )  / fracs[i].err() * 100.,4) << " %)" );  
            }
            covTot += cov[0];
            covTot_frac += cov[1];        
            covs.push_back(cov[0]);
            covs_frac.push_back(cov[1]);
        }   
    
    }
        
    /*
    // Bkg  
    auto sysFilesBkg = NamedParameter<string>("sysFilesBkg", vector<string>()).getVector(); // basename, number of files
    auto sysMethodBkg = NamedParameter<string>("sysMethodBkg", "Diff");

    vector<FitResult> fitsBkg;  
    for(unsigned int i = 1 ; i <= stoi(sysFilesBkg[1]) ; ++i){
        FitResult fr(sysFilesBkg[0] + to_string(i)+".txt");
        EnsureRing( fr, starterFit );    
        fitsBkg.push_back(fr);
    }

    TMatrixD covBkg(params.size(),params.size());
    if(sysMethodBkg=="Diff"){
        covBkg = diffMethod( fitsBkg, starterFit);
    }
    else if(sysMethodBkg=="Comb"){
        covBkg = combMethod( fitsBkg);
    }
    else if(sysMethodBkg=="Cov"){
        covBkg = covMethod( fitsBkg, starterFit);
    }

    // Print output  
    INFO(endl << "Bkg systematic");
    for( unsigned int i = 0 ; i < params.size() ; ++i){
        
        INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
            << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
            << setw(8) <<  round( params[i]->err() ,4) << " ± " 
            << setw(8) <<  round(  sqrt( covBkg(i,i)  ),4) << " ( " 
            << setw(8) << round( sqrt(covBkg(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
        
    }  
    
    // FF  
    auto sysFilesFF = NamedParameter<string>("sysFilesFF", vector<string>()).getVector(); // basename, number of files
    auto sysMethodFF = NamedParameter<string>("sysMethodFF", "Diff");

    vector<FitResult> fitsFF;  
    for(unsigned int i = 1 ; i <= stoi(sysFilesFF[1]) ; ++i){
        FitResult fr(sysFilesFF[0] + to_string(i)+".txt");
        EnsureRing( fr, starterFit );    
        fitsFF.push_back(fr);
    }

    TMatrixD covFF(params.size(),params.size());
    if(sysMethodFF=="Diff"){
        covFF = diffMethod( fitsFF, starterFit);
    }
    else if(sysMethodFF=="Comb"){
        covFF = combMethod( fitsFF);
    }
    else if(sysMethodFF=="Cov"){
        covFF = covMethod( fitsFF, starterFit);
    }

    // Print output  
    INFO(endl << "FF systematic");
    for( unsigned int i = 0 ; i < params.size() ; ++i){
        
        INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
            << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
            << setw(8) <<  round( params[i]->err() ,4) << " ± " 
            << setw(8) <<  round(  sqrt( covFF(i,i)  ),4) << " ( " 
            << setw(8) << round( sqrt(covFF(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
        
    }  
    
    // NR  
    auto sysFilesNR = NamedParameter<string>("sysFilesNR", vector<string>()).getVector(); // basename, number of files
    auto sysMethodNR = NamedParameter<string>("sysMethodNR", "Diff");

    if(sysFilesNR.size()>0){

        vector<FitResult> fitsNR;  
        for(unsigned int i = 1 ; i <= stoi(sysFilesNR[1]) ; ++i){
            FitResult fr(sysFilesNR[0] + to_string(i)+".txt");
            EnsureRing( fr, starterFit );    
            fitsNR.push_back(fr);
        }

        TMatrixD covNR(params.size(),params.size());
        if(sysMethodNR=="Diff"){
            covNR = diffMethod( fitsNR, starterFit);
        }
        else if(sysMethodNR=="Comb"){
            covNR = combMethod( fitsNR);
        }
        else if(sysMethodNR=="Cov"){
            covNR = covMethod( fitsNR, starterFit);
        }

        // Print output
        INFO(endl << "NR systematic");
        for( unsigned int i = 0 ; i < params.size() ; ++i){
            INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
                << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
                << setw(8) <<  round( params[i]->err() ,4) << " ± " 
                << setw(8) <<  round(  sqrt( covNR(i,i)  ),4) << " ( " 
                << setw(8) << round( sqrt(covNR(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
            
        }  
        covTot += covNR; 
    }
    
    // Resonances  
    auto sysFilesRes = NamedParameter<string>("sysFilesRes", vector<string>()).getVector(); // basename, number of files
    auto sysMethodRes = NamedParameter<string>("sysMethodRes", "Diff");

    vector<FitResult> fitsRes;  
    for(unsigned int i = 0 ; i < stoi(sysFilesRes[1]) ; ++i){
        FitResult fr(sysFilesRes[0] + to_string(i)+".txt");
        EnsureRing( fr, starterFit );    
        fitsRes.push_back(fr);
    }

    TMatrixD covRes(params.size(),params.size());
    if(sysMethodRes=="Diff"){
        covRes = diffMethod( fitsRes, starterFit);
    }
    else if(sysMethodRes=="Comb"){
        covRes = combMethod( fitsRes);
    }
    else if(sysMethodRes=="Cov"){
        covRes = covMethod( fitsRes, starterFit);
    }

    // Print output
    INFO(endl << "Res systematic");
    for( unsigned int i = 0 ; i < params.size() ; ++i){
        
        INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
            << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
            << setw(8) <<  round( params[i]->err() ,4) << " ± " 
            << setw(8) <<  round(  sqrt( covRes(i,i)  ),4) << " ( " 
            << setw(8) << round( sqrt(covRes(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
        
    }  
    
    // LS  
    auto sysFilesLS = NamedParameter<string>("sysFilesLS", vector<string>()).getVector(); // basename, number of files
    auto sysMethodLS = NamedParameter<string>("sysMethodLS", "Diff");

    vector<FitResult> fitsLS;  
    for(unsigned int i = 1 ; i <= stoi(sysFilesLS[1]) ; ++i){
        FitResult fr(sysFilesLS[0] + to_string(i)+".txt");
        EnsureRing( fr, starterFit );    
        fitsLS.push_back(fr);
    }

    TMatrixD covLS(params.size(),params.size());
    if(sysMethodLS=="Diff"){
        covLS = diffMethod( fitsLS, starterFit);
    }
    else if(sysMethodLS=="Comb"){
        covLS = combMethod( fitsLS);
    }
    else if(sysMethodLS=="Cov"){
        covLS = covMethod( fitsLS, starterFit);
    }

    // Print output
    INFO(endl << "LS systematic");
    for( unsigned int i = 0 ; i < params.size() ; ++i){
        
        INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
            << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
            << setw(8) <<  round( params[i]->err() ,4) << " ± " 
            << setw(8) <<  round(  sqrt( covLS(i,i)  ),4) << " ( " 
            << setw(8) << round( sqrt(covLS(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
        
    }  
    
    // Alt models  
    auto sysFilesAltAmp = NamedParameter<string>("sysFilesAltAmp", vector<string>()).getVector(); // basename, number of files
    auto sysMethodAltAmp = NamedParameter<string>("sysMethodAltAmp", "Diff");
    TMatrixD covAltAmp(params.size(),params.size());

    if(sysFilesAltAmp.size()>0){
        vector<FitResult> fitsAltAmp;  
        for(unsigned int i = 0 ; i < stoi(sysFilesAltAmp[1]) ; ++i){
            FitResult fr(sysFilesAltAmp[0] + to_string(i)+".txt");
            EnsureRing( fr, starterFit );    
            fitsAltAmp.push_back(fr);
        }

        if(sysMethodAltAmp=="Diff"){
            covAltAmp = diffMethod( fitsAltAmp, starterFit);
        }
        else if(sysMethodAltAmp=="Comb"){
            covAltAmp = combMethod( fitsAltAmp);
        }
        else if(sysMethodAltAmp=="Cov"){
            covAltAmp = covMethod( fitsAltAmp, starterFit);
        }

        // Print output
        INFO(endl << "AltAmp systematic");
        for( unsigned int i = 0 ; i < params.size() ; ++i){
            INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
                << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
                << setw(8) <<  round( params[i]->err() ,4) << " ± " 
                << setw(8) <<  round(  sqrt( covAltAmp(i,i)  ),4) << " ( " 
                << setw(8) << round( sqrt(covAltAmp(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
            
        }  
        covTot += covAltAmp;
    }
    
    */
    
    // Total
    cout << "Total systematic" << endl;
    //covTot.Print();
    
    cout <<  "Fit parameters " << endl;    
    for( unsigned int i = 0 ; i < params.size() ; ++i){
            INFO( setw(60) << left << params[i]->name() << " " << setprecision(4) 
                << setw(10) <<  right << round(  params[i]->mean(),4) << " ± " 
                << setw(8) <<  round( params[i]->err() ,4) << " ± " 
                << setw(8) <<  round(  sqrt( covTot(i,i)  ),4) << " ( " 
                << setw(8) << round( sqrt(covTot(i,i) )  / params[i]->err() * 100.,4) << " %)" );  
    }    
    cout <<  "Fit fractions " << endl;    
    for( unsigned int i = 0 ; i < fracs.size() ; ++i){
          if(fracs[i].err()<0.0001)continue;
          INFO( setw(60) << left << fracs[i].name() << " " << setprecision(4) 
              << setw(10) <<  right << round(  fracs[i].val() * 100,4) << " ± " 
              << setw(8) <<  round( fracs[i].err() * 100,4) << " ± " 
              << setw(8) <<  round(  sqrt( covTot_frac(i,i)  ) * 100,4) << " ( " 
              << setw(8) << round( sqrt(covTot_frac(i,i) )  / fracs[i].err() * 100.,4) << " %)" );  
    }
    
    AmpGen::MinuitParameterSet output_params; 
    for( unsigned int i = 0 ; i < params.size() ; ++i){
      AmpGen::MinuitParameter* param = new MinuitParameter(params[i]->name(), params[i]->flag(), params[i]->mean(), sqrt( covTot(i,i)  ) , params[i]->minInit(), params[i]->maxInit() );
      output_params.add( param);      
    }
    AmpGen::FitResult output( output_params, covTot );
    output.writeToFile(outDir+"log_sys.txt");

    // Create tables

    // Total systematics table in terms of sigma_stat   
    ofstream SummaryFile,SummaryFile2;
    SummaryFile.open(outDir+"sys_summary_table.tex",ofstream::trunc);
    SummaryFile2.open(outDir+"sys_summary_table2.tex",ofstream::trunc);

    SummaryFile << "\\begin{tabular}{l " ;
    for(int i =0 ; i < covs.size() ; i++) SummaryFile << " c " ;
    SummaryFile << " | c }" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Fit Parameter & " ;
    for(int i =0 ; i < covs.size() ; i++)  SummaryFile << sysNames[i] << " & " ;
    SummaryFile << " Total " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";
    
    SummaryFile2 << "\\begin{tabular}{l " ;
    for(int i =0 ; i < covs.size() ; i++) SummaryFile2 << " c " ;
    SummaryFile2 << " | c }" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "Fit Parameter & " ;
    for(int i =0 ; i < covs.size() ; i++)  SummaryFile2 << sysNames[i] << " & " ;
    SummaryFile2 << " Total " << " \\\\ " << "\n";
    SummaryFile2 << "\\hline" << "\n";

    for(int i =0 ; i < params.size() ; i++){
        if( params[i]->name().find( "B+" ) != string::npos ){
            double tot = 0.;
            SummaryFile << fixed << setprecision(2) << "$" << starterFit.latexName(params[i]->name()) << "$ & " ;
            for(int j =0 ; j <covs.size() ; j++){
                tot += (covs[j])[i][i];
                SummaryFile << sqrt((covs[j])[i][i])/params[i]->err() << " & ";  
            }
            SummaryFile << sqrt(tot)/params[i]->err() << " \\\\ " << "\n"; 
        }
        else {
            double tot = 0.;
            SummaryFile2 << fixed << setprecision(2) << "$" << starterFit.latexName(params[i]->name()) << "$ & " ;
            for(int j =0 ; j <covs.size() ; j++){
                tot += (covs[j])[i][i];
                SummaryFile2 << sqrt((covs[j])[i][i])/params[i]->err() << " & ";  
            }
            SummaryFile2 << sqrt(tot)/params[i]->err() << " \\\\ " << "\n"; 
        }
    }

    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";

    
    // Total fractions systematics table in terms of sigma_stat   
    ofstream FracSummaryFile,FracSummaryFile2;
    FracSummaryFile.open(outDir+"sys_fracSummary_table.tex",ofstream::trunc);
    FracSummaryFile2.open(outDir+"sys_fracSummary_table2.tex",ofstream::trunc);
    
    FracSummaryFile << "\\begin{tabular}{l " ;
    for(int i =0 ; i < covs_frac.size() ; i++) FracSummaryFile << " c " ;
    FracSummaryFile << " | c }" << "\n";
    FracSummaryFile << "\\hline" << "\n";
    FracSummaryFile << "\\hline" << "\n";
    FracSummaryFile << "Fit fraction & " ;
    for(int i =0 ; i < covs_frac.size() ; i++)  FracSummaryFile << sysNames[i] << " & " ;
    FracSummaryFile << " Total " << " \\\\ " << "\n";
    FracSummaryFile << "\\hline" << "\n";
    
    FracSummaryFile2 << "\\begin{tabular}{l " ;
    for(int i =0 ; i < covs_frac.size() ; i++) FracSummaryFile2 << " c " ;
    FracSummaryFile2 << " | c }" << "\n";
    FracSummaryFile2 << "\\hline" << "\n";
    FracSummaryFile2 << "\\hline" << "\n";
    FracSummaryFile2 << "Fit fraction & " ;
    for(int i =0 ; i < covs_frac.size() ; i++)  FracSummaryFile2 << sysNames[i] << " & " ;
    FracSummaryFile2 << " Total " << " \\\\ " << "\n";
    FracSummaryFile2 << "\\hline" << "\n";
    
    ofstream FracResultFile,FracResultFile2;
    FracResultFile.open(outDir+"sys_fracResult_table.tex",ofstream::trunc);
    FracResultFile2.open(outDir+"sys_fracResult_table2.tex",ofstream::trunc);

    FracResultFile << "\\begin{tabular}{l c} " << "\n";
    FracResultFile << "\\hline" << "\n";
    FracResultFile << "\\hline" << "\n";
    FracResultFile << "Decay channel &  Fit fraction " << " \\\\ " << "\n";
    FracResultFile << "\\hline" << "\n";

    for(int i =0 ; i < fracs.size() ; i++){
        if( fracs[i].err()<0.0001)continue;
        if( fracs[i].name().find( "B+" ) != string::npos ){
            double tot = 0.;
            FracSummaryFile << fixed << setprecision(2) << "$" << starterFit.latexName(fracs[i].name()) << "$ & " ;
            for(int j =0 ; j <covs_frac.size() ; j++){
                tot += (covs_frac[j])[i][i];
                FracSummaryFile << sqrt((covs_frac[j])[i][i])/fracs[i].err() << " & ";  
            }
            FracSummaryFile << sqrt(tot)/fracs[i].err() << " \\\\ " << "\n"; 
            
            FracResultFile << fixed << setprecision(2) << "$" << starterFit.latexName(fracs[i].name()) << "$ & $" ;
            FracResultFile << fracs[i].val() * 100. << " \\pm " ;
            FracResultFile << fracs[i].err() * 100. << " \\pm " ;
            if(tot>0)FracResultFile << sqrt(tot) * 100. ;
            FracResultFile << "$ \\\\ " << "\n";  
        }
        else {
            double tot = 0.;
            FracSummaryFile2 << fixed << setprecision(2) << "$" << starterFit.latexName(fracs[i].name()) << "$ & " ;
            for(int j =0 ; j <covs_frac.size() ; j++){
                tot += (covs_frac[j])[i][i];
                FracSummaryFile2 << sqrt((covs_frac[j])[i][i])/fracs[i].err() << " & ";  
            }
            FracSummaryFile2 << sqrt(tot)/fracs[i].err() << " \\\\ " << "\n"; 
            
            FracResultFile2 << fixed << setprecision(2) << "$" << starterFit.latexName(fracs[i].name()) << "$ & $" ;
            FracResultFile2 << fracs[i].val() * 100. << " \\pm " ;
            FracResultFile2 << fracs[i].err() * 100. << " \\pm " ;
            if(tot>0)FracResultFile2 << sqrt(tot) * 100. ;
            FracResultFile2 << "$ \\\\ " << "\n"; 
        }
    }

    FracSummaryFile << "\\hline" << "\n";
    FracSummaryFile << "\\hline" << "\n";
    FracSummaryFile << "\\end{tabular}" << "\n";
    
    FracSummaryFile2 << "\\hline" << "\n";
    FracSummaryFile2 << "\\hline" << "\n";
    FracSummaryFile2 << "\\end{tabular}" << "\n";
    
    FracResultFile << "\\hline" << "\n";
    FracResultFile << "\\hline" << "\n";
    FracResultFile << "\\end{tabular}" << "\n";
    
    FracResultFile2 << "\\hline" << "\n";
    FracResultFile2 << "\\hline" << "\n";
    FracResultFile2 << "\\end{tabular}" << "\n";
}


//starterFit.setCov((int)starterFit.parameters().size(),(int)starterFit.parameters().size() ,cov);
//TMatrixD cov_reduced =  starterFit.getReducedCovariance();
//cov_reduced.Print();

int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  analyzeResults();

  return 0;
}
