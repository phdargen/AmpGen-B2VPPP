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


vector<TMatrixD> sampleVarMethod( const vector<FitResult>& fits, const FitResult& starterFit, const FitResult& toyFit, const string& mode = "Bias"){

    // Fit parameters    
    auto params = starterFit.floating();
    TMatrixD covMatrix( params.size(), params.size() );
    
    auto params_toy = toyFit.mps();

    for( unsigned int i = 0 ; i < params.size(); ++i ){
      string name_i = params[i]->name();
      double mean = 0;
      double maxDiff = 0;
      double pull_mean = 0;
      for(auto& fit : fits){
          auto params_fit = fit.floating();
          for(auto& p: params_fit){
              if(name_i == p->name() || name_i == replaceAll(p->name(),".BL","") ){
                  mean += p->mean();
                  maxDiff = abs(p->mean() - params[i]->mean()) > maxDiff ? abs(p->mean() - params[i]->mean()) : maxDiff;
                  if(params_toy->find(name_i))pull_mean += (params_toy->find(name_i)->mean()-p->mean())/p->err();
                  else ERROR(name_i << " not found in toy log file");
              }
          }
      }
      if(mean==0){
          ERROR("Parameter " << name_i << " not found ");
          covMatrix(i,i) = 0;
          continue;
      }  
        
      mean/=double(fits.size());
      pull_mean/=double(fits.size());

      double var = 0;
      double pull_var = 0;
      for(auto& fit : fits){
          auto params_fit = fit.floating();
          for(auto& p: params_fit){
              if(name_i == p->name() || name_i == replaceAll(p->name(),".BL","")  ){
                  var += pow(p->mean() - mean,2);
                  if(params_toy->find(name_i))pull_var += pow((params_toy->find(name_i)->mean()-p->mean())/p->err() - pull_mean,2);
                  else ERROR(name_i << " not found in toy log file");
              }
          }
      }        
      var/=double(fits.size()-1.);
      pull_var/=double(fits.size()-1.);
        
      if(mode=="Bias")covMatrix(i,i) = pow(mean - params[i]->mean(),2 );
      if(mode=="Var") covMatrix(i,i) = var;
      if(mode=="MaxDiff") covMatrix(i,i) = pow(maxDiff,2);
        
      if(mode=="PullBias")covMatrix(i,i) = pow(pull_mean * params[i]->err(),2);
      if(mode=="PullBias")INFO(name_i << ": pull = " << pull_mean << " +/- " << sqrt(pull_var) );
    }

    // Fit fractions
    auto fracs = starterFit.fitFractions();
    TMatrixD covMatrix_frac( fracs.size(), fracs.size() );

    auto fracsToy = toyFit.fitFractions();
    
    for( unsigned int i = 0 ; i < fracs.size(); ++i ){
        string name_i = fracs[i].name();          
        double mean = 0;
        double maxDiff = 0;
        double pull_mean = 0;
        for(auto& fit : fits){
            auto fracs_fit = fit.fitFractions();
            for(auto& f: fracs_fit){
                if(name_i == f.name() || name_i == replaceAll(f.name(),".BL","")  ){
                    mean += f.val();
                    maxDiff = abs(f.val() - fracs[i].val()) > maxDiff ? abs(f.val() - fracs[i].val()) : maxDiff;
                    double frac_toy_val = -1;
                    for(auto& f_toy: fracsToy){
                        if(name_i == f_toy.name() ) frac_toy_val = f_toy.val();
                        else if(name_i == replaceAll(f_toy.name(),".BL","") ) frac_toy_val = f_toy.val();
                    }
                    if(frac_toy_val<0)ERROR(name_i << " not found in log file");
                    pull_mean += (f.val() - frac_toy_val)/f.err();
                }
            }
        }
        if(mean==0){
            ERROR("Fraction " << name_i << " not found ");
            covMatrix_frac(i,i) = 0;
            continue;
        }  
        mean/=double(fits.size());
        pull_mean/=double(fits.size());

        double var = 0;
        double pull_var = 0;
        for(auto& fit : fits){
            auto fracs_fit = fit.fitFractions();
            for(auto& f: fracs_fit){
                if(name_i == f.name() ){
                    var += pow(f.val() - mean,2);
                    double frac_toy_val = -1;
                    for(auto& f_toy: fracsToy){
                        if(name_i == f_toy.name() || name_i == replaceAll(f.name(),".BL","")  ) frac_toy_val = f_toy.val();
                        else if(name_i == replaceAll(f_toy.name(),".BL","")  ) frac_toy_val = f_toy.val();
                    }
                    if(frac_toy_val<0)ERROR(name_i << " not found in log file");
                    pull_var += pow((f.val() - frac_toy_val)/f.err() - pull_mean,2);
                }
            }
        }
        var/=double(fits.size()-1.);
        pull_var/=double(fits.size()-1.);
        
        if(mode=="Bias")covMatrix_frac(i,i) = pow(mean - fracs[i].val(),2);
        if(mode=="Var") covMatrix_frac(i,i) = var;
        if(mode=="MaxDiff") covMatrix_frac(i,i) = pow(maxDiff,2);
        
        if(mode=="PullBias")covMatrix_frac(i,i) = pow(pull_mean * fracs[i].err(),2);
        if(mode=="PullBias")INFO(name_i << ": pull = " << pull_mean << " +/- " << sqrt(pull_var) );
    }
    
    return vector<TMatrixD>({covMatrix,covMatrix_frac}); 
}

vector<TMatrixD> diffMethod(const FitResult& fit1, const FitResult& fit2, const FitResult& starterFit){

    // Fit parameters    
    auto params = starterFit.floating();
    TMatrixD covMatrix( params.size(), params.size() );

    for( unsigned int i = 0 ; i < params.size(); ++i ){
      string name_i = params[i]->name();

      double val1 = 0;
      auto params_fit1 = fit1.floating();
      for(auto& p: params_fit1) 
          if(name_i == p->name() ) val1 = p->mean();    

      double val2 = 0;
      auto params_fit2 = fit2.floating();
      for(auto& p: params_fit2) 
          if(name_i == p->name() ) val2 = p->mean();    
      
      if(val1==0 || val2==0){
          ERROR("Parameter " << name_i << " not found ");
          covMatrix(i,i) = 0;
          continue;
      }          
      //INFO(name_i << " val1 = " << val1 << " val2 = " << val2 << " diff " <<  (val1-val2)/2. );  
      covMatrix(i,i) = pow((val1-val2)/2.,2);
    }  

    // Fit fractions
    auto fracs = starterFit.fitFractions();
    TMatrixD covMatrix_frac( fracs.size(), fracs.size() );

    for( unsigned int i = 0 ; i < fracs.size(); ++i ){
        string name_i = fracs[i].name();          

        double val1 = 0;
        auto fracs_fit1 = fit1.fitFractions();
        for(auto& f: fracs_fit1) 
            if(name_i == f.name() ) val1 = f.val();    

        double val2 = 0;
        auto fracs_fit2 = fit2.fitFractions();
        for(auto& f: fracs_fit2) 
            if(name_i == f.name() ) val2 = f.val();    
        
        if(val1==0 || val2==0){
            ERROR("Fit fraction " << name_i << " not found ");
            covMatrix_frac(i,i) = 0;
            continue;
        }          
        covMatrix_frac(i,i) = pow((val1-val2)/2.,2);
    }
    
    return vector<TMatrixD>({covMatrix,covMatrix_frac}); 
}


string removeLineshapeMods(string name){
    name =  replaceAll(name,"[GSpline]","");
    name =  replaceAll(name,"[GSpline.BL]","");
    name =  replaceAll(name,"[D;GSpline.BL]","");
    name =  replaceAll(name,"[D;GSpline]","");
    name =  replaceAll(name,"[GounarisSakurai.Omega.BL]","");
    name =  replaceAll(name,"[GounarisSakurai.Omega]","");

    return name;
}

void analyzeResults(){

    INFO("Doing systematic analysis");    
    string baseLineFit = NamedParameter<string>("baseLineFit", "log.txt");
    string outDir = NamedParameter<string>("outDir", "sys/out/");

    FitResult starterFit(baseLineFit);  
    starterFit.print();
    auto params = starterFit.floating();
    auto fracs = starterFit.fitFractions();
    auto interferenceFracs = starterFit.interferenceFractions();

    vector<FitResult> allModelFits;
    allModelFits.push_back(starterFit);
    
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
          if(sys=="AltAmp")allModelFits.push_back(fr);
        }

        vector<TMatrixD> cov({TMatrixD(params.size(),params.size()),TMatrixD(fracs.size(),fracs.size())});
        if(sysFiles.size()>0){
            cout << "Calculating " <<  sys << " systematic with method " << sysMethod << endl;
            
            if(sys=="Res" && sysMethod=="Diff"){
                TMatrixD covMatrix( params.size(), params.size() );
                TMatrixD covMatrix_frac( fracs.size(), fracs.size() );
                for(unsigned int i=0; i < fits.size()/2 ; i++ ){
                    auto tmp = diffMethod(fits[i],fits[i+fits.size()/2],starterFit);
                    covMatrix += tmp[0];
                    covMatrix_frac += tmp[1];   
                    /*
                    for( unsigned int j = 0 ; j < params.size() ; ++j)
                        if(j==6){
                            cout << i << " " << i+fits.size()/2 << endl;
                            cout << sqrt(tmp[0][j][j])/params[j]->err() << endl;
                        }
                     */
                }
                cov = vector<TMatrixD>({covMatrix,covMatrix_frac});
            }
            
            else if(sys=="Bkg" && sysMethod=="Diff"){
                TMatrixD covMatrix( params.size(), params.size() );
                TMatrixD covMatrix_frac( fracs.size(), fracs.size() );
                
                auto tmp = diffMethod(fits[0],fits[1],starterFit);
                covMatrix += tmp[0];
                covMatrix_frac += tmp[1];   
                
                fits.erase(fits.begin());
                fits.erase(fits.begin());                
                tmp = sampleVarMethod(fits,starterFit,starterFit,"Bias");
                covMatrix += tmp[0];
                covMatrix_frac += tmp[1];   

                cov = vector<TMatrixD>({covMatrix,covMatrix_frac});
            }
                
            else if(sys=="Toys"){
                string baseLineToy = NamedParameter<string>("baseLineToy", "");
                FitResult starterToy(baseLineToy);
                cov = sampleVarMethod(fits,starterFit,starterToy,sysMethod);
            }
            
            else cov = sampleVarMethod(fits,starterFit,starterFit,sysMethod);
            
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
    
    // Fit parameter results table
    ofstream ResultFile,ResultFile2;
    ResultFile.open(outDir+"result_table.tex",ofstream::trunc);
    ResultFile2.open(outDir+"result_table2.tex",ofstream::trunc);

    ResultFile << "\\begin{tabular}{l r r} " << "\n";
    ResultFile << "\\hline" << "\n";
    ResultFile << "\\hline" << "\n";
    ResultFile << "Amplitude coupling & Amp & Phase  " << " \\\\ " << "\n";
    ResultFile << "\\hline" << "\n";
    
    ResultFile2 << "\\begin{tabular}{l r} " << "\n";
    ResultFile2 << "\\hline" << "\n";
    ResultFile2 << "\\hline" << "\n";
    ResultFile2 << "Fit Parameter &  Value " << " \\\\ " << "\n";
    ResultFile2 << "\\hline" << "\n";

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
      
    for(int i =0 ; i < params.size() ; i++){

        if( params[i]->name().find( "_Re" ) != string::npos || params[i]->name().find( "_Im" ) != string::npos ){
            double tot_1 = 0.;
            double tot_2 = 0.;
            for(int j =0 ; j <covs.size() ; j++)tot_1 += (covs[j])[i][i];
 
            if(params[i]->name().find( "_Im" )!= string::npos){
                ResultFile << fixed << setprecision(1) << "$" << replaceAll(starterFit.latexName(params[i]->name()),"Phase","") << "$ & $-$ & $" ;
                ResultFile << params[i]->mean()  << " \\pm " ;
                ResultFile << params[i]->err() ;
                if(tot_1>0)ResultFile << " \\pm " << sqrt(tot_1)  ;
            }
            
            else if(params[i]->name().find( "_Re" )!= string::npos){
                ResultFile << fixed << setprecision(2) << "$" << replaceAll(starterFit.latexName(params[i]->name()),"Amp","") << "$ & $" ;
                ResultFile << params[i]->mean()  << " \\pm " ;
                ResultFile << params[i]->err() ;
                if(tot_1>0)ResultFile << " \\pm " << sqrt(tot_1)  ;
                ResultFile << " $ & $ ";
                if(params[i+1]->name().find( "_Im" )!= string::npos){
                    for(int j =0 ; j <covs.size() ; j++)tot_2 += (covs[j])[i+1][i+1];

                    ResultFile << fixed << setprecision(1);
                    ResultFile << params[i+1]->mean()  << " \\pm " ;
                    ResultFile << params[i+1]->err() ;
                    if(tot_2>0)ResultFile << " \\pm " << sqrt(tot_2)  ;
                    i++;
                }
               
            }
            
            ResultFile << "$ \\\\ " << "\n";
        }
        
        else{
            double tot = 0.;
            for(int j =0 ; j <covs.size() ; j++)tot += (covs[j])[i][i];
            ResultFile2 << fixed << setprecision(2) << "$" << starterFit.latexName(params[i]->name()) << "$ & $" ;
            ResultFile2 << params[i]->mean()  << " \\pm " ;
            ResultFile2 << params[i]->err() ;
            if(tot>0)ResultFile2 << " \\pm " << sqrt(tot)  ;
            ResultFile2 << "$ \\\\ " << "\n";
        }
        
    }

    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";
    
    ResultFile << "\\hline" << "\n";
    ResultFile << "\\hline" << "\n";
    ResultFile << "\\end{tabular}" << "\n";
    
    ResultFile2 << "\\hline" << "\n";
    ResultFile2 << "\\hline" << "\n";
    ResultFile2 << "\\end{tabular}" << "\n";

    
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
    
    FracResultFile2 << "\\begin{tabular}{l c} " << "\n";
    FracResultFile2 << "\\hline" << "\n";
    FracResultFile2 << "\\hline" << "\n";
    FracResultFile2 << "Decay channel &  Fit fraction " << " \\\\ " << "\n";
    FracResultFile2 << "\\hline" << "\n";

    for(int i =0 ; i < fracs.size() ; i++){
        //if( fracs[i].err()<0.0001)continue;
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
            FracResultFile << fracs[i].err() * 100.  ;
            if(tot>0)FracResultFile << " \\pm " << sqrt(tot) * 100. ;
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
            FracResultFile2 << fracs[i].err() * 100.  ;
            if(tot>0)FracResultFile2 << " \\pm " << sqrt(tot) * 100. ;
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
    
    // Interference Fractions
    ofstream IntFracResultFile;
    IntFracResultFile.open(outDir+"intFracResult_table.tex",ofstream::trunc);

    IntFracResultFile << "\\begin{tabular}{l l r} " << "\n";
    IntFracResultFile << "\\hline" << "\n";
    IntFracResultFile << "\\hline" << "\n";
    IntFracResultFile << "Decay channel $i$ & Decay channel $j$ &  $IF_{ij} [ \\% ] $ " << " \\\\ " << "\n";
    IntFracResultFile << "\\hline" << "\n";
    
    for(int i =0 ; i < interferenceFracs.size() ; i++){
            if(abs(interferenceFracs[i].val() * 100.) < 0.1)continue;
            IntFracResultFile << fixed << setprecision(2) << "$" << starterFit.latexName(interferenceFracs[i].name()) << "$ & $" ;
            IntFracResultFile << interferenceFracs[i].val() * 100. << " \\pm " ;
            IntFracResultFile << interferenceFracs[i].err() * 100.  ;
            IntFracResultFile << "$ \\\\ " << "\n";
    }

    IntFracResultFile << "\\hline" << "\n";
    IntFracResultFile << "\\hline" << "\n";
    IntFracResultFile << "\\end{tabular}" << "\n";
    
    // Alt model table
    ofstream altModelsFile;
    altModelsFile.open(outDir+"altModels_table.tex",ofstream::trunc);

    altModelsFile << "\\begin{tabular}{l";
    for(int i=0;i<allModelFits.size();i++) altModelsFile << " r";
    altModelsFile << " } " << "\n";
    altModelsFile << "\\hline" << "\n";
    altModelsFile << "\\hline" << "\n";
    altModelsFile << "Decay channel $F_{i} [ \\% ] $ ";
    for(int i=0;i<allModelFits.size();i++) altModelsFile << " & Model " << i ;
    altModelsFile << " \\\\ " << "\n" << "\\hline" << "\n";
    
    vector<FitFraction> allAmps;
    MinuitParameterSet allPars;
    for(int i=0;i<allModelFits.size();i++){
        auto fracs_i = allModelFits[i].fitFractions();
        auto params_i = allModelFits[i].floating();

        for(auto &f : fracs_i){
            bool foundAmp = false;
            for(auto& f_allAmps : allAmps){
                if(removeLineshapeMods(f.name()) == removeLineshapeMods(f_allAmps.name())){
                    foundAmp = true;
                    break;
                }
            }
            if(!foundAmp){allAmps.emplace_back(removeLineshapeMods(f.name()),f.val(),f.err()); INFO(f.name());}
        }
        
        for(auto &p : params_i){
            allPars.addOrGet(p->name(),p->flag(),p->mean(),p->err());
        }
        
    }
    
    INFO("allAmps.size() = " << allAmps.size());
    //std::sort(allAmps.begin(), allAmps.end());
    //std::reverse(allAmps.begin(), allAmps.end());
    
    // Fit fractions
    for(int i=0;i<allAmps.size();i++){

        if( allAmps[i].name().find( "kMatrix" ) != string::npos || allAmps[i].name().find( "FOCUS" ) != string::npos || allAmps[i].name().find( "Sum_KPi" ) != string::npos || allAmps[i].name().find( "Sum_PiPi" ) != string::npos ) continue;
            
        altModelsFile << fixed << setprecision(2) << "$" << starterFit.latexName(allAmps[i].name()) << "$ " ;
                
        for(int j=0;j<allModelFits.size();j++){

            auto fracs_j = allModelFits[j].fitFractions();
            bool foundAmpName = false;

            for(auto &f : fracs_j){
                if(removeLineshapeMods(f.name())==removeLineshapeMods(allAmps[i].name())){
                    altModelsFile << " & $ " << f.val() * 100. << " \\pm " ;
                    altModelsFile << f.err() * 100.  ;
                    altModelsFile << " $ ";
                    foundAmpName = true;
                    break;
                }
            }
            if(!foundAmpName)altModelsFile << " & $-$ ";
        }
        
        altModelsFile << " \\\\ " << "\n";
    }
    
    // Fit parameters
    altModelsFile << "\\hline" << "\n";
    for(int i=0;i<allPars.size();i++){
            
        if( allPars[i]->name().find( "_Re" ) != string::npos || allPars[i]->name().find( "_Im" ) != string::npos ) continue;
        
        double scale = 1;
        if( allPars[i]->name().find( "_mass" ) != string::npos || allPars[i]->name().find( "_width" ) != string::npos )scale=1000.;
        altModelsFile << fixed << setprecision(1) << "$" << starterFit.latexName(allPars[i]->name()) << "$ " ;
                
        for(int j=0;j<allModelFits.size();j++){

            auto params_j = allModelFits[j].floating();
            bool foundPar = false;

            for(auto &p : params_j){
                if(p->name()==allPars[i]->name()){
                    altModelsFile << " & $ " << p->mean() * scale  << " \\pm " ;
                    altModelsFile << p->err() * scale  ;
                    altModelsFile << " $ ";
                    foundPar = true;
                    break;
                }
            }
            if(!foundPar)altModelsFile << " & $-$ ";
        }
        
        altModelsFile << " \\\\ " << "\n";
    }
    
    
    // Fit quality
    altModelsFile << "\\hline" << "\n";
    
    altModelsFile << fixed << setprecision(1) << "$ -2nLL $ " ;
    for(int j=0;j<allModelFits.size();j++){
        double dLL = allModelFits[j].LL()-allModelFits[0].LL();
        double sigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(abs(dLL),abs(allModelFits[j].nParam()-allModelFits[0].nParam()))) ;
        altModelsFile << " & $ " << dLL << " (" << sigma <<  "\\sigma ) $ ";
    }
    altModelsFile << " \\\\ " << "\n";
    
    altModelsFile << fixed << setprecision(2) << "$ N_{par} $ " ;
    for(int j=0;j<allModelFits.size();j++) altModelsFile << " & $ " << allModelFits[j].nParam() << " $ ";
    altModelsFile << " \\\\ " << "\n";
    
    altModelsFile << fixed << setprecision(2) << "$ \\chi^2/\\nu $ " ;
    for(int j=0;j<allModelFits.size();j++) altModelsFile << " & $ " << allModelFits[j].chi2()/allModelFits[j].dof() << " $ ";
    altModelsFile << " \\\\ " << "\n";
    
    altModelsFile << fixed << setprecision(2) << "$ \\chi^2/bin $ " ;
    for(int j=0;j<allModelFits.size();j++) altModelsFile << " & $ " << allModelFits[j].chi2()/allModelFits[j].nBins() << " $ ";
    altModelsFile << " \\\\ " << "\n";

    // End of table
    altModelsFile << "\\hline" << "\n";
    altModelsFile << "\\hline" << "\n";
    altModelsFile << "\\end{tabular}" << "\n";
}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  analyzeResults();

  return 0;
}
