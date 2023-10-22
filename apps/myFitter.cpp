#include <chrono>
#include <ctime>
#include <sstream>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <algorithm>   
#include <vector>       
#include <random>       
#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

#include "AmpGen/Kinematics.h"
#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/BkgPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"
#include "AmpGen/ErrorPropagator.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "AmpGen/LHCbStyle.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/Factory.h"
#if ENABLE_AVX2
    #include "AmpGen/EventListSIMD.h"
    using EventList_type = AmpGen::EventListSIMD;
#else
    #include "AmpGen/EventList.h"
    using EventList_type = AmpGen::EventList; 
#endif
#include "AmpGen/HyperPlot/HyperHistogram.h"

#include <TGraph.h>
#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

using namespace std;
using namespace AmpGen;

namespace AmpGen { 
  make_enum(pdfTypes, CoherentSum, PolarisedSum, PolarisedSumBkgBDT, IncoherentSum, FixedLib)
  make_enum(phspTypes, PhaseSpace, RecursivePhaseSpace, TreePhaseSpace)
}  

void makePlotWeightFile(PolarisedSum& sig, IncoherentSum& bkg, const EventList_type& eventsPlotMC, const LinearErrorPropagator& ep, TRandom3& rand){
    
    string outDir = NamedParameter<std::string>("outDir", ".");    
    string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName", "Fit_weights.root");  
    FitWeightFileName = outDir + "/" + FitWeightFileName;

    TFile* weight_file = TFile::Open(FitWeightFileName.c_str(),"RECREATE");
    weight_file->cd();
    TTree* weight_tree = new TTree("DalitzEventList","DalitzEventList");
    
    auto plot_amps = NamedParameter<string>("plot_amps", std::vector<string>(),"amplitudes to plot" ).getVector();
    auto plot_weights = NamedParameter<string>("plot_weights", std::vector<string>(),"plot weight names" ).getVector();
    
    vector<double> weights(plot_amps.size(),0.);
    vector<TBranch*> branches; 
    vector<vector<unsigned>> indices;

    auto weightFunction = sig.componentEvaluator();
    //auto evaluator_sig = sig.evaluator();

    // Amplitude component weights
    for(unsigned int i = 0; i < plot_amps.size(); i++){
        INFO( "Plotting amp " << plot_amps[i] << " with weight " <<  plot_weights[i] );

        auto branch = weight_tree->Branch( plot_weights[i].c_str(),&weights[i]);
        branches.push_back(branch);
        
        stringstream ss( plot_amps[i] );
        std::vector<std::string> selectAmps;
        while( ss.good() ){
            std::string substr;
            getline( ss, substr, '_' );
            selectAmps.push_back( AmpGen::programatic_name(substr) );
        }
        
        vector<unsigned> index;
        unsigned counter = 0;
        for( auto& key : weightFunction.keys ){
            stringstream ss( key );
            std::vector<std::string> selectKeys;
            while( ss.good() ){
                std::string substr;
                getline( ss, substr, '!' );
                selectKeys.push_back( substr );
                //INFO("Found " <<substr);
            }
            if(selectKeys.size()!=2)throw "ERROR";
            for(int a = 0; a<selectAmps.size();a++)for(int b = 0; b<selectAmps.size();b++){
                if ( ( (selectKeys[0].find(selectAmps[a]) != std::string::npos) && (selectKeys[1].find(selectAmps[b]) != std::string::npos) ) || plot_amps[i] == "") {
                    index.push_back(counter);
                    //if(selectAmp != "")INFO("Found " << key);
                }
            }
            counter++;
        }
        indices.push_back(index);
    }
    
    // Fit weights
    vector<double> weightsFit(3,0.);
    auto branch_w_full = weight_tree->Branch( "weight",&weightsFit[0]);
    branches.push_back(branch_w_full);
    auto branch_w_sig = weight_tree->Branch( "weight_sig",&weightsFit[1]);
    branches.push_back(branch_w_sig);
    auto branch_w_bkg = weight_tree->Branch( "weight_bkg",&weightsFit[2]);
    branches.push_back(branch_w_bkg);
    
    // Error bands
    auto params = ep.params();
    unsigned int nPermErrorBands   = NamedParameter<unsigned int>("nPermErrorBands",0);  
    auto gep = GaussErrorPropagator( ep.cov(), params, &rand );
    vector<double> weightsErr(nPermErrorBands,0.);
    vector<vector<double>> weightsErrVec; 

    // Sys error
    string sysCovFileName = NamedParameter<std::string>("sysCovFileName", "");  
    FitResult fr(sysCovFileName);
    auto params_sys = fr.floating();
    bool doSysErrBand = params_sys.size() > 0 ? true : false;
    INFO("SysCovFile has " << params_sys.size() << " parameters vs " << params.size() << " fit parameters " );

    for(unsigned int i=0; i < params_sys.size(); ++i){
        if(params[i]->name() != params_sys[i]->name() ) {
            ERROR("SysCovFile incosistent " << params.size() << " " << params_sys.size() << " " << params[i]->name() << " " << params_sys[i]->name() );
            doSysErrBand = false;
            break;
        }
    }
    auto gep_tot = GaussErrorPropagator( doSysErrBand ? ep.cov()+fr.getReducedCovariance() : ep.cov(), params, &rand );
    vector<double> weightsErrTot(nPermErrorBands,0.);
    vector<vector<double>> weightsErrTotVec; 
        
    if( nPermErrorBands > 0 ){
        for(int i = 0; i < nPermErrorBands; i++){
            auto branch = weight_tree->Branch( ("weightErr_"+to_string(i)).c_str(),&weightsErr[i]);     
            branches.push_back(branch);
            if(doSysErrBand){
                auto branch_tot = weight_tree->Branch( ("weightErrTot_"+to_string(i)).c_str(),&weightsErrTot[i]);     
                branches.push_back(branch_tot);
            }
        }
        for(int i = 0; i < nPermErrorBands; i++){
            //INFO("sig.norm()" << sig.norm());
            //INFO("sig.prob_unnormalisedNoCache()" << sig.getValNoCache(eventsPlotMC[1]) << endl);
            INFO("Doing error band perturb " << i);
            gep.reset();
            gep.perturb();
            sig.reset();   
            sig.prepare();
            vector<double> weightsErrTmp;
            for( const auto& evt : eventsPlotMC ){
                weightsErrTmp.push_back( ( sig.getValNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ) ; 
            }
            weightsErrVec.push_back(weightsErrTmp);

            if(doSysErrBand){
                gep_tot.reset();
                gep_tot.perturb();
                sig.reset();   
                sig.prepare();
                vector<double> weightsErrTotTmp;
                for( const auto& evt : eventsPlotMC ){
                    weightsErrTotTmp.push_back( ( sig.getValNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ) ; 
                }
                weightsErrTotVec.push_back(weightsErrTotTmp);
            }
        }
        gep.reset();
        if(doSysErrBand)gep_tot.reset();
        sig.reset();   
        sig.prepare();
    }

    unsigned counter = 0;
    for( const auto& evt : eventsPlotMC ){
        
        //weightsFit[1] = evt.weight() * evaluator_sig(evt) / evt.genPdf() ; 
        weightsFit[1] = ( sig.getValNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ; 
        weightsFit[2] = ( bkg.prob_unnormalisedNoCache(evt) * bkg.getWeight() / bkg.norm() ) * evt.weight() / evt.genPdf() ; 
        //weightsFit[2] = evt.weight() * evaluator_bkg(evt) / evt.genPdf() ;  // does not work ???
        weightsFit[0] = weightsFit[1] + weightsFit[2]; 
        
        auto weightFun = weightFunction(evt);
        for(int i = 0; i < plot_amps.size(); i++){
            weights[i] = 0;
            for( unsigned j = 0 ; j != indices[i].size(); ++j ) weights[i] += ( weightFun[indices[i][j]] * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ; 
        }

        for(int i = 0; i < nPermErrorBands; i++){
            weightsErr[i] = weightsErrVec[i][counter] + weightsFit[2] ;
            if(doSysErrBand)weightsErrTot[i] = weightsErrTotVec[i][counter] + weightsFit[2];
        }
        
        counter++;        
        weight_tree->Fill();
    }

    weight_tree->Write();
    weight_file->Write();
    weight_file->Close();
    std::cout << "Fit weights saved" << std::endl;    
}

void makePlotWeightFile(PolarisedSum& sig, BkgPDF& bkg, const EventList_type& eventsPlotMC, const LinearErrorPropagator& ep, TRandom3& rand){
    
    string outDir = NamedParameter<std::string>("outDir", ".");
    string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName", "Fit_weights.root");
    FitWeightFileName = outDir + "/" + FitWeightFileName;

    TFile* weight_file = TFile::Open(FitWeightFileName.c_str(),"RECREATE");
    weight_file->cd();
    TTree* weight_tree = new TTree("DalitzEventList","DalitzEventList");
    
    auto plot_amps = NamedParameter<string>("plot_amps", std::vector<string>(),"amplitudes to plot" ).getVector();
    auto plot_weights = NamedParameter<string>("plot_weights", std::vector<string>(),"plot weight names" ).getVector();
    
    vector<double> weights(plot_amps.size(),0.);
    vector<TBranch*> branches;
    vector<vector<unsigned>> indices;

    auto weightFunction = sig.componentEvaluator();
    //auto evaluator_sig = sig.evaluator();

    // Amplitude component weights
    for(unsigned int i = 0; i < plot_amps.size(); i++){
        INFO( "Plotting amp " << plot_amps[i] << " with weight " <<  plot_weights[i] );

        auto branch = weight_tree->Branch( plot_weights[i].c_str(),&weights[i]);
        branches.push_back(branch);
        
        stringstream ss( plot_amps[i] );
        std::vector<std::string> selectAmps;
        while( ss.good() ){
            std::string substr;
            getline( ss, substr, '_' );
            selectAmps.push_back( AmpGen::programatic_name(substr) );
        }
        
        vector<unsigned> index;
        unsigned counter = 0;
        for( auto& key : weightFunction.keys ){
            stringstream ss( key );
            std::vector<std::string> selectKeys;
            while( ss.good() ){
                std::string substr;
                getline( ss, substr, '!' );
                selectKeys.push_back( substr );
                //INFO("Found " <<substr);
            }
            if(selectKeys.size()!=2)throw "ERROR";
            for(int a = 0; a<selectAmps.size();a++)for(int b = 0; b<selectAmps.size();b++){
                if ( ( (selectKeys[0].find(selectAmps[a]) != std::string::npos) && (selectKeys[1].find(selectAmps[b]) != std::string::npos) ) || plot_amps[i] == "") {
                    index.push_back(counter);
                    //if(selectAmp != "")INFO("Found " << key);
                }
            }
            counter++;
        }
        indices.push_back(index);
    }
    
    // Fit weights
    vector<double> weightsFit(3,0.);
    auto branch_w_full = weight_tree->Branch( "weight",&weightsFit[0]);
    branches.push_back(branch_w_full);
    auto branch_w_sig = weight_tree->Branch( "weight_sig",&weightsFit[1]);
    branches.push_back(branch_w_sig);
    auto branch_w_bkg = weight_tree->Branch( "weight_bkg",&weightsFit[2]);
    branches.push_back(branch_w_bkg);
    
    // Error bands
    auto params = ep.params();
    unsigned int nPermErrorBands   = NamedParameter<unsigned int>("nPermErrorBands",0);
    auto gep = GaussErrorPropagator( ep.cov(), params, &rand );
    vector<double> weightsErr(nPermErrorBands,0.);
    vector<vector<double>> weightsErrVec;

    // Sys error
    string sysCovFileName = NamedParameter<std::string>("sysCovFileName", "");
    FitResult fr(sysCovFileName);
    auto params_sys = fr.floating();
    bool doSysErrBand = params_sys.size() > 0 ? true : false;
    INFO("SysCovFile has " << params_sys.size() << " parameters vs " << params.size() << " fit parameters " );

    for(unsigned int i=0; i < params_sys.size(); ++i){
        if(params[i]->name() != params_sys[i]->name() ) {
            ERROR("SysCovFile incosistent " << params.size() << " " << params_sys.size() << " " << params[i]->name() << " " << params_sys[i]->name() );
            doSysErrBand = false;
            break;
        }
    }
    auto gep_tot = GaussErrorPropagator( doSysErrBand ? ep.cov()+fr.getReducedCovariance() : ep.cov(), params, &rand );
    vector<double> weightsErrTot(nPermErrorBands,0.);
    vector<vector<double>> weightsErrTotVec;
        
    if( nPermErrorBands > 0 ){
        for(int i = 0; i < nPermErrorBands; i++){
            auto branch = weight_tree->Branch( ("weightErr_"+to_string(i)).c_str(),&weightsErr[i]);
            branches.push_back(branch);
            if(doSysErrBand){
                auto branch_tot = weight_tree->Branch( ("weightErrTot_"+to_string(i)).c_str(),&weightsErrTot[i]);
                branches.push_back(branch_tot);
            }
        }
        for(int i = 0; i < nPermErrorBands; i++){
            //INFO("sig.norm()" << sig.norm());
            //INFO("sig.prob_unnormalisedNoCache()" << sig.getValNoCache(eventsPlotMC[1]) << endl);
            INFO("Doing error band perturb " << i);
            gep.reset();
            gep.perturb();
            sig.reset();
            sig.prepare();
            vector<double> weightsErrTmp;
            for( const auto& evt : eventsPlotMC ){
                weightsErrTmp.push_back( ( sig.getValNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ) ;
            }
            weightsErrVec.push_back(weightsErrTmp);

            if(doSysErrBand){
                gep_tot.reset();
                gep_tot.perturb();
                sig.reset();
                sig.prepare();
                vector<double> weightsErrTotTmp;
                for( const auto& evt : eventsPlotMC ){
                    weightsErrTotTmp.push_back( ( sig.getValNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ) ;
                }
                weightsErrTotVec.push_back(weightsErrTotTmp);
            }
        }
        gep.reset();
        if(doSysErrBand)gep_tot.reset();
        sig.reset();
        sig.prepare();
    }

    unsigned counter = 0;
    for( const auto& evt : eventsPlotMC ){
        
        //weightsFit[1] = evt.weight() * evaluator_sig(evt) / evt.genPdf() ;
        weightsFit[1] = ( sig.getValNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ;
        weightsFit[2] = ( bkg.prob_unnormalisedNoCache(evt) * bkg.getWeight() / bkg.norm() ) * evt.weight() / evt.genPdf() ;
        //weightsFit[2] = evt.weight() * evaluator_bkg(evt) / evt.genPdf() ;  // does not work ???
        weightsFit[0] = weightsFit[1] + weightsFit[2];
        
        auto weightFun = weightFunction(evt);
        for(int i = 0; i < plot_amps.size(); i++){
            weights[i] = 0;
            for( unsigned j = 0 ; j != indices[i].size(); ++j ) weights[i] += ( weightFun[indices[i][j]] * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ;
        }

        for(int i = 0; i < nPermErrorBands; i++){
            weightsErr[i] = weightsErrVec[i][counter] + weightsFit[2] ;
            if(doSysErrBand)weightsErrTot[i] = weightsErrTotVec[i][counter] + weightsFit[2];
        }
        
        counter++;
        weight_tree->Fill();
    }

    weight_tree->Write();
    weight_file->Write();
    weight_file->Close();
    std::cout << "Fit weights saved" << std::endl;
}

void makePlotWeightFile(IncoherentSum& sig, const EventList_type& eventsPlotMC, const LinearErrorPropagator& ep, TRandom3& rand){
    string outDir = NamedParameter<std::string>("outDir", ".");    
    string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName", "Fit_weights.root");  
    FitWeightFileName = outDir + "/" + FitWeightFileName;

    TFile* weight_file = TFile::Open(FitWeightFileName.c_str(),"RECREATE");
    weight_file->cd();
    TTree* weight_tree = new TTree("DalitzEventList","DalitzEventList");
    
    double weight;    
    vector<TBranch*> branches; 
    branches.push_back(weight_tree->Branch( "weight", &weight));

    // Error bands
    unsigned int nPermErrorBands   = NamedParameter<unsigned int>("nPermErrorBands",0);  
    auto gep = GaussErrorPropagator( ep.cov(), ep.params(), &rand );
    vector<double> weightsErr(nPermErrorBands,0.);
    vector<vector<double>> weightsErrVec; 

    if( nPermErrorBands > 0 ){
        for(int i = 0; i < nPermErrorBands; i++){
            auto branch = weight_tree->Branch( ("weightErr_"+to_string(i)).c_str(),&weightsErr[i]);     
            branches.push_back(branch);
        }
        for(int i = 0; i < nPermErrorBands; i++){
            //INFO("sig.norm()" << sig.norm());
            //INFO("sig.prob_unnormalisedNoCache()" << sig.prob_unnormalisedNoCache(eventsPlotMC[1]) << endl);
            gep.perturb();
            sig.reset(false);   
            sig.prepare();
            vector<double> weightsErrTmp;
            for( const auto& evt : eventsPlotMC ){
                weightsErrTmp.push_back( ( sig.prob_unnormalisedNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ) ; 
            }
            weightsErrVec.push_back(weightsErrTmp);
        }
        gep.reset();
        sig.reset(false);   
        sig.prepare();
    }
    
    unsigned counter = 0;
    for( const auto& evt : eventsPlotMC ){
        weight = ( sig.prob_unnormalisedNoCache(evt) * sig.getWeight() / sig.norm() ) * evt.weight() / evt.genPdf() ; 
        
        for(int i = 0; i < nPermErrorBands; i++){
            weightsErr[i] = weightsErrVec[i][counter];
        }
        counter++;
        
        weight_tree->Fill();
    }

    weight_tree->Write();
    weight_file->Write();
    weight_file->Close();
    std::cout << "Fit weights saved" << std::endl;    
}

std::vector<ThreeBodyCalculator> threeBodyCalculators( MinuitParameterSet& mps )
{
  std::vector<std::string> threeBodiesToIntegrate = NamedParameter<std::string>( "ThreeBodiesToIntegrate" ).getVector();
  std::vector<ThreeBodyCalculator> calculators;
  for ( auto& v : threeBodiesToIntegrate ){
      auto spline_params = NamedParameter<double>( v + "::Spline").getVector();
      size_t nBins;
      double min, max;
      if( spline_params.size() == 3 ){
          nBins = size_t( spline_params[0] );
          min   =         spline_params[1] ; 
          max   =         spline_params[2];
      }
      else {
          nBins = NamedParameter<double>( v + "::Spline::N"  , 0. );
          min   = NamedParameter<double>( v + "::Spline::Min", 0. );
          max   = NamedParameter<double>( v + "::Spline::Max", 0. );
      }
      calculators.emplace_back( v, mps,nBins,min,max);
      INFO(v);
      INFO(nBins << " " << min << " " << max  << endl);

  }
  return calculators;
}

void randomizeStartingPoint( MinuitParameterSet& MPS, TRandom3& rand, bool SplineOnly = false )
{
  for (auto& param : MPS) {
    if ( param->isFree() == 0 ) continue;
    if ( param->name().find( "cut_dim" ) != std::string::npos) continue;
    if ( SplineOnly && param->name().find( "::Spline" ) == std::string::npos ) continue;
    
    if (param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos 
        //|| param->name().find( "_mass" )  != std::string::npos  || param->name().find( "_width" ) != std::string::npos
        //|| param->name().find( "_alpha" ) != std::string::npos  || param->name().find( "_beta" ) != std::string::npos
    ){         
        
        double min = param->minInit();
        double max = param->maxInit();
        //if(param->name().find( "_Re" ) != std::string::npos){
           //min = 0;
           // max = 2;
        //}
        if((min == 0 && max == 0) || min == max){
            INFO(min << " " << max );
            min = param->mean()-100*param->stepInit();
            max = param->mean()+100*param->stepInit();
            INFO(min << " " << max );
        } 
        
        double new_val = rand.Uniform(min+5*param->stepInit(),max-5*param->stepInit()); 
        param->setInit(new_val);
        param->setCurrentFitVal(new_val);
        INFO( param->name() << " = " << param->mean() );
    }
  }
}

void perturb( MinuitParameterSet& MPS, double sigma = 1)
{
    for ( auto& param : MPS ) {
        if ( !param->isFree() ) continue;
        if ( param->name().find( "cut_dim" ) != std::string::npos) continue;
        if ( param->name().find( "_Re" ) == std::string::npos && param->name().find( "_Im" ) == std::string::npos ) continue;
            
        double new_val = gRandom->Gaus( param->mean(), param->err()* sigma);
        param->setInit(new_val);
        param->setCurrentFitVal( new_val );
    }
}

void vary( MinuitParameterSet& MPS, string& paramName, double sigma = 1)
{
    for ( auto& param : MPS ) {

        if ( param->name().find( paramName ) == std::string::npos ) continue;
            
        double new_val = param->mean() + param->err() * sigma;
        INFO("Setting parameter " << paramName << " from " << param->mean() << " to new val = " << new_val);
        
        param->setInit(new_val);
        param->setCurrentFitVal( new_val );
        return;
    }
    ERROR("Paramater " << paramName << " not found");
}

void varyGauss( MinuitParameterSet& MPS, vector<string>& paramNames, double sigma = 1)
{
    for ( auto& param : MPS ) {

        for ( auto& paramName : paramNames ) {
        
            if ( param->name().find( paramName ) == std::string::npos ) continue;
                
            double new_val = gRandom->Gaus( param->mean(), param->err()* sigma);
            INFO("Setting parameter " << paramName << " from " << param->mean() << " to new val = " << new_val);
            
            param->setInit(new_val);
            param->setCurrentFitVal( new_val );
        }
        
    }
}

void scan( MinuitParameterSet& MPS, string name, double min, double max, int step, int nSteps)
{
    if(step>=nSteps && name.find( "_mass" ) != std::string::npos){
        for ( auto& param : MPS ){
            TString width_name = TString(name).ReplaceAll("_mass","_width");
            if ( param->name().find(width_name) == std::string::npos ) continue;
            
            double new_val = param->mean() * (step/nSteps) * 2.;
            param->setInit(new_val);
            param->setCurrentFitVal(new_val);
            param->fix();
            
            cout << "Set " << width_name << " to " << new_val << endl;
            step = step % nSteps;
        }
    }
    
    for ( auto& param : MPS ) {
        if ( param->name().find( name ) == std::string::npos ) continue;
            
        double new_val = min + (max-min)*(double)step/(double)nSteps;
        param->setInit(new_val);
        param->setCurrentFitVal( new_val );
        param->fix();
        cout << "Scaning over " << name << ": set val to " << new_val << endl; 
    }
}

void setSplineVals(MinuitParameterSet& mps, string& head){
    auto fixBoundaries = NamedParameter<int>( "setSplineVals::fixBoundaries", 0);
    auto spline_params = NamedParameter<double>( head + "::Spline").getVector();
    int nBins = int( spline_params[0] );
    double min = spline_params.size() == 4 ? spline_params[1] * spline_params[1] : spline_params[1];
    double max = spline_params.size() == 4 ? spline_params[2] * spline_params[2] : spline_params[2];

    INFO("Set spline vals for " << head);
    double m = mps.find(head+"_mass")->mean();
    double gamma = mps.find(head+"_width")->mean();
    INFO("Use BW with m = " << m << " g= " << gamma);
    
    int from = fixBoundaries==1 ? 1 : 0;
    int to = fixBoundaries==1 ? nBins-1 : nBins;
    for(int i= from; i<to; i++){
        double s = min + (max-min)* i/((double)nBins-1.);
        //complex<double> BW = -complex<double>(0,1) * m * gamma/(m*m - s -  complex<double>(0,1) * m * gamma);
        complex<double> BW =  m * gamma/(m*m - s -  complex<double>(0,1) * sqrt(s) * gamma);
        mps.find(head+"::Spline::Re::"+to_string(i))->setInit(abs(BW));
        mps.find(head+"::Spline::Im::"+to_string(i))->setInit(arg(BW)*180./3.141);
    }
    
    if(fixBoundaries==2){
        mps.find(head+"::Spline::Re::"+to_string(0))->setInit(0.);
        //mps.find(head+"::Spline::Im::"+to_string(0))->setInit( mps.find(head+"::Spline::Im::"+to_string(1))->mean() );
        mps.find(head+"::Spline::Re::"+to_string(0))->fix();
        //mps.find(head+"::Spline::Im::"+to_string(0))->fix();

        mps.find(head+"::Spline::Re::"+to_string(nBins-1))->setInit(0.);
        //mps.find(head+"::Spline::Im::"+to_string(nBins-1))->setInit( mps.find(head+"::Spline::Im::"+to_string(nBins-2))->mean() );
        mps.find(head+"::Spline::Re::"+to_string(nBins-1))->fix();
        //mps.find(head+"::Spline::Im::"+to_string(nBins-1))->fix();
    }
    
}

void sanityChecks(MinuitParameterSet& mps){

   for(int i=0;i<mps.size();++i)if(mps[i]->name().find( "cut_dim" ) != std::string::npos){ mps.unregister( mps.at(i)); i=0; }

   // check 
   for(int i=0;i<mps.size();++i)
       if(mps[i]->name().find( "::Spline" ) != std::string::npos)
           if(mps[i]->name().find( "::Re" ) == std::string::npos && mps[i]->name().find( "::Im" ) == std::string::npos){ mps.unregister( mps.at(i)); i=0; }
    
   string head = NamedParameter<std::string>("Head","");
   int nBins = 0;
   if(head!=""){
       auto spline_params = NamedParameter<double>( head + "::Spline").getVector();
       nBins = int( spline_params[0] );
    }
    for(int i=0;i<mps.size();i++){
               bool found = false;                          
               
               auto param = mps[i];
               if(!param->isFree())continue;
               if(param->name().find( "::Spline::") != std::string::npos){   
                   if(param->name().find( head + "::Spline::") != std::string::npos) found = true;
                   TString name(param->name());
                   name.ReplaceAll(head + "::Spline::Re::","");
                   name.ReplaceAll(head + "::Spline::Im::","");
                   if(atoi(name)>=nBins)found = false;
               }
               else found = true;
               
               if(found == false){    
                   INFO("Fitting " << mps[i]->name() << " but there is no matching amplitude, remove it from MinuitParameterSet");
                   mps.unregister(mps[i]);  
                   i=0;                
               }
   }
    
   for(int i=0;i<mps.size();i++){
       if((int)mps[i]->flag()==3)continue;
       if((mps[i]->name().find( "_mass" ) != std::string::npos || mps[i]->name().find( "_width" ) != std::string::npos || mps[i]->name().find( "_alpha" ) != std::string::npos || mps[i]->name().find( "_beta" ) != std::string::npos ) ){
           TString name(mps[i]->name());
           name.ReplaceAll("_mass","");
           name.ReplaceAll("_width","");     
           name.ReplaceAll("_alpha","");     
           name.ReplaceAll("_beta","");     
           
           bool found = false;           
           for(int j=0;j<mps.size();j++){               
               if(  mps[j]->name().find(name) != std::string::npos && mps[j]->name().find( "_Re" ) != std::string::npos ) found=true;
           }  
           if(found == false){    
               INFO("Fitting " << mps[i]->name() << " but there is no matching amplitude, remove it from MinuitParameterSet");
               mps.unregister(mps[i]);  
               i=0;                
           }
        }
   }  
}

void checkAmps(PolarisedSum& sig, MinuitParameterSet& mps){
    for(int i=0;i<mps.size();++i){
        string name = mps[i]->name();
        if(name.find( "::pole::" ) != std::string::npos)continue;
        if(name.find( "Inco" ) != std::string::npos)continue;
        if(name.find( "_Re" ) != std::string::npos || name.find( "_Im" ) != std::string::npos) {            
            vector<string>name_split = split(name,',');
            //INFO(name_split[0]);
            //INFO(name_split[1]);
            //auto mod = NamedParameter<std::string>("Particle::DefaultModifier","");
            if(name.find( "GSpline.BL" ) != std::string::npos)continue;
            if(name.find( "Omega.BL" ) != std::string::npos)continue;
            if(sig.findAmp(name_split[0])==0){
                INFO("Removing " << name);
                mps.unregister(mps[i]);
                i=0;
            }
        }
    }
}

struct phsp_cut {
    phsp_cut(std::vector<unsigned int> dim, std::vector<double> limits, bool invertCut = false):_dim(dim),_limits(limits),_invertCut(invertCut){}
    bool operator()(const Event& evt){
        if(_dim.size()==0)return true;
        if(sqrt(evt.s(_dim)) > _limits[0] && sqrt(evt.s(_dim)) < _limits[1] )return !_invertCut;
        else return _invertCut;
    }
    private:
      std::vector<unsigned int> _dim;
      std::vector<double> _limits;
      bool _invertCut;
};

struct rnd_cut {
    rnd_cut(double fraction,bool random, int N):_fraction(fraction),_random(random),_counter(0),_N(N){
        INFO("Keeping " << _fraction * _N << " events (" << _fraction*100 << " %)" );
    }
    bool operator()(const Event& evt){
        
        if(_random){        
            if( gRandom->Rndm()>_fraction)return true;
            else return false;
        }
        else {
            if(_counter > (1-_fraction) * _N ){
                INFO(_counter);
                return true;
                }
            else{
                _counter++;
                return false;
                }
        }
    }
    private:
      double _fraction;
      bool _random;
      int _counter;
      int _N;
};

void addGaussianConstraint( Minimiser& mini, MinuitParameterSet& mps )
{
    std::vector<std::string> llConfigs = NamedParameter<std::string>( "GaussianConstraint",std::vector<std::string>() ).getVector();
    if(llConfigs.size()==0)return;
    
    for ( auto& ll_config : llConfigs ) {
        
        for(int i=0;i<mps.size();++i){
            string name = mps[i]->name();
            if(name == ll_config){
                auto ll_term = new GaussianConstraint();
                ll_term->configure( ll_config, mps );
                mini.addExtendedTerm( ll_term );
                continue;
            }
        }
    }
}

void addLASSO( Minimiser& mini, PolarisedSum& sig, MinuitParameterSet& mps, double lambda )
{
    auto ll_term = new LASSO();
    ll_term->configure( "lambda " + to_string(lambda), sig, mps );
    mini.addExtendedTerm( ll_term );
    INFO("Added LASSO term with lambda = " << lambda);
}

template <typename likelihoodType>
FitResult* doFit( likelihoodType&& likelihood, EventList_type& data, EventList_type& mc, MinuitParameterSet& MPS, PolarisedSum& sig, double lambda = -1 )
{
    auto time_wall = std::chrono::high_resolution_clock::now();
    auto time      = std::clock();
    /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
    Minimiser mini( likelihood, &MPS );
    addGaussianConstraint( mini, MPS );
    if(lambda > 0)addLASSO( mini, sig, MPS, lambda );

    auto threeBodyShapes     = threeBodyCalculators( MPS );
    unsigned int updateWidth = NamedParameter<unsigned int>( "UpdateWidth", 0 );
    
    //for ( auto& shape : threeBodyShapes )shape.debug(2,0.5);
    
    //std::vector<TGraph*> rw_old;
    //for( auto& shape : threeBodyShapes ) rw_old.push_back(shape.widthGraph(1.));
    //for( auto& graph : rw_old) graph->Write();
    
    if ( updateWidth ) {
        for ( auto& shape : threeBodyShapes ) shape.updateRunningWidth( MPS );
    }
    unsigned int nIterations            = NamedParameter<unsigned int>( "nIterationsWidth", 0 );
    std::vector<std::string> SlowParams = NamedParameter<std::string>( "Release", "" ).getVector();
    std::vector<MinuitParameter*> slowParamPtrs;
    if ( nIterations != 0 ) {
        for ( auto& param : SlowParams ) {
            auto it      = MPS.find( param );
            if ( it != nullptr ) {
                slowParamPtrs.push_back( it );
                it->fix();
            } else {
                WARNING( "Trying to release non-existent parameter: " << param );
            }
        }
    }
    INFO( "Fitting PDF, iterating width " << " " << nIterations + 1 << " times" );
    for ( unsigned int iteration = 0; iteration < nIterations + 1; ++iteration ) {
        mini.doFit();
        if ( iteration == 0 && nIterations != 0 ) {
            for ( auto& shape : threeBodyShapes ) shape.updateRunningWidth( MPS );
            for ( auto& param : slowParamPtrs ) param->setFree(); /// release the parameter ///
        }
    }
    
    //std::vector<TGraph*> rw_new;
    //for( auto& shape : threeBodyShapes ) rw_new.push_back(shape.widthGraph(1.));
    //for( auto& graph : rw_new) graph->Write();
    
    unsigned int nReTries = NamedParameter<unsigned int>( "Fit::nReTries", 5 );
    int nTries = 1;
    while(mini.status()>0 && nTries-1 < nReTries){
            INFO("Fit not converged, try again with small perturbation");
            perturb(MPS, nTries);
            mini.doFit();
            nTries++;
    }
    
    FitResult* fr = new FitResult(mini);
    
    auto twall_end  = std::chrono::high_resolution_clock::now();
    double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
    double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
    INFO("Done after " << nTries << " iteration");
    if(mini.status() == 0)INFO("Fit conveged");
    else if(mini.status() == 1)INFO("Valid minimum, Covar was made pos def");
    else INFO("Fit failed");
    INFO( "Wall time = " << tWall / 1000. );
    INFO( "CPU  time = " << time_cpu );
    
    return fr;
}

double cosTheta( const Event& evt ){
    
    TLorentzVector p0 = pFromEvent( evt, 0 ); //psi
    
    TLorentzVector p1 = pFromEvent( evt, 2 ); //pip
    TLorentzVector p2 = pFromEvent( evt, 3 ); //pim
    TLorentzVector p3 = pFromEvent( evt, 1 ); //Kp
    
    TLorentzVector pR = p1 + p2 + p3;
    p0.Boost( -pR.BoostVector() );
    p1.Boost( -pR.BoostVector() );
    p2.Boost( -pR.BoostVector() );
    p3.Boost( -pR.BoostVector() );
    
    TVector3 ez = -( p0.Vect() ).Unit();
    TVector3 n = ( p1.Vect().Cross( p2.Vect() ) ).Unit();
    
    return ez.Dot(n);    
}

double chi( const Event& evt ){
    
    //EventType B+ psi(2S)0 K+ pi+ pi- 
    
    TLorentzVector p0 = pFromEvent( evt, 0 ); //psi
    TLorentzVector p1 = pFromEvent( evt, 2 ); //pip
    TLorentzVector p2 = pFromEvent( evt, 3 ); //pim
    TLorentzVector p3 = pFromEvent( evt, 1 ); //Kp
    
    TLorentzVector pR = p1 + p2 + p3;
    p0.Boost( -pR.BoostVector() );
    p1.Boost( -pR.BoostVector() );
    p2.Boost( -pR.BoostVector() );
    p3.Boost( -pR.BoostVector() );
    
    TVector3 ez = -( p0.Vect() ).Unit();
    TVector3 n = ( p1.Vect().Cross( p2.Vect() ) ).Unit();
    
    double cosChi = - ( ( n.Cross(p1.Vect()) ).Unit() ).Dot( ( n.Cross(p0.Vect()) ).Unit() );
    double sinChi = - ( ( n.Cross(p1.Vect()) ).Unit() ).Cross( ( n.Cross(p0.Vect()) ).Unit() ).Dot(n);
    
    return TMath::ATan2( sinChi, cosChi );
}

inline double cosThetaMuAngle(const Event& evt){

    //EventType B+ K+ pi+ pi- mu+ mu-
    TLorentzVector p0 = pFromEvent( evt, 0 );
    TLorentzVector p1 = pFromEvent( evt, 1 );
    TLorentzVector p2 = pFromEvent( evt, 2 );
    TLorentzVector p3 = pFromEvent( evt, 3 );
    TLorentzVector p4 = pFromEvent( evt, 4 );
    
    TLorentzVector pR = p3 + p4;
    p0.Boost( -pR.BoostVector() );
    p1.Boost( -pR.BoostVector() );
    p2.Boost( -pR.BoostVector() );
    p3.Boost( -pR.BoostVector() );
    p4.Boost( -pR.BoostVector() );
    
    TVector3 pKs = (p0+p1+p2).Vect().Unit();
    
    return -(pKs).Dot(p3.Vect().Unit());
}

inline double chiMuAngle(const Event& evt){
    
    //EventType B+ K+ pi+ pi- mu+ mu-
    TLorentzVector p0 = pFromEvent( evt, 0 );
    TLorentzVector p1 = pFromEvent( evt, 1 );
    TLorentzVector p2 = pFromEvent( evt, 2 );
    TLorentzVector p3 = pFromEvent( evt, 3 );
    TLorentzVector p4 = pFromEvent( evt, 4 );
    
    TLorentzVector pB = p0 + p1 + p2 + p3 + p4;
    p0.Boost( -pB.BoostVector() );
    p1.Boost( -pB.BoostVector() );
    p2.Boost( -pB.BoostVector() );
    p3.Boost( -pB.BoostVector() );
    p4.Boost( -pB.BoostVector() );
    
    TVector3 pKs = (p0+p1+p2).Vect().Unit();
    TVector3 pPsi = (p3+p4).Vect().Unit();
    
    TVector3 aK = p0.Vect() - p0.Vect().Dot(pKs) * pKs;
    TVector3 aMu = p3.Vect() - p3.Vect().Dot(pPsi) * pPsi;
    
    double cos = aK.Dot(aMu);
    double sin = pPsi.Cross(aK).Dot(aMu);
    
    return TMath::ATan2( sin, cos );
}

vector<double> getChi2(const EventList_type& dataEvents, const EventList_type& mcEvents, const std::function<double( const Event& )>& fcn, const int dim = 5, const int minEventsPerBin = 25, double minBinWidth = 0., int mode = 0 ){
    
    EventType pdg = dataEvents.eventType();
    //EventType B+ psi(2S)0 K+ pi+ pi- 
    vector<unsigned int> m123{1,2,3};
    vector<unsigned int> m13{1,3};
    vector<unsigned int> m23{2,3};
    vector<unsigned int> m023{0,2,3};
    vector<unsigned int> m02{0,2};
    vector<unsigned int> m03{0,3};
    vector<unsigned int> m01{0,1};
    vector<unsigned int> m013{0,1,3};
    vector<unsigned int> m012{0,1,2};

    //return {evt.s( 1, 2, 3 ), evt.s( 0, 1 ), evt.s( 0, 2 ), evt.s( 2, 3 ), evt.s( 0, 1, 2 )};

    HyperPointSet points( dim );
    HyperPoint min(dim);
    HyperPoint max(dim);

    if(dim==7){
        min = HyperPoint( pdg.minmax(m123).first, pdg.minmax(m13).first, pdg.minmax(m23).first, pdg.minmax(m023).first, pdg.minmax(m02).first, -1, -3.141 );
        max = HyperPoint( pdg.minmax(m123).second, pdg.minmax(m13).second, pdg.minmax(m23).second, pdg.minmax(m023).second, pdg.minmax(m02).second, 1, 3.141 );
    }
    else if(dim==5 && mode ==0){
        min = HyperPoint( pdg.minmax(m123).first, pdg.minmax(m13).first, pdg.minmax(m23).first, pdg.minmax(m023).first, pdg.minmax(m02).first );
        max = HyperPoint( pdg.minmax(m123).second, pdg.minmax(m13).second, pdg.minmax(m23).second, pdg.minmax(m023).second, pdg.minmax(m02).second );        
    }
    else if(dim==5){
        min = HyperPoint( pow(pdg.minmax(m123).first,2), pow(pdg.minmax(m01).first,2), pow(pdg.minmax(m02).first,2), pow(pdg.minmax(m23).first,2), pow(pdg.minmax(m012).first,2) );
        max = HyperPoint( pow(pdg.minmax(m123).second,2), pow(pdg.minmax(m01).second,2), pow(pdg.minmax(m02).second,2), pow(pdg.minmax(m23).second,2), pow(pdg.minmax(m012).second,2) );        
    }
    
    HyperCuboid limits(min, max );    
    
    for( auto& evt : dataEvents ){
        HyperPoint point( dim );
        if(mode ==0){
            point.at(0)= sqrt(evt.s(m123));
            point.at(1)= sqrt(evt.s(m13));
            point.at(2)= sqrt(evt.s(m23));
            point.at(3)= sqrt(evt.s(m023));
            point.at(4)= sqrt(evt.s(m02));
        }
        else{
            point.at(0)= evt.s(m123);
            point.at(1)= evt.s(m01);
            point.at(2)= evt.s(m02);
            point.at(3)= evt.s(m23);
            point.at(4)= evt.s(m012);
        }
        if(dim==7){
            point.at(5)= cosTheta(evt);
            point.at(6)= chi(evt);
        }
        point.addWeight(evt.weight());
        points.push_back(point);
    }

    HyperPointSet pointsMC( dim);
    for( auto& evt : mcEvents ){
        HyperPoint point( dim );
        if(mode ==0){
            point.at(0)= sqrt(evt.s(m123));
            point.at(1)= sqrt(evt.s(m13));
            point.at(2)= sqrt(evt.s(m23));
            point.at(3)= sqrt(evt.s(m023));
            point.at(4)= sqrt(evt.s(m02));
        }
        else{
            point.at(0)= evt.s(m123);
            point.at(1)= evt.s(m01);
            point.at(2)= evt.s(m02);
            point.at(3)= evt.s(m23);
            point.at(4)= evt.s(m012);
        }
        if(dim==7){
            point.at(5)= cosTheta(evt);
            point.at(6)= chi(evt);            
        }
        point.addWeight(fcn( evt ) * evt.weight() / evt.genPdf());
        pointsMC.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                            /*** Name of the binning algorithm you want to use     */
                            HyperBinningAlgorithms::SMART_MULTI,
                            /***  The minimum number of events allowed in each bin */
                            /***  from the HyperPointSet provided (points1)        */
                            AlgOption::MinBinContent      (minEventsPerBin),                    
                            /*** This minimum bin width allowed. Can also pass a   */
                            /*** HyperPoint if you would like different min bin    */
                            /*** widths for each dimension                         */
                            AlgOption::MinBinWidth        (minBinWidth),                                                 
                            /*** If you want to use the sum of weights rather than */
                            /*** the number of events, set this to true.           */    
                            AlgOption::UseWeights         (true),                         
                            /*** Some algorithms use a random number generator. Set*/
                            /*** the seed here                                     */
                            AlgOption::RandomSeed         (1),                         
                            /*** What dimesnion would you like to split first? Only*/
                            /*** applies to certain algortihms                     */
                            AlgOption::StartDimension     (0)
                            /*** What dimesnions would you like to bin in?         */
                            //AlgOption::BinningDimensions  (binningDims),                      
                            /*** Setting this option will make the agorithm draw   */
                            /*** the binning scheme at each iteration              */
                            //AlgOption::DrawAlgorithm("Algorithm")                 
                            );
        //    dataHist.save("histData.root");
    //     HyperHistogram binningHist("histData.root",5);    
    //     HyperHistogram dataHist( binningHist.getBinning() );
    //     dataHist.fill(points); 
    
    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    mcHist.normalise(dataHist.integral());
    
    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();    
    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    
    vector<double> vals;
    vals.push_back(chi2);
    vals.push_back((double)nBins);
    
    return vals;
}

vector<double> getChi2MuMu(const EventList_type& dataEvents, const EventList_type& mcEvents, const std::function<double( const Event& )>& fcn, const std::function<double( const Event& )>& fcn_bkg, const int dim = 7, const int minEventsPerBin = 25, double minBinWidth = 0., int mode = 0, unsigned int filter_plotN = 0 ){
    
    EventType pdg = dataEvents.eventType();
    //EventType B+ K+ pi+ pi- mu+ mu-

    vector<unsigned int> m012{0,1,2};
    vector<unsigned int> m02{0,2};
    vector<unsigned int> m12{1,2};
    vector<unsigned int> m3412{3,4,1,2};
    vector<unsigned int> m341{3,4,1};
    vector<unsigned int> m342{3,4,2};
    vector<unsigned int> m340{3,4,0};
    vector<unsigned int> m3402{3,4,0,2};
    vector<unsigned int> m3401{3,4,0,1};
    vector<unsigned int> m01{0,1};
    vector<unsigned int> m34{3,4};
    
    HyperPointSet points( dim );
    HyperPoint min(dim);
    HyperPoint max(dim);

    if(dim==7 && mode==1){
        min = HyperPoint( pdg.minmax(m012).first, pdg.minmax(m02).first, pdg.minmax(m12).first, pdg.minmax(m3412).first, pdg.minmax(m341).first, pdg.minmax(m340).first, pdg.minmax(m3402).first);
        max = HyperPoint( pdg.minmax(m012).second, pdg.minmax(m02).second, pdg.minmax(m12).second, pdg.minmax(m3412).second, pdg.minmax(m341).second, pdg.minmax(m340).second, pdg.minmax(m3402).second );
    }
    else if(dim==7){
        min = HyperPoint( pdg.minmax(m012).first, pdg.minmax(m02).first, pdg.minmax(m12).first, pdg.minmax(m3412).first, pdg.minmax(m341).first, -1, -3.141 );
        max = HyperPoint( pdg.minmax(m012).second, pdg.minmax(m02).second, pdg.minmax(m12).second, pdg.minmax(m3412).second, pdg.minmax(m341).second, 1, 3.141 );
    }
    else if(dim==5 && mode ==0){
        min = HyperPoint( pdg.minmax(m012).first, pdg.minmax(m02).first, pdg.minmax(m12).first, pdg.minmax(m3412).first, pdg.minmax(m341).first );
        max = HyperPoint( pdg.minmax(m012).second, pdg.minmax(m02).second, pdg.minmax(m12).second, pdg.minmax(m3412).second, pdg.minmax(m341).second );
    }
    else if(dim==5){
        min = HyperPoint( pow(pdg.minmax(m012).first,2), pow(pdg.minmax(m02).first,2), pow(pdg.minmax(m12).first,2), pow(pdg.minmax(m3412).first,2), pow(pdg.minmax(m341).first,2) );
        max = HyperPoint( pow(pdg.minmax(m012).second,2), pow(pdg.minmax(m02).second,2), pow(pdg.minmax(m12).second,2), pow(pdg.minmax(m3412).second,2), pow(pdg.minmax(m341).second,2) );
    }
    
    HyperCuboid limits(min, max );
    
    auto invertCut_plot = NamedParameter<bool>("invertCut_plot"+to_string(filter_plotN), 0,"Invert cut logic");
    auto cut_dim_plot = NamedParameter<unsigned int>("cut_dim_plot"+to_string(filter_plotN), std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot = NamedParameter<double>("cut_limits_plot"+to_string(filter_plotN), std::vector<double>(),"cut window" ).getVector();
    phsp_cut filter_plot(cut_dim_plot,cut_limits_plot,invertCut_plot);
    
    for( auto& evt : dataEvents ){
        
        if(filter_plotN>0)if(!filter_plot(evt))continue;
        
        HyperPoint point( dim );
        if(mode == 1 && dim==5){
            point.at(0)= evt.s(m012);
            point.at(1)= evt.s(m02);
            point.at(2)= evt.s(m12);
            point.at(3)= evt.s(m3412);
            point.at(4)= evt.s(m341);
        }
        else{
            point.at(0)= sqrt(evt.s(m012));
            point.at(1)= sqrt(evt.s(m02));
            point.at(2)= sqrt(evt.s(m12));
            point.at(3)= sqrt(evt.s(m3412));
            point.at(4)= sqrt(evt.s(m341));
        }
        if(dim==7 && mode==0){
            point.at(5)= cosThetaMuAngle(evt);
            point.at(6)= chiMuAngle(evt);
        }
        else if(dim==7 && mode==1){
            point.at(5)= sqrt(evt.s(m340));
            point.at(6)= sqrt(evt.s(m3402));
        }
        point.addWeight(evt.weight());
        points.push_back(point);
    }

    HyperPointSet pointsMC( dim);
    for( auto& evt : mcEvents ){
        
        if(filter_plotN>0)if(!filter_plot(evt))continue;
        
        HyperPoint point( dim );
        if(mode == 1 && dim==5){
            point.at(0)= evt.s(m012);
            point.at(1)= evt.s(m02);
            point.at(2)= evt.s(m12);
            point.at(3)= evt.s(m3412);
            point.at(4)= evt.s(m341);
        }
        else{
            point.at(0)= sqrt(evt.s(m012));
            point.at(1)= sqrt(evt.s(m02));
            point.at(2)= sqrt(evt.s(m12));
            point.at(3)= sqrt(evt.s(m3412));
            point.at(4)= sqrt(evt.s(m341));
        }
        if(dim==7 && mode==0){
            point.at(5)= cosThetaMuAngle(evt);
            point.at(6)= chiMuAngle(evt);
        }
        else if(dim==7 && mode==1){
            point.at(5)= sqrt(evt.s(m340));
            point.at(6)= sqrt(evt.s(m3402));
        }
        point.addWeight( ( fcn( evt ) + fcn_bkg( evt ) ) * evt.weight() / evt.genPdf());
        pointsMC.push_back(point);
    }

    HyperHistogram dataHist(limits, points,
                            /*** Name of the binning algorithm you want to use     */
                            HyperBinningAlgorithms::SMART_MULTI,
                            /***  The minimum number of events allowed in each bin */
                            /***  from the HyperPointSet provided (points1)        */
                            AlgOption::MinBinContent      (minEventsPerBin),
                            /*** This minimum bin width allowed. Can also pass a   */
                            /*** HyperPoint if you would like different min bin    */
                            /*** widths for each dimension                         */
                            AlgOption::MinBinWidth        (minBinWidth),
                            /*** If you want to use the sum of weights rather than */
                            /*** the number of events, set this to true.           */
                            AlgOption::UseWeights         (true),
                            /*** Some algorithms use a random number generator. Set*/
                            /*** the seed here                                     */
                            AlgOption::RandomSeed         (1),
                            /*** What dimesnion would you like to split first? Only*/
                            /*** applies to certain algortihms                     */
                            AlgOption::StartDimension     (0)
                            /*** What dimesnions would you like to bin in?         */
                            //AlgOption::BinningDimensions  (binningDims),
                            /*** Setting this option will make the agorithm draw   */
                            /*** the binning scheme at each iteration              */
                            //AlgOption::DrawAlgorithm("Algorithm")
                            );
    
    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC);
    mcHist.normalise(dataHist.integral());
    
    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();
    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    
    vector<double> vals;
    vals.push_back(chi2);
    vals.push_back((double)nBins);
    
    return vals;
}

int main( int argc, char* argv[])
{
  time_t startTime = time(0);
  OptionsParser::setArgs( argc, argv );

  std::string dataFile = NamedParameter<std::string>("DataSample","", "Name of file containing data sample to fit." );
  std::string weightData = NamedParameter<std::string>("weightData", "");  
  if(weightData==" " || weightData=="noWeight")  weightData="";
  std::string intFile = NamedParameter<std::string>("IntegrationSample","","Name of file containing events to use for MC integration.");
  std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
  std::string phspFile = NamedParameter<std::string>("phspFile","","phspFile");

  std::string outDir = NamedParameter<std::string>("outDir", ".");
  std::string logFile = NamedParameter<std::string>("LogFile", "log.txt", "Name of the output log file");
  std::string tableFile = NamedParameter<std::string>("TableFile", "table.tex", "Name of the output log file");
  std::string modelFile = NamedParameter<std::string>("ModelFile", "model.txt", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("ResultsFile", "result.root", "Name of the output plot file");
    
  auto MinEventsChi2 = NamedParameter<Int_t>("MinEventsChi2", 15, "MinEventsChi2" );
  auto chi2Mode = NamedParameter<Int_t>("chi2Mode", 0);
    
  std::string doSystematic = NamedParameter<std::string>("doSystematic", "Baseline");  
  cout << "Doing systematic " << doSystematic << endl;  
  
  auto normAmps = NamedParameter<bool>("normAmps", 1);
  auto excludeNorm = NamedParameter<std::string>("excludeNorm", std::vector<std::string>()).getVector();
  auto combineNorm = NamedParameter<std::string>("combineNorm", std::vector<std::string>()).getVector();

  int status = 0;
  if(outDir != ".") status &= system( ("mkdir -p " + outDir).c_str() );
  if( status != 0 ) ERROR("Building OutDir directory failed");      
  logFile = outDir + "/" + logFile;
  tableFile = outDir + "/" + tableFile;
  modelFile = outDir + "/" + modelFile;
  plotFile = outDir + "/" + plotFile;
    
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()).getVector();
  auto bExtraNames = NamedParameter<std::string>("ExtraBranches", std::vector<std::string>()).getVector();
  auto bNamesMC = NamedParameter<std::string>("BranchesMC", std::vector<std::string>()).getVector();
  auto bExtraNamesMC = NamedParameter<std::string>("ExtraBranchesMC", std::vector<std::string>()).getVector();
  //if(bNamesMC.size()==0)bNamesMC=bNames;
  //if(bExtraNamesMC.size()==0)bExtraNamesMC=bExtraNames;
  auto pNames = NamedParameter<std::string>("EventType" , "").getVector(); 
  
  int      nThreads = NamedParameter<int>     ("nCores"    , 8           , "Number of threads to use" );
  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
  auto        nBins    = NamedParameter<Int_t>     ("nBins"      , 50           , "Number of bins" );

  if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

  TRandom3 rndm;
  seed =  atoi(argv[1]);
  cout << "Using random seed = " << seed << endl;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("LogFile: " << logFile << "; Plots: " << plotFile );
  
#ifdef _OPENMP
  nThreads = nThreads < 0 ? min(omp_get_num_procs(),(int)thread::hardware_concurrency()) : nThreads;  
  INFO( "Hardware_concurrency " << thread::hardware_concurrency());
  INFO( "omp_get_num_procs " << omp_get_num_procs() );
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  MinuitParameterSet MPS;
  MPS.loadFromStream();

  // Add amps from list
  auto addAmpList = NamedParameter<string>("addAmpList","");
  if(addAmpList != ""){
        MinuitParameterSet* addAmp = new MinuitParameterSet();
        addAmp->loadFromFile(addAmpList);
        INFO("Found " << addAmp->size()/2 << " amplitudes in " << addAmpList);
        int amp_index = atoi(argv[1])*2;
        if(amp_index < addAmp->size()-1){ 
            MinuitParameter* amp_re= addAmp->at(amp_index);
            MinuitParameter* amp_im= addAmp->at( replaceAll(addAmp->at(amp_index)->name(),"_Re","_Im") );
            if(amp_re == 0 || amp_im == 0 ){
                ERROR("Could not add amplitude: " << replaceAll(addAmp->at(amp_index)->name(),"_Re","") );
                return 0;
            }
            MPS.addOrGet( amp_re->name(), amp_re->flag(), amp_re->meanInit(), amp_re->stepInit(), amp_re->minInit(), amp_re->maxInit()  );
            MPS.addOrGet( amp_im->name(), amp_im->flag(), amp_im->meanInit(), amp_im->stepInit(), amp_im->minInit(), amp_im->maxInit()  );
            INFO("Added amplitude: " << replaceAll(addAmp->at(amp_index)->name(),"_Re","") );
        }
        else ERROR("Not enough amplitudes in list");
        delete addAmp;
  }
    
  //Remove amps from list
  vector<string> removeAmps = NamedParameter<string>("removeAmps",vector<string>() ).getVector();
  for(int n=0; n<removeAmps.size();n++){
      for(int i=0;i<MPS.size();i++){
          if(MPS[i]->name().find(removeAmps[n]) != std::string::npos){
              INFO("Removing parameter " << MPS[i]->name());
              MPS.unregister(MPS[i]);
              i=-1;
          }
      }
  }
    
  //Randomize start values   
  auto randomizeStartVals = NamedParameter<int>("randomizeStartVals", 0);
  if(randomizeStartVals==1) randomizeStartingPoint(MPS,rndm);
  if(randomizeStartVals>1) perturb(MPS,randomizeStartVals-1);

  //Scan parameter
  auto scanParam = NamedParameter<std::string>("scanParam", std::vector<std::string>()).getVector(); // name, min, max, nSteps
  if(scanParam.size()==4)scan( MPS, scanParam[0], std::stod(scanParam[1]), std::stod(scanParam[2]), seed, std::stoi(scanParam[3]) );
        
  auto scale_transform = [bExtraNames](auto& event){ for( size_t x = 0 ; x < event.size()-bExtraNames.size(); ++x ) event[x] /= 1000.; };
  auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
  auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
  auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
  auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();
  auto invertCut2 = NamedParameter<bool>("invertCut2", 0,"Invert cut logic");
  auto cut_dim2 = NamedParameter<unsigned int>("cut_dim2", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
  auto cut_limits2 = NamedParameter<double>("cut_limits2", std::vector<double>(),"cut window" ).getVector();
  auto invertCut3 = NamedParameter<bool>("invertCut3", 0,"Invert cut logic");
  auto cut_dim3 = NamedParameter<unsigned int>("cut_dim3", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
  auto cut_limits3 = NamedParameter<double>("cut_limits3", std::vector<double>(),"cut window" ).getVector();
  phsp_cut filter(cut_dim,cut_limits,invertCut);
  phsp_cut filter2(cut_dim2,cut_limits2,invertCut2);
  phsp_cut filter3(cut_dim3,cut_limits3,invertCut3);

  EventType evtType(pNames);
  vector<size_t> entryList,entryListMC,entryListPlotMC;
 
  auto doBootstrap = NamedParameter<int>("doBootstrap", 0);
  auto N_bootstrap = NamedParameter<int>("N_bootstrap", -1);
  auto maxDataEvents = NamedParameter<int>("maxDataEvents", -1);

  auto doBootstrapMC = NamedParameter<int>("doBootstrapMC", 0);
  auto N_bootstrapMC = NamedParameter<int>("N_bootstrapMC", -1);
  auto maxIntEvents = NamedParameter<int>("maxIntEvents", -1);
  
  if(useFilter==1){
        EventList_type dummyEvents(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));
        EventList_type dummyEventsMC(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
        
        if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") dummyEvents.transform( scale_transform );
        if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") dummyEventsMC.transform( scale_transform );
        
        for(int i=0; i< dummyEvents.size(); i++ )if(filter(dummyEvents[i]) && filter2(dummyEvents[i]) && filter3(dummyEvents[i])){
            if(entryList.size()<maxDataEvents)entryList.push_back(i);
        }
        for(int i=0; i< dummyEventsMC.size(); i++ )if(filter(dummyEventsMC[i])&& filter2(dummyEventsMC[i]) && filter3(dummyEventsMC[i])){
            if(entryListMC.size()<maxIntEvents)entryListMC.push_back(i);
            entryListPlotMC.push_back(i);
        }
  }
  else{
      if(doBootstrap){
          EventList_type dummyEvents(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));
          int N_sample = dummyEvents.size();
          if(N_bootstrap == -1)N_bootstrap = N_sample;
          
          vector<int> b_indices;
          while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(rndm.Uniform(0,N_sample-1)));
          sort(b_indices.begin(), b_indices.end());
          for(unsigned int i=0; i<b_indices.size(); i++)entryList.push_back(b_indices[i]);
      }
      
      if(doBootstrapMC){
          EventList_type dummyEvents(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
          int N_sample = dummyEvents.size();
          if(N_bootstrap == -1)N_bootstrap = maxIntEvents > 0 ? maxIntEvents : N_sample;
          
          vector<int> b_indices;
          while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(rndm.Uniform(0,N_sample-1)));
          sort(b_indices.begin(), b_indices.end());
          for(unsigned int i=0; i<b_indices.size(); i++)entryListMC.push_back(b_indices[i]);
      }
  }
    
  if(entryList.size()==0 && maxDataEvents>0)for(int i = 0; i < maxDataEvents; i++)entryList.push_back(i);
  if(entryListMC.size()==0 && maxIntEvents>0)for(int i = 0; i < maxIntEvents; i++)entryListMC.push_back(i);

  EventList_type events(dataFile, evtType, Branches(bNames), ExtraBranches(bExtraNames), GetGenPdf(false), WeightBranch(weightData),EntryList(entryList));
  EventList_type eventsMC(intFile, evtType, Branches(bNamesMC), ExtraBranches(bExtraNamesMC), WeightBranch(weightMC), GetGenPdf(true),EntryList(entryListMC));
  double nSig = events.sumWeights();

  if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
    INFO("Changing data units from MeV -> GeV");
    events.transform( scale_transform );
  }
  if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
    INFO("Changing MC units from MeV -> GeV");
    eventsMC.transform( scale_transform );
  }

  auto pdfType = NamedParameter<pdfTypes>( "Type", pdfTypes::PolarisedSum);
  if(pdfType==pdfTypes::PolarisedSum){
        
          vector<string> replaceAmpAwithB = NamedParameter<string>("replaceAmpAwithB",vector<string>({""}) ).getVector();
          vector<FitResult*> replaceAmpAwithBResults;
          for(int n=0; n<replaceAmpAwithB.size();n++){
          
              // Signal pdf
              PolarisedSum sig(evtType, MPS);
              sig.setMC( eventsMC );
              checkAmps(sig, MPS);
              
              // Bkg pdf
              IncoherentSum bkg( evtType, MPS, "Inco" );
              bkg.setMC( eventsMC );
              
              // Bkg pdf from BDT
              BkgPDF bkgBDT( evtType );
              bkgBDT.setMC( eventsMC );
              
              sig.setWeight( MPS["f_sig"] );
              bkg.setWeight( MPS["f_bkg"] );
              bkgBDT.setWeight( MPS["f_bkg"] );
              
              sanityChecks(MPS);
              cout << "Number of amplitudes = " << sig.numAmps() << endl;
              
              // Set spline start vals
              string head = NamedParameter<std::string>("Head","");
              if(head!="")setSplineVals(MPS,head);
              
              if(doSystematic=="Res"){
                  std::vector<std::string> paramsToVary = NamedParameter<std::string>( "ParamsToVary",std::vector<std::string>() ).getVector();
                  int paramIndex = seed < paramsToVary.size() ? seed : seed - paramsToVary.size() ;
                  if(paramIndex >= paramsToVary.size()) FATAL("Nothing to vary");
                  
                  INFO("Have " << paramsToVary.size() << " params to vary +/- 1 sigma " );
                  vary(MPS,paramsToVary[ paramIndex ], seed < paramsToVary.size() ? -1 : +1);
                  doSystematic += paramsToVary[ paramIndex ] + seed < paramsToVary.size() ? (string) "_m" : (string) "_p"  + "1sigma" ;
              }
              
              if(doSystematic=="ResAll"){
                  std::vector<std::string> paramsToVary = NamedParameter<std::string>( "ParamsToVary",std::vector<std::string>() ).getVector();
                  INFO("Have " << paramsToVary.size() << " params to vary within uncertainties" );
                  varyGauss(MPS, paramsToVary );
              }
              
              unsigned int useBkgBDT = NamedParameter<unsigned int>("useBkgBDT",1);
              if(useBkgBDT){
                  if(bExtraNames.size()==0) ERROR("No bkgPdf branch set for data");
                  if(bExtraNamesMC.size()==0) ERROR("No bkgPdf branch set for MC");
              }
              auto ll = make_likelihood(events, sig, bkg) ;
              auto ll_bdt = make_likelihood(events, sig, bkgBDT) ;
              sig.prepare();
              bkg.prepare();
              bkgBDT.prepare();
              
              if(normAmps){
                  if(phspFile == ""){
                      sig.normaliseAmps(excludeNorm);
                  }
                  else{
                      auto bNamesPhsp = NamedParameter<std::string>("BranchesPhsp", std::vector<std::string>()).getVector();
                      std::string weightPhsp = NamedParameter<std::string>("weightPhsp", "weight");
                      EventList_type eventsPhspMC = EventList_type(phspFile, evtType, GetGenPdf(true), Branches(bNamesPhsp), WeightBranch(weightPhsp));
                      sig.setMC( eventsPhspMC );
                      sig.prepare();
                      sig.normaliseAmps(excludeNorm);
                      sig.setMC( eventsMC );
                      sig.prepare();
                  }
              }
              
              if(useBkgBDT){
                  events[0].print();
                  cout << "sig.getValNoCache(evt) = " <<  sig.getValNoCache(events[0]) << endl;
                  cout << "bkgBDT.prob_unnormalisedNoCache(events[0]) = " <<  bkgBDT.prob_unnormalisedNoCache(events[0]) << endl;
                  eventsMC[0].print();
                  cout << "sig.getValNoCache(evt) = " <<  sig.getValNoCache(eventsMC[0]) << endl;
                  cout << "bkgBDT.prob_unnormalisedNoCache(eventsMC[0]) = " <<  bkgBDT.prob_unnormalisedNoCache(eventsMC[0]) << endl;
                  
                  cout << "bkgBDT.norm() = " << bkgBDT.norm() << endl;
                  cout << "bkg.norm() = " << bkg.norm() << endl;
                  cout << "sig.norm() = " << sig.norm() << endl;
                  
                  double norm_sum(0.),norm_sumBkg(0.), norm_sumBkgBDT(0.);
                  for(auto& evt : events){
                      norm_sumBkgBDT += bkgBDT(evt);
                      norm_sumBkg += bkg.prob_unnormalisedNoCache(evt) * bkg.getWeight()/bkg.norm();
                      norm_sum += sig(evt);
                      if(TMath::IsNaN(sig.getValNoCache(evt)))evt.print();
                  }
                  cout << "norm_sumBDT = " << norm_sumBkgBDT << endl;
                  cout << "norm_sumBkg = " << norm_sumBkg << endl;
                  cout << "norm_sum = " << norm_sum << endl;
                  //throw "";
              }
              
              int badEventsData = 0;
              for(auto& evt : events){
                  if(TMath::IsNaN(sig.getValNoCache(evt))){
                      //evt.print();
                      badEventsData++;
                  }
              }
              WARNING("badEventsData = " << badEventsData);
              
              int badEventsMC = 0;
              for(auto& evt : eventsMC){
                  if(TMath::IsNaN(sig.getValNoCache(evt))){
                      //evt.print();
                      badEventsMC++;
                  }
              }
              WARNING("badEventsMC = " << badEventsMC);
              
              // Do fit
              auto nFits = NamedParameter<int>("nFits", 1);
              auto performFit = NamedParameter<bool>("doFit", 1,"doFit");
              auto performLASSO = NamedParameter<bool>("doLASSO", 0,"doLASSO");
              double lambda = performLASSO ? seed + 0.1  : -1;
              
              if(!performFit){
                  for(int i=0;i<MPS.size();i++){
                      MPS[i]->fix();
                  }
              }
              
              double min_LL = 9999999;
              double minChi2 = 999;
              double minChi2PerBin = 999;
              double minChi2Dalitz = 999;
              int nParams = 0;
              
              FitResult* fr_best;
              for (unsigned int i = 0; i< nFits; ++i) {
                  
                  cout << "==============================================" << endl;
                  INFO("Start fit number " << i);
                  
                  if(i>0){
                      if(scanParam.size()==4)perturb(MPS,i);
                      else randomizeStartingPoint(MPS,rndm);
                      
                      sig.setMC( eventsMC );
                      sig.prepare();
                      //sig.normaliseAmps(excludeNorm);
                  }
                  FitResult* fr;
                  if(useBkgBDT) fr = doFit(ll_bdt, events, eventsMC, MPS, sig, lambda );
                  else fr = doFit(ll, events, eventsMC, MPS, sig, lambda );
                  
                  if(fr->LL()>min_LL || TMath::IsNaN(fr->LL())){
                      INFO("Fit did not improve: LL = " << fr->LL() << " ; min_LL = " << min_LL);
                      continue;
                  }
                  else {
                      INFO("Fit did improve: LL = " << fr->LL() << " ; min_LL = " << min_LL);
                      min_LL=fr->LL();
                      fr_best = fr;
                  }
                  
                  auto ep = fr->getErrorPropagator();
                  auto fitFractions = sig.fitFractions( ep );
                  fr->addFractions( fitFractions );
                  auto interferenceFractions = sig.interferenceFractions( ep );
                  fr->addInterferenceFractions( interferenceFractions );
                  
                  vector<double> thresholds{0.01,0.05,0.1,0.5,1,2};
                  vector<double> numFracAboveThresholds = sig.numFracAboveThreshold(thresholds);
                  
                  // Plot spline
                  if(head!="")fr->plotSpline(head,outDir);
                  
                  // Estimate the chi2
                  auto evaluator_sig = sig.evaluator();
                  auto evaluator_bkg = useBkgBDT ? bkgBDT.evaluator() : bkg.evaluator();
                  
                  INFO("Calculating chi2 ...");
                  vector<double> chi2Hyper = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 0 );
                  INFO("Calculating chi2 Dalitz ...");
                  vector<double> chi2HyperDalitz = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 1 );

                  auto filter_plotN = NamedParameter<unsigned int>("filter_plotN", std::vector<unsigned int>()).getVector();
                  if(filter_plotN.size()>0){
                      INFO("Calculating chi2 in phsp slices ...");
                      for(auto& cutN : filter_plotN){
                          INFO("phsp slice " << to_string(cutN));
                          vector<double> chi2HyperN = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 0, cutN );
                      }
                  }

                  minChi2 = chi2Hyper[0]/(chi2Hyper[1]-1.-(double)nParams);
                  minChi2PerBin = chi2Hyper[0]/(chi2Hyper[1]-1.);
                  minChi2Dalitz = chi2HyperDalitz[0]/(chi2HyperDalitz[1]-1.-(double)nParams);
                  fr->addChi2(chi2Hyper[0],chi2Hyper[1]);
                  
                  nParams = fr->floating().size();
                  cout << "Number of amplitudes = " << sig.numAmps() << endl;
                  cout << "Number of parameters = " << nParams << endl;
                  
                  int fixParamsOptionsFile = NamedParameter<int>("fixParamsOptionsFile",0);
                  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
                  fr->setSystematic(doSystematic);
                  fr->print();
                  fr->printToLatexTable(tableFile);
                  fr->writeToRootFile( output, seed, 0, fr->LL(),  sig.numAmps(), nSig, thresholds, numFracAboveThresholds );
                  if(replaceAmpAwithB.size()>0 && n>0){
                      fr->writeToFile(replaceAll(logFile,".txt","_" + to_string(n-1) + ".txt"));
                      fr->writeToOptionsFile(replaceAll(modelFile,".txt","_" +  to_string(n-1) + ".txt"), 0);
                  }
                  else{
                      fr->writeToFile(logFile);
                      fr->writeToOptionsFile(modelFile, fixParamsOptionsFile);
                  }
                  output->cd();
                  output->Close();
                  
                  unsigned int saveWeights   = NamedParameter<unsigned int>("saveWeights",1);
                  if( saveWeights ){
                      EventList_type eventsPlotMC;
                      if(maxIntEvents == -1) eventsPlotMC = eventsMC;
                      else {
                          eventsPlotMC = EventList_type(intFile, evtType, Branches(bNamesMC), ExtraBranches(bExtraNamesMC), WeightBranch(weightMC), GetGenPdf(true), EntryList(entryListPlotMC));
                      }
                      sig.setMC( eventsPlotMC );
                      sig.prepare();
                      if(useBkgBDT){
                          bkgBDT.setMC( eventsPlotMC );
                          bkgBDT.prepare();
                          makePlotWeightFile(sig, bkgBDT, eventsPlotMC, ep, rndm);
                      }
                      else{
                          bkg.setMC( eventsPlotMC );
                          bkg.prepare();
                          makePlotWeightFile(sig, bkg, eventsPlotMC, ep, rndm);
                      }
                      //Chi2Estimator chi2Plot( events, eventsPlotMC, evaluator, MinEvents(MinEventsChi2), Dim(3*(pNames.size()-1)-7));
                      //fr->addChi2( chi2Plot.chi2(), chi2Plot.nBins() );
                      //fr->print();
                  }
                  
              }
              
              replaceAmpAwithBResults.push_back(fr_best);
              
              vector<string> paramsToRelease = NamedParameter<string>("ParamsToRelease",vector<string>() ).getVector();
              vector<FitResult*> paramsToReleaseResults;
              vector<int> paramsToReleaseNParams;
              for(auto& param: paramsToRelease){
                  auto it = MPS.find( param );
                  if(it == nullptr) continue;
                  it->setFree();
                  auto it_width = MPS.find( replaceAll(it->name(),"_mass","_width") );
                  if(it_width != nullptr){
                      it_width->setFree();
                      INFO("Refit with released " << replaceAll(it->name(),"_mass","") << " mass and width" );
                  }
                  else INFO("Refit with released " << it->name() );
                  FitResult* fr;
                  if(useBkgBDT) fr = doFit(ll_bdt, events, eventsMC, MPS, sig );
                  else fr = doFit(ll, events, eventsMC, MPS, sig );
                  
                  auto evaluator_sig = sig.evaluator();
                  auto evaluator_bkg = useBkgBDT ? bkgBDT.evaluator() : bkg.evaluator();
                  vector<double> chi2Hyper = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 0 );
                  vector<double> chi2HyperDalitz = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 1 );
                  
                  int dParams = abs((int)fr->floating().size() - nParams);
                  double LL_sig = sqrt(2) * TMath::ErfcInverse(TMath::Prob( abs(fr->LL() - min_LL),dParams ));
                  if(abs(fr->LL() - min_LL)>0 && TMath::Prob(abs(fr->LL() - min_LL),dParams) == 0 ) LL_sig = 100;
                  if(fr->LL() > min_LL) LL_sig *= -1.;
                  INFO( "LL = " << fr->LL() );
                  INFO( "LL improvement = " << fr->LL() - min_LL << " (" << LL_sig << " sigma)" );
                  INFO( "Chi2 improvement = " << minChi2 - chi2Hyper[0]/(chi2Hyper[1]-1.-(double)fr->floating().size()) );
                  INFO( "Chi2PerBin improvement = " << minChi2PerBin - chi2Hyper[0]/(chi2Hyper[1]-1.) );
                  INFO( "Chi2Dalitz improvement = " << minChi2Dalitz - chi2HyperDalitz[0]/(chi2HyperDalitz[1]-1.-(double)fr->floating().size()) );
                  
                  fr->addChi2(chi2Hyper[0],chi2Hyper[1]);
                  fr->print();
                  paramsToReleaseResults.push_back(fr);
                  paramsToReleaseNParams.push_back((int)fr->floating().size());
                  
                  it->resetToInit();
                  it->fix();
                  if(it_width != nullptr){
                      it_width->resetToInit();
                      it_width->fix();
                  }
              }
              
              if(paramsToRelease.size()>0){
                  INFO("----SUMMARY::ParamsToRelease----");
                  for(int i=0;i<paramsToRelease.size();i++){
                      FitResult* fr = paramsToReleaseResults[i];
                      int dParams = abs(paramsToReleaseNParams[i] - nParams);
                      double LL_sig = sqrt(2) * TMath::ErfcInverse(TMath::Prob(abs(fr->LL() - min_LL),dParams));
                      if(abs(fr->LL() - min_LL)>0 && TMath::Prob(abs(fr->LL() - min_LL),dParams) == 0 ) LL_sig = 100;
                      if(fr->LL() > min_LL) LL_sig *= -1.;
                      double dChi2 = minChi2 - fr->chi2()/((double)fr->nBins()-1.-(double)paramsToReleaseNParams[i]);
                      double dChi2PerBin = minChi2PerBin - fr->chi2()/((double)fr->nBins()-1.);
                      INFO("Param= " << paramsToRelease[i] << ":: dLL = " << fr->LL() - min_LL << " (" << LL_sig << " sigma) ; dChi2 = " << dChi2 << " ; dChi2PerBin = " << dChi2PerBin );
                  }
              }
              
              vector<string> ampsToRemove = NamedParameter<string>("AmpsToRemove",vector<string>() ).getVector();
              vector<FitResult*> ampsToRemoveResults;
              vector<int> ampsToRemoveNParams;
              int counterAmpsToRemove = 0;
              int scale_dParams = 1;

              for(auto& amp: ampsToRemove){
                  
                  vector<MinuitParameter*> removedParams;
                  vector<bool> removedParamIsFixed;
                  
                  for(int i=0;i<MPS.size();i++){
                      auto param = MPS[i];
                      if(param->name().find(amp) != std::string::npos)
                          if(param->name().find("_Re") != std::string::npos || param->name().find("_Im") != std::string::npos
                             || param->name().find("_mass") != std::string::npos || param->name().find("_width") != std::string::npos
                             ){
                              removedParams.push_back(param);
                              removedParamIsFixed.push_back(param->isFixed());
                              if(param->name().find("_mass") != std::string::npos || param->name().find("_width") != std::string::npos)scale_dParams=2;
                              else param->setCurrentFitVal(0.);
                              param->fix();
                              INFO("Refit and fix parameter " << param->name() << " to 0");
                          }
                  }
                  
                  FitResult* fr;
                  if(useBkgBDT) fr = doFit(ll_bdt, events, eventsMC, MPS, sig );
                  else fr = doFit(ll, events, eventsMC, MPS, sig );
                  
                  auto evaluator_sig = sig.evaluator();
                  auto evaluator_bkg = useBkgBDT ? bkgBDT.evaluator() : bkg.evaluator();
                  vector<double> chi2Hyper = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 0 );
                  vector<double> chi2HyperDalitz = getChi2MuMu(events, eventsMC, evaluator_sig, evaluator_bkg, 7, MinEventsChi2, 0., 1 );
                  
                  int dParams = abs((int)fr->floating().size() - nParams);
                  double LL_sig = sqrt(2) * TMath::ErfcInverse(TMath::Prob(abs(fr->LL() - min_LL),scale_dParams*dParams));
                  if(abs(fr->LL() - min_LL)>0 && TMath::Prob(abs(fr->LL() - min_LL),dParams) == 0 ) LL_sig = 100;
                  if(fr->LL() > min_LL) LL_sig *= -1.;
                  INFO( "LL = " << fr->LL() );
                  INFO( "LL improvement = " << fr->LL() - min_LL << " (" << LL_sig << " sigma)" );
                  INFO( "Chi2 improvement = " << minChi2 - chi2Hyper[0]/(chi2Hyper[1]-1.-(double)nParams) );
                  INFO( "Chi2PerBin improvement = " << minChi2PerBin - chi2Hyper[0]/(chi2Hyper[1]-1.) );
                  INFO( "Chi2Dalitz improvement = " << minChi2Dalitz - chi2HyperDalitz[0]/(chi2HyperDalitz[1]-1.-(double)nParams) );
                  
                  fr->addChi2(chi2Hyper[0],chi2Hyper[1]);
                  ampsToRemoveResults.push_back(fr);
                  ampsToRemoveNParams.push_back((int)fr->floating().size());
                  
                  fr->print();
                  fr->writeToFile(replaceAll(logFile,".txt","_" + to_string(counterAmpsToRemove) + ".txt"));
                  fr->writeToOptionsFile(replaceAll(modelFile,".txt","_" +  to_string(counterAmpsToRemove) + ".txt"), 0);
                  counterAmpsToRemove++;
                  
                  for(int i=0;i<removedParams.size();i++){
                      removedParams[i]->resetToInit();
                      if(!removedParamIsFixed[i])removedParams[i]->setFree();
                  }
                  
              }
              
              if(ampsToRemove.size()>0){
                  INFO("----SUMMARY::AmpsToRemove----");
                  for(int i=0;i<ampsToRemove.size();i++){
                      FitResult* fr = ampsToRemoveResults[i];
                      int dParams = abs(ampsToRemoveNParams[i] - nParams);
                      double LL_sig = sqrt(2) * TMath::ErfcInverse(TMath::Prob(abs(fr->LL() - min_LL),scale_dParams*dParams));
                      if(abs(fr->LL() - min_LL)>0 && TMath::Prob(abs(fr->LL() - min_LL),dParams) == 0 ) LL_sig = 100;
                      if(fr->LL() > min_LL) LL_sig *= -1.;
                      double dChi2 = minChi2 - fr->chi2()/((double)fr->nBins()-1.-(double)ampsToRemoveNParams[i]);
                      double dChi2PerBin = minChi2PerBin - fr->chi2()/((double)fr->nBins()-1.);
                      INFO("Amp= " << ampsToRemove[i] << ":: dLL = " << fr->LL() - min_LL << " (" << LL_sig << " sigma) ; dChi2 = " << dChi2 << " ; dChi2PerBin = " << dChi2PerBin  );
                  }
              }
              
              //Replace amp A with B and refit
              if(replaceAmpAwithB.size()>=2 && n+1 < replaceAmpAwithB.size()){
                  for(int i=0;i<MPS.size();++i){
                      auto param = MPS[i];
                      //INFO("param " << i << " = " << param->name());
                      //if(param->name().find("B+[D]{Z(4240)A+{psi(2S)0,pi+},K*(892)0{K+,pi-}}") != std::string::npos)INFO("Found " << param->name());
                      if(param->name().find(replaceAmpAwithB[n]) != std::string::npos){
                          auto name = param->name();
                          auto new_name = replaceAll(name,replaceAmpAwithB[n],replaceAmpAwithB[n+1]);
                          
                          // fix L
                          if(new_name.find("B+[P]{Z(T)")!= std::string::npos)
                              new_name = replaceAll(new_name,"[P]","");
                          if(new_name.find("B+[D]{Z(T)")!= std::string::npos)
                              new_name = replaceAll(new_name,"[D]","");
                          if(new_name.find("B+[P]{Z(P)")!= std::string::npos)
                              new_name = replaceAll(new_name,"[P]","");
                          if(new_name.find("B+[D]{Z(P)")!= std::string::npos)
                              new_name = replaceAll(new_name,"[D]","");
                          if(new_name.find("B+[P]{rhoOmega20,Zs(P)")!= std::string::npos)
                              new_name = replaceAll(new_name,"[P]","");
                          if(new_name.find("B+[D]{rhoOmega20,Zs(P)")!= std::string::npos)
                              new_name = replaceAll(new_name,"[D]","");
                          
                          // fix particle ordering
                          if(new_name.find("Z(P)+{psi(2S)0,pi+},K*(892)0{K+,pi-}")!= std::string::npos)
                              new_name = replaceAll(new_name,"Z(P)+{psi(2S)0,pi+},K*(892)0{K+,pi-}","K*(892)0{K+,pi-},Z(P)+{psi(2S)0,pi+}");
                          if(new_name.find("Z(P)+{psi(2S)0,pi+},KPi40")!= std::string::npos)
                              new_name = replaceAll(new_name,"Z(P)+{psi(2S)0,pi+},KPi40","KPi40,Z(P)+{psi(2S)0,pi+}");

                          if(new_name.find("K*(892)0{K+,pi-},Z(T)+{psi(2S)0,pi+}")!= std::string::npos)
                              new_name = replaceAll(new_name,"K*(892)0{K+,pi-},Z(T)+{psi(2S)0,pi+}","Z(T)+{psi(2S)0,pi+},K*(892)0{K+,pi-}");
                          if(new_name.find("KPi40,Z(T)+{psi(2S)0,pi+}")!= std::string::npos)
                              new_name = replaceAll(new_name,"KPi40,Z(T)+{psi(2S)0,pi+}","Z(T)+{psi(2S)0,pi+},KPi40");

                          if(new_name.find("K*(892)0{K+,pi-},Z(PT)+{psi(2S)0,pi+}")!= std::string::npos)
                              new_name = replaceAll(new_name,"K*(892)0{K+,pi-},Z(PT)+{psi(2S)0,pi+}","Z(PT)+{psi(2S)0,pi+},K*(892)0{K+,pi-}");
                          if(new_name.find("KPi40,Z(PT)+{psi(2S)0,pi+}")!= std::string::npos)
                              new_name = replaceAll(new_name,"KPi40,Z(PT)+{psi(2S)0,pi+}","Z(PT)+{psi(2S)0,pi+},KPi40");
                          
                          if(new_name.find("Zs(P)+{psi(2S)0,K+},PiPi40")!= std::string::npos)
                              new_name = replaceAll(new_name,"Zs(P)+{psi(2S)0,K+},PiPi40","PiPi40,Zs(P)+{psi(2S)0,K+}");
                          if(new_name.find("Zs(P)+{psi(2S)0,K+},rhoOmega20")!= std::string::npos)
                              new_name = replaceAll(new_name,"Zs(P)+{psi(2S)0,K+},rhoOmega20","rhoOmega20,Zs(P)+{psi(2S)0,K+}");

                          if(new_name.find("B+{rhoOmega20,Zs(T)+{psi(2S)0,K+}}")!= std::string::npos)
                              new_name = replaceAll(new_name,"B+{rhoOmega20,Zs(T)+{psi(2S)0,K+}}","B+{Zs(T)+{psi(2S)0,K+},rhoOmega20}");
                          if(new_name.find("B+{PiPi40,Zs(T)+{psi(2S)0,K+}}")!= std::string::npos)
                              new_name = replaceAll(new_name,"B+{PiPi40,Zs(T)+{psi(2S)0,K+}}","B+{Zs(T)+{psi(2S)0,K+},PiPi40}");

                          if(new_name.find("B+{rhoOmega20,Zs(PT)+{psi(2S)0,K+}}")!= std::string::npos)
                              new_name = replaceAll(new_name,"B+{rhoOmega20,Zs(PT)+{psi(2S)0,K+}}","B+{Zs(PT)+{psi(2S)0,K+},rhoOmega20}");
                          if(new_name.find("B+{PiPi40,Zs(PT)+{psi(2S)0,K+}}")!= std::string::npos)
                              new_name = replaceAll(new_name,"B+{PiPi40,Zs(PT)+{psi(2S)0,K+}}","B+{Zs(PT)+{psi(2S)0,K+},PiPi40}");
                          
                          if(MPS.find(new_name) != nullptr)MPS.unregister(MPS.find(new_name));
                          
                          //check for not allowed decays
                          if(new_name.find("X(P)")!= std::string::npos && (new_name.find("Z")!= std::string::npos || new_name.find("PiPi")!= std::string::npos) )WARNING(new_name << " not allowed");
                          else if(new_name.find("Xs(P)")!= std::string::npos && (new_name.find("Z")!= std::string::npos || new_name.find("KPi")!= std::string::npos) )WARNING(new_name << " not allowed");
                          else if(new_name.find("X(S)")!= std::string::npos && (new_name.find("Z(4055)V")!= std::string::npos || new_name.find("Z(4100)")!= std::string::npos || new_name.find("Z(V)")!= std::string::npos ) )WARNING(new_name << " not allowed");
                          else if(new_name.find("X(P)")!= std::string::npos && ( new_name.find("Z(4055)A")!= std::string::npos || new_name.find("Z(A)")!= std::string::npos || new_name.find("Z(4240)A")!= std::string::npos || new_name.find("Z(4430)")!= std::string::npos ) )WARNING(new_name << " not allowed");
                          
                          else if(new_name.find("X(V)")!= std::string::npos && new_name.find("Z(P)")!= std::string::npos  )WARNING(new_name << " not allowed");
                          else if(new_name.find("X(4685)")!= std::string::npos && new_name.find("Z(P)")!= std::string::npos  )WARNING(new_name << " not allowed");

                          else if( (new_name.find("X(4500)")!= std::string::npos || new_name.find("X(4700)") != std::string::npos || new_name.find("X(S)") != std::string::npos) && new_name.find("Z(V)")!= std::string::npos )WARNING(new_name << " not allowed");
                          
                          else if( (new_name.find("Xs(A)")!= std::string::npos || new_name.find("Xs2(A)") != std::string::npos) && new_name.find("Zs(V)")!= std::string::npos )WARNING(new_name << " not allowed");
                          else if( (new_name.find("Xs(A)")!= std::string::npos || new_name.find("Xs2(A)") != std::string::npos) && new_name.find("Z(V)")!= std::string::npos )WARNING(new_name << " not allowed");
                          else if( (new_name.find("Xs(A)")!= std::string::npos || new_name.find("Xs2(A)") != std::string::npos) && new_name.find("Zs(P)")!= std::string::npos )WARNING(new_name << " not allowed");
                          else if( (new_name.find("Xs(A)")!= std::string::npos || new_name.find("Xs2(A)") != std::string::npos) && new_name.find("Z(P)")!= std::string::npos )WARNING(new_name << " not allowed");
                          
                          else if( (new_name.find("Xs3(P)")!= std::string::npos || new_name.find("Xs3(PT)") != std::string::npos)
                                  && ( new_name.find("Z(4240)A")!= std::string::npos || new_name.find("Z(4430)")!= std::string::npos || new_name.find("Zs(4000)")!= std::string::npos || new_name.find("Zs(4220)A")!= std::string::npos  ) )WARNING(new_name << " not allowed");

                          else{
                              MPS.add(new_name,param->flag(),param->meanInit(),param->stepInit(),param->minInit(),param->maxInit());
                              INFO("Replaced param " << name << " with " << new_name);
                          }
                          MPS.unregister(param);
                          i=-1;
                      }
                  }
              }
              
          }
      
          if(replaceAmpAwithB.size()>=2){
              INFO("----SUMMARY::replaceAmpAwithB----");
              FitResult* fr_base = replaceAmpAwithBResults[0];
              INFO("Replaced amp = " << replaceAmpAwithB[0]);
              fr_base->print();
              for(int i=1;i<replaceAmpAwithB.size();i++){
                  FitResult* fr = replaceAmpAwithBResults[i];
                  double dLL = fr->LL() - fr_base->LL();
                  double LL_sig = sqrt(2) * TMath::ErfcInverse(TMath::Prob(abs(dLL),1));
                  if(abs(dLL)>0 && TMath::Prob(abs(dLL),1) == 0 ) LL_sig = 100;
                  if(dLL>0) LL_sig *= -1.;
                  double dChi2 = fr_base->chi2()/((double)fr_base->nBins()-1.-fr_base->nParam()) - fr->chi2()/((double)fr->nBins()-1.-fr->nParam());
                  double dChi2PerBin = fr_base->chi2()/((double)fr_base->nBins()-1.) - fr->chi2()/((double)fr->nBins()-1.);
                  INFO("Amp= " << replaceAmpAwithB[i] << ":: dLL = " << dLL << " (" << LL_sig << " sigma) ; dChi2 = " << dChi2 << " ; dChi2PerBin = " << dChi2PerBin  );
              }
          }
  //
  }

  else if(pdfType==pdfTypes::IncoherentSum){
        
        // Signal pdf
        IncoherentSum sig(evtType, MPS, "Inco");
        sig.setMC( eventsMC );
        //checkAmps(sig, MPS);
        sanityChecks(MPS);
        //cout << "Number of amplitudes = " << sig.numAmps() << endl;
            
        auto ll = make_likelihood(events, sig);
        sig.prepare();
        //sig.normaliseAmps();
            
        // Do fit
        auto nFits = NamedParameter<int>("nFits", 1);
        double min_LL = 99999;
    
        for (unsigned int i = 0; i< nFits; ++i) {

            cout << "==============================================" << endl;
            INFO("Start fit number " << i);

            if(i>0){
                randomizeStartingPoint(MPS,rndm);
                sig.setMC( eventsMC );
                //sig.prepare();
            }
            FitResult* fr;// = doFit(ll, events, eventsMC, MPS );
            
            if(fr->LL()>min_LL || TMath::IsNaN(fr->LL())){
                INFO("Fit did not improve: LL = " << fr->LL() << " ; min_LL = " << min_LL);
                //continue;
            }
            else {
                INFO("Fit did improve: LL = " << fr->LL() << " ; min_LL = " << min_LL);
                min_LL=fr->LL();
            }
            
            auto ep = fr->getErrorPropagator();
            auto fitFractions = sig.fitFractions( ep );
            fr->addFractions( fitFractions );
                        
            // Estimate the chi2
            auto evaluator = sig.evaluator();
            if(pNames.size()==6){
                vector<double> chi2Hyper = getChi2MuMu(events, eventsMC, evaluator, evaluator, 7, MinEventsChi2, 0. );
                fr->addChi2( chi2Hyper[0], chi2Hyper[1] );
            }
            else{
                Chi2Estimator chi2( events, eventsMC, evaluator, MinEvents(MinEventsChi2), Dim(3*(pNames.size()-1)-7));
                fr->addChi2( chi2.chi2(), chi2.nBins() );
            }
            
            int fixParamsOptionsFile = NamedParameter<int>("fixParamsOptionsFile",0);
            TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
            fr->print();
            fr->writeToFile(logFile);
            fr->printToLatexTable(tableFile);
            fr->writeToOptionsFile(modelFile, fixParamsOptionsFile);
            fr->writeToRootFile( output, seed, 0, fr->LL(), 0, nSig, {0}, {0} );
            output->cd();
            output->Close();

            unsigned int saveWeights   = NamedParameter<unsigned int>("saveWeights",1);
            if( saveWeights ){
                EventList_type eventsPlotMC;
                if(maxIntEvents == -1) eventsPlotMC = eventsMC;
                else {
                    eventsPlotMC = EventList_type(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true),EntryList(entryListPlotMC));
                }
                sig.setMC( eventsPlotMC );
                //sig.prepare();
                makePlotWeightFile(sig,eventsPlotMC, ep, rndm);
            }
            
    }

}
      
      
//  int status = 0;
//  if(outDir != ".") status &= system( ("mkdir -p " + outDir).c_str() );
//  if( status != 0 ) ERROR("Building OutDir directory failed");  
//  else if(outDir != "."){
//        system(("mv " + logFile + " " + outDir + "/").c_str());
//        system(("mv " + tableFile + " " + outDir + "/").c_str());
//        system(("mv " + modelFile + " " + outDir + "/").c_str());
//        system(("mv " + plotFile + " " + outDir + "/").c_str());
//        const string FitWeightFileName = NamedParameter<string>("FitWeightFileName", "Fit_weights.root");  
//        system(("mv " + FitWeightFileName + " " + outDir + "/").c_str());
//  }

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
}
