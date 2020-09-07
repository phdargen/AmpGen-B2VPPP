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

#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
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

#if ENABLE_AVX2
#include "AmpGen/EventListSIMD.h"
using EventList_type = AmpGen::EventListSIMD;
#else
#include "AmpGen/EventList.h"
using EventList_type = AmpGen::EventList; 
#endif

#include <TGraph.h>
#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>

using namespace std;
using namespace AmpGen;

void makePlotWeightFile(PolarisedSum& sig, const EventList& eventsPlotMC){
    const std::string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName", "Fit_weights.root");  
    TFile* weight_file = TFile::Open(FitWeightFileName.c_str(),"RECREATE");
    weight_file->cd();
    TTree* weight_tree = new TTree("DalitzEventList","DalitzEventList");
    
    auto plot_amps = NamedParameter<string>("plot_amps", std::vector<string>(),"amplitudes to plot" ).getVector();
    auto plot_weights = NamedParameter<string>("plot_weights", std::vector<string>(),"plot weight names" ).getVector();
    plot_amps.push_back("");
    plot_weights.push_back("weight");
    
    vector<double> weights(plot_amps.size(),0.);
    vector<TBranch*> branches; 
    vector<vector<unsigned>> indices;

    auto weightFunction = sig.componentEvaluator(&eventsPlotMC);
    
    for(int i = 0; i < plot_amps.size(); i++){
        cout << "Plotting amp " << plot_amps[i] << " with weight " <<  plot_weights[i] << endl;

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
                getline( ss, substr, 'x' );
                selectKeys.push_back( substr );
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

    for( const auto& evt : eventsPlotMC ){
        auto weightFun = weightFunction(evt);
        for(int i = 0; i < plot_amps.size(); i++){
            weights[i] = 0;
            for( unsigned j = 0 ; j != indices[i].size(); ++j ) weights[i] += evt.weight() * weightFun[indices[i][j]] / evt.genPdf() ; 
        }
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
  }
  return calculators;
}

void randomizeStartingPoint( MinuitParameterSet& MPS, TRandom3& rand, bool SplineOnly = false )
{
  for (auto& param : MPS) {
    if ( param->isFree() == 0 ) continue;
    if ( SplineOnly && param->name().find( "::Spline" ) == std::string::npos ) continue;
    if (param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos ||
    param->name().find( "_mass" ) != std::string::npos || param->name().find( "_width" ) != std::string::npos ){        
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
        double new_val = gRandom->Gaus( param->mean(), param->err()* sigma);
        param->setInit(new_val);
        param->setCurrentFitVal( new_val );
    }
}

// void removeRandomAmps( MinuitParameterSet& mps, TRandom3& rand, int n = 1 )
// {
//   const int N = count_amplitudes(mps);
//   cout << "Removing " << n << " random amplitudes out of " << N << endl;

//   vector<int> removeMe;
//   while (removeMe.size()<n){
//     int rnd = rand.Rndm()*N+1;
//     if(!(std::find(removeMe.begin(), removeMe.end(),rnd)!=removeMe.end()))removeMe.push_back(rnd);
//   }
//   vector<MinuitParameter*> removeMeParam;

//   for(int i = 0; i<removeMe.size();i++)cout << removeMe[i] << endl;

//   for(int i = 0; i<n;i++){
//     unsigned int counter = 0;
//     for ( auto param = mps.cbegin(); param != mps.cend(); ++param ) {
//       if ( ( *param )->name().find( "_Re" ) == std::string::npos ) continue;
//       counter++;
//       cout << counter << " " << ( *param )->name()  << endl; 
//       if(counter==removeMe[i]){
//         TString name(( *param )->name());
//         name.ReplaceAll("_Re","_Im");
//         removeMeParam.push_back(*param);
//         for(auto param_im = mps.cbegin(); param_im != mps.cend(); ++param_im ){
//           if ( ( *param_im )->name().find( name ) != std::string::npos ) removeMeParam.push_back(*param_im);
//         }
//       }
//     }
//   }
  
//   for(int i = 0 ; i< removeMeParam.size();i++){
//       cout << "Removing parameter " << removeMeParam[i]->name() << endl;
//       mps.unregister(removeMeParam[i]);
//   }

// }

// void removeRandomAmps( MinuitParameterSet& mps, int seed, vector<string> ampList, int n =1 )
// {
//   const int N = count_amplitudes(mps);
//   cout << "Removing " << n << " random amplitudes out of the list: " << endl;
//   for(int i = 0; i<ampList.size();i++)cout << ampList[i] << endl;

//   vector<int> candidateList;
//   for(int i=0;i<mps.size();i++){
//       cout << mps[i]->name() << endl;
//       for(int j = 0; j<ampList.size();j++){
//         if (mps[i]->name().find( ampList[j] ) != std::string::npos && mps[i]->name().find("_Re") != std::string::npos) candidateList.push_back(i);
//       }
//   }

//   cout << "Found candidates " << endl;
//   for(int i = 0; i<candidateList.size();i++)cout << candidateList[i] << endl;

//   std::default_random_engine rand(seed);
//   std::shuffle(candidateList.begin(), candidateList.end(), rand);

//   cout << "Found candidates (rnd order)" << endl;
//   for(int i = 0; i<candidateList.size();i++)cout << candidateList[i] << endl;
  
//   vector<MinuitParameter*> removeMeParam;
//   n = n < candidateList.size() ? n : candidateList.size();
//   for(int i = 0; i<n;i++){
//         TString name(mps[candidateList[i]]->name());
//         name.ReplaceAll("_Re","_Im");
//         removeMeParam.push_back(mps[candidateList[i]]);
//         for(auto param_im = mps.cbegin(); param_im != mps.cend(); ++param_im ){
//           if ( ( *param_im )->name().find( name ) != std::string::npos ) removeMeParam.push_back(*param_im);
//         }
//   }
  
//   for(int i = 0 ; i< removeMeParam.size();i++){
//       cout << "Removing parameter " << removeMeParam[i]->name() << endl;
//       mps.unregister(removeMeParam[i]);
//   }

// }

// void sanityChecks(MinuitParameterSet& mps){

//   vector<int> checkList;

//   for(int i=0;i<mps.size();i++){
//       if(mps[i]->isFixed())continue;
//       if((mps[i]->name().find( "_mass" ) != std::string::npos || mps[i]->name().find( "_width" ) != std::string::npos)) checkList.push_back(i);
//   }
  
//   for(int i= 0; i<checkList.size();i++){
//       bool found = false;
//       TString name(mps[checkList[i]]->name());
//       name.ReplaceAll("_mass","");
//       name.ReplaceAll("_width","");     

//       for(int j=0;j<mps.size();j++){
//           if(mps[j]->name().find(name) != std::string::npos && mps[j]->name().find( "_Re" ) != std::string::npos) found = true;
//       }  

//       if(found == false){
//         cout << "Fitting " << mps[checkList[i]]->name() << " but there is no matching amplitude, set to constant" << endl;
//         mps[checkList[i]]->fix();
//       }
//   }

// }


struct phsp_cut {
    phsp_cut(std::vector<unsigned int> dim, std::vector<double> limits, bool invertCut = false):_dim(dim),_limits(limits),_invertCut(invertCut){}
    bool operator()(const Event& evt){
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

template <typename likelihoodType>
FitResult* doFit( likelihoodType&& likelihood, EventList& data, EventList& mc, MinuitParameterSet& MPS )
{
    auto time_wall = std::chrono::high_resolution_clock::now();
    auto time      = std::clock();
    /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
    Minimiser mini( likelihood, &MPS );
    
    auto threeBodyShapes     = threeBodyCalculators( MPS );
    unsigned int updateWidth = NamedParameter<unsigned int>( "UpdateWidth", 0 );
    
    std::vector<TGraph*> rw_old;
    for( auto& shape : threeBodyShapes ) rw_old.push_back(shape.widthGraph(1.));
    for( auto& graph : rw_old) graph->Write();
    
    if ( updateWidth ) {
        for ( auto& shape : threeBodyShapes ) shape.updateRunningWidth( MPS );
    }
    unsigned int nIterations            = NamedParameter<unsigned int>( "nIterations", 0 );
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
    INFO( "Fitting PDF, iterating " << " " << nIterations + 1 << " times" );
    for ( unsigned int iteration = 0; iteration < nIterations + 1; ++iteration ) {
        mini.doFit();
        if ( iteration == 0 && nIterations != 0 ) {
            for ( auto& shape : threeBodyShapes ) shape.updateRunningWidth( MPS );
            for ( auto& param : slowParamPtrs ) param->setFree(); /// release the parameter ///
        }
    }
    
    std::vector<TGraph*> rw_new;
    for( auto& shape : threeBodyShapes ) rw_new.push_back(shape.widthGraph(1.));
    for( auto& graph : rw_new) graph->Write();
    
    int nTries = 1;
    unsigned int nReTries = NamedParameter<unsigned int>( "nReTries", 5 );
    while(mini.status()>0 && nTries-1 < nReTries){
        INFO("Fit not converged, try again with small perturbation");
        perturb(MPS,nTries);
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
    
    /* Make the plots for the different components in the PDF, i.e. the signal and backgrounds. 
     The structure assumed the PDF is some SumPDF<eventListType, pdfType1, pdfType2,... >. */
    //auto      nBins     = NamedParameter<Int_t>     ("nBins"      , 50           , "Number of bins" );
    //unsigned int counter = 1;
    //for_each(likelihood.pdfs(), [&](auto& pdf){
        //auto pfx = PlotOptions::Prefix("Model");
        //auto mc_plot3 = mc.makeDefaultProjections(WeightFunction(pdf),pfx,PlotOptions::Bins(nBins) );
        //for( auto& plot : mc_plot3 )
        //{
            //plot->Scale( ( data.integral() * pdf.getWeight() ) / plot->Integral() );
            //plot->Write();
        //}
        //counter++;
    //});
    
    /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
    //Chi2Estimator chi2( data, mc, likelihood, MinEvents(25) );
    //chi2.writeBinningToFile("chi2_binning.txt");
    //fr->addChi2( chi2.chi2(), chi2.nBins() );
    
    //fr->print();
    return fr;
}

int main( int argc, char* argv[] )
{
  time_t startTime = time(0);
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "model.txt", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("ResultsFile"     , "results.root", "Name of the output plot file");
  
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();

  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
  auto        nBins    = NamedParameter<Int_t>     ("nBins"      , 50           , "Number of bins" );

  if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

  TRandom3 rndm;
  seed +=  atoi(argv[1]);
  cout << "Using random seed = " << seed << endl;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("LogFile: " << logFile << "; Plots: " << plotFile );
  
#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  auto removeRandom = NamedParameter<int>("removeRandomAmps", 0);
  auto removeRandomAmpsList = NamedParameter<string>("removeRandomAmpsList", std::vector<string>()).getVector();
  //if(removeRandom>0) removeRandomAmps(MPS,seed,removeRandomAmpsList,removeRandom);
  auto randomizeStartVals = NamedParameter<bool>("randomizeStartVals", 0);
  if(randomizeStartVals) randomizeStartingPoint(MPS,rndm);
  //sanityChecks(MPS);

  EventType evtType(pNames);
  EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch("weight"));
  
  auto maxIntEvents = NamedParameter<int>("maxIntEvents", -1);  
  vector<size_t> entryListMC;
  for(int i = 0; i < maxIntEvents; i++)entryListMC.push_back(i);
  EventList eventsMC(intFile, evtType, Branches(bNames), WeightBranch("weight"), GetGenPdf(true),EntryList(entryListMC));
  
  auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
  auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
  auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
  auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();

  phsp_cut filter(cut_dim,cut_limits,invertCut);
  if(useFilter==1)events.filter(filter);
  if(useFilter==1)eventsMC.filter(filter);

  //auto integratorEventFraction = NamedParameter<double>("integratorEventFraction", 1);
  //rnd_cut filter_rnd(integratorEventFraction,false,eventsMC.size());
  //if(integratorEventFraction < 1)eventsMC.filter(filter_rnd);

    /*
  PolarisedSum sig_norm(evtType, MPS);
  sig_norm.setMC( eventsMC );   
  auto ll_norm = make_likelihood(events, sig_norm);
  sig_norm.prepare();
  cout << sig_norm(events[0]) << endl;
*/
    
  PolarisedSum sig(evtType, MPS);
  //for ( unsigned int i = 0; i < sig.numAmps(); ++i )sig.scaleCoupling(i,1./sqrt(sig_norm.norm(i,i).real()));
  //for ( unsigned int i = 0; i < sig.numAmps(); ++i )sig.scaleCoupling(i,2.);
  sig.setMC( eventsMC );
  cout << "Number of amplitudes = " << sig.numAmps() << endl;
//  for ( unsigned int i = 0; i < sig.numAmps(); ++i )sig.scaleCoupling(i,1./sqrt(sig_norm.norm(i,i).real()));
    auto ll = make_likelihood(events, sig);
    sig.prepare();
    cout << ll(events[0]) << endl;
    sig.normalizeAmps();
  //sig.reset();
  //sig.prepare();
  //sig.updateNorms();

/*
    sig.prepare();
    cout << ll(events[0]) << endl;
        
    for ( unsigned int i = 0; i < sig.numAmps(); ++i ) {
        cout <<  sig[i].decayDescriptor() << endl;    
        cout <<  sig.norm(i,i).real() << endl;
        cout <<  sig[i].coefficient << endl << endl;    
        //sig.scaleCoupling(i,2.);
    }
    //sig.reset();
    
    sig.normalizeAmps();
    cout << endl << "scaled" << endl;
    //sig.prepare();
//    cout << ll(events[0]) << endl;
//
    for ( unsigned int i = 0; i < sig.numAmps(); ++i ) {
        cout <<  sig[i].decayDescriptor() << endl;    
        cout <<  sig.norm(i,i).real()  << endl;
        cout <<  sig[i].coefficient << endl << endl;          
    }
    
    //throw "";
*/
    
  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
  /* Do the fit and return the fit results, which can be written to the log and contains the 
     covariance matrix, fit parameters, and other observables such as fit fractions */
  FitResult* fr = doFit(ll, events, eventsMC, MPS );

  /* Calculate the `fit fractions` using the signal model and the error propagator (i.e. 
     fit results + covariance matrix) of the fit result, and write them to a file. 
   */
  auto fitFractions = sig.fitFractions( fr->getErrorPropagator() );   
  fr->addFractions( fitFractions );

  double sumFractions(0);
  for( auto& f : fitFractions ){
      if(f.name()=="Sum_B+")sumFractions = f.val();
  }
    
  /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
  auto evaluator = sig.evaluator(&eventsMC);
  auto MinEventsChi2 = NamedParameter<Int_t>("MinEventsChi2", 15, "MinEventsChi2" );
  Chi2Estimator chi2( events, eventsMC, evaluator, MinEvents(MinEventsChi2) );
  //chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2.chi2(), chi2.nBins() );

  fr->print();
  fr->writeToOptionsFile( logFile );
  fr->writeToRootFile( output, seed, sig.numAmps(),sumFractions );
  output->cd();
  output->Close();

  unsigned int saveWeights   = NamedParameter<unsigned int>("saveWeights",1);  
  if( saveWeights ){
      EventList eventsPlotMC;
      if(maxIntEvents == -1) eventsPlotMC = eventsMC;
      else {
        eventsPlotMC = EventList(intFile, evtType, Branches(bNames), WeightBranch("weight"), GetGenPdf(true));
        if(useFilter==1)eventsPlotMC.filter(filter);
          auto plotEventFraction = NamedParameter<double>("plotEventFraction", 1);
          rnd_cut filter_rnd2(plotEventFraction,true,eventsMC.size());
          if(plotEventFraction < 1)eventsMC.filter(filter_rnd2);
      }
      sig.setMC( eventsPlotMC );
      sig.prepare();

      makePlotWeightFile(sig,eventsPlotMC);  
  }

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
}

