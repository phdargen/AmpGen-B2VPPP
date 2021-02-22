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
  make_enum(pdfTypes, CoherentSum, PolarisedSum, FixedLib)
  make_enum(phspTypes, PhaseSpace, RecursivePhaseSpace, TreePhaseSpace)
}  

void makePlotWeightFile(PolarisedSum& sig, const EventList_type& eventsPlotMC){
    
    string outDir = NamedParameter<std::string>("outDir", ".");    
    string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName", "Fit_weights.root");  
    FitWeightFileName = outDir + "/" + FitWeightFileName;

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

    auto weightFunction = sig.componentEvaluator();
    
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

void makePlotWeightFile(CoherentSum& sig, const EventList_type& eventsPlotMC){
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

    auto weightFunction = sig.componentEvaluator();
    
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
    /*
    if (param->name().find( "_Re" ) != std::string::npos || param->name().find( "_Im" ) != std::string::npos ||
    param->name().find( "_mass" ) != std::string::npos || param->name().find( "_width" ) != std::string::npos
        || param->name().find( "_alpha" ) != std::string::npos
        || param->name().find( "_beta" ) != std::string::npos){         */
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
    //}
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

void sanityChecks(MinuitParameterSet& mps){

   for(int i=0;i<mps.size();++i)if(mps[i]->name().find( "cut_dim" ) != std::string::npos){ mps.unregister( mps.at(i)); i=0; }

   for(int i=0;i<mps.size();++i)if(mps[i]->name().find( "::Spline" ) != std::string::npos)if(mps[i]->name().find( "::Spline::" ) == std::string::npos){ mps.unregister( mps.at(i)); i=0; }
    
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
       if(!mps[i]->isFree())continue;
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
        if(name.find( "_Re" ) != std::string::npos || name.find( "_Im" ) != std::string::npos) {            
            vector<string>name_split = split(name,',');
            //INFO(name_split[0]);
            //INFO(name_split[1]);
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

template <typename likelihoodType>
FitResult* doFit( likelihoodType&& likelihood, EventList_type& data, EventList_type& mc, MinuitParameterSet& MPS )
{
    auto time_wall = std::chrono::high_resolution_clock::now();
    auto time      = std::clock();
    /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
    Minimiser mini( likelihood, &MPS );
    addGaussianConstraint( mini, MPS );

    auto threeBodyShapes     = threeBodyCalculators( MPS );
    unsigned int updateWidth = NamedParameter<unsigned int>( "UpdateWidth", 0 );
    
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


vector<double> getChi2MuMu(const EventList_type& dataEvents, const EventList_type& mcEvents, const std::function<double( const Event& )>& fcn, const int dim = 7, const int minEventsPerBin = 25, double minBinWidth = 0., int mode = 0 ){
    
    EventType pdg = dataEvents.eventType();
 
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

    if(dim==7){
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
    
    for( auto& evt : dataEvents ){
        HyperPoint point( dim );
        if(mode ==0){
            point.at(0)= sqrt(evt.s(m012));
            point.at(1)= sqrt(evt.s(m02));
            point.at(2)= sqrt(evt.s(m12));
            point.at(3)= sqrt(evt.s(m3412));
            point.at(4)= sqrt(evt.s(m341));
        }
        else{
            point.at(0)= evt.s(m012);
            point.at(1)= evt.s(m02);
            point.at(2)= evt.s(m12);
            point.at(3)= evt.s(m3412);
            point.at(4)= evt.s(m341);
        }
        if(dim==7){
            point.at(5)= cosThetaMuAngle(evt);
            point.at(6)= chiMuAngle(evt);
        }
        point.addWeight(evt.weight());
        points.push_back(point);
    }

    HyperPointSet pointsMC( dim);
    for( auto& evt : mcEvents ){
        HyperPoint point( dim );
        if(mode ==0){
            point.at(0)= sqrt(evt.s(m012));
            point.at(1)= sqrt(evt.s(m02));
            point.at(2)= sqrt(evt.s(m12));
            point.at(3)= sqrt(evt.s(m3412));
            point.at(4)= sqrt(evt.s(m341));
        }
        else{
            point.at(0)= evt.s(m012);
            point.at(1)= evt.s(m02);
            point.at(2)= evt.s(m12);
            point.at(3)= evt.s(m3412);
            point.at(4)= evt.s(m341);
        }
        if(dim==7){
            point.at(5)= cosThetaMuAngle(evt);
            point.at(6)= chiMuAngle(evt);
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
  std::string weightData = NamedParameter<std::string>("weightData", "weight");  
  std::string intFile = NamedParameter<std::string>("IntegrationSample","","Name of file containing events to use for MC integration.");
  std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
  std::string phspFile = NamedParameter<std::string>("phspFile","","phspFile");

  std::string outDir = NamedParameter<std::string>("outDir", ".");
  std::string logFile = NamedParameter<std::string>("LogFile", "log.txt", "Name of the output log file");
  std::string tableFile = NamedParameter<std::string>("TableFile", "table.tex", "Name of the output log file");
  std::string modelFile = NamedParameter<std::string>("ModelFile", "model.txt", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("ResultsFile", "result.root", "Name of the output plot file");
  
  auto normAmps = NamedParameter<bool>("normAmps", 1);
  auto excludeNorm = NamedParameter<std::string>("excludeNorm", std::vector<std::string>()).getVector();

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
  if(bNamesMC.size()==0)bNamesMC=bNames;
  if(bExtraNamesMC.size()==0)bExtraNamesMC=bExtraNames;
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
    
  auto randomizeStartVals = NamedParameter<bool>("randomizeStartVals", 0);
  if(randomizeStartVals) randomizeStartingPoint(MPS,rndm);
  //sanityChecks(MPS);
    
  auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
  auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
  auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
  auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
  auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();
  phsp_cut filter(cut_dim,cut_limits,invertCut);

  EventType evtType(pNames);
  vector<size_t> entryList,entryListMC,entryListPlotMC;
  auto maxIntEvents = NamedParameter<int>("maxIntEvents", -1);  
  if(useFilter==1){
      EventList_type dummyEvents(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));
      EventList_type dummyEventsMC(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
    
      if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") dummyEvents.transform( scale_transform );      
      if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") dummyEventsMC.transform( scale_transform );
      
      for(int i=0; i< dummyEvents.size(); i++ )if(filter(dummyEvents[i]))entryList.push_back(i);
      for(int i=0; i< dummyEventsMC.size(); i++ )if(filter(dummyEventsMC[i])){
          if(entryListMC.size()<maxIntEvents)entryListMC.push_back(i);
          entryListPlotMC.push_back(i);
       }
  }
  else for(int i = 0; i < maxIntEvents; i++)entryListMC.push_back(i);
    
  EventList_type events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData),EntryList(entryList));
//    EventList_type events(dataFile, evtType, Branches(bNames), ExtraBranches(bExtraNames), GetGenPdf(false), WeightBranch(weightData),EntryList(entryList));
    
  double nSig = events.sumWeights();
  //EventList_type eventsMC(intFile, evtType, Branches(bNamesMC), ExtraBranches(bExtraNamesMC), WeightBranch(weightMC), GetGenPdf(true),EntryList(entryListMC));
  EventList_type eventsMC(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true),EntryList(entryListMC));
    
  if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
    INFO("Changing data units from MeV -> GeV");
    events.transform( scale_transform );
  }
  if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
    INFO("Changing MC units from MeV -> GeV");
    eventsMC.transform( scale_transform );
  }
  
  auto pdfType = NamedParameter<pdfTypes>( "Type", pdfTypes::PolarisedSum);
  if(pdfType==pdfTypes::CoherentSum){
        PolarisedSum sig(evtType, MPS);
        sig.setMC( eventsMC );
        IncoherentSum bkg( evtType, MPS, "Inco" );
        bkg.setMC( eventsMC );

        sig.setWeight( MPS["f_sig"] );
        bkg.setWeight( MPS["f_bkg"] );

        cout << "Number of amplitudes = " << sig.size() << endl;    
        TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();

        auto ll =  make_likelihood(events, sig, bkg);
        sig.prepare();
        sig.normaliseAmps();

        FitResult* fr = doFit(ll, events, eventsMC, MPS );
        MPS.print();
        auto fitFractions = sig.fitFractions( fr->getErrorPropagator() );   
        fr->addFractions( fitFractions );

        double sumFractions(0);
        for( auto& f : fitFractions ){
            if(f.name()=="Sum_B+" || f.name()=="Sum_B0" )sumFractions = f.val();
        }            
        // Estimate the chi2 
        auto evaluator = sig.evaluator();
        auto MinEventsChi2 = NamedParameter<Int_t>("MinEventsChi2", 15, "MinEventsChi2" );
        Chi2Estimator chi2( events, eventsMC, evaluator, MinEvents(MinEventsChi2), Dim(3*(pNames.size()-1)-7));
        fr->addChi2( chi2.chi2(), chi2.nBins() );

        fr->print();
        fr->writeToFile( logFile );
        fr->printToLatexTable(tableFile);
        fr->writeToOptionsFile( modelFile );
        fr->writeToRootFile( output, seed, 0, sig.size(), nSig );
        output->cd();
        output->Close();

        unsigned int saveWeights   = NamedParameter<unsigned int>("saveWeights",1);  
        if( saveWeights ){
            EventList_type eventsPlotMC;
            if(maxIntEvents == -1) eventsPlotMC = eventsMC;
            else {
                eventsPlotMC = EventList_type(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
            }
            sig.setMC( eventsPlotMC );
            sig.prepare();
            makePlotWeightFile(sig,eventsPlotMC);  
        }
  } 
  else if(pdfType==pdfTypes::PolarisedSum){  
        // Signal pdf
        PolarisedSum sig(evtType, MPS);
        sig.setMC( eventsMC );
        checkAmps(sig, MPS);
        sanityChecks(MPS);
        cout << "Number of amplitudes = " << sig.numAmps() << endl;
            
        auto ll = make_likelihood(events, sig);
        sig.prepare();
      
        if(normAmps){
            if(phspFile == ""){
                sig.normaliseAmps(excludeNorm);
            }
            else{
                auto bNamesPhsp = NamedParameter<std::string>("BranchesPhsp", std::vector<std::string>()).getVector();
                std::string weightPhsp = NamedParameter<std::string>("weightPhsp", "weight");
                EventList_type eventsPhspMC = EventList_type(phspFile, evtType, GetGenPdf(true),Branches(bNamesPhsp), WeightBranch(weightPhsp));
                sig.setMC( eventsPhspMC );
                sig.prepare();
                sig.normaliseAmps(excludeNorm);

                sig.setMC( eventsMC );
                sig.prepare();
            }
        }
            
        // Do fit          
        auto nFits = NamedParameter<int>("nFits", 1);  
        double min_LL = 99999;
    
        for (unsigned int i = 0; i< nFits; ++i) {

            cout << "==============================================" << endl;
            INFO("Start fit number " << i);

            if(i>0){
                randomizeStartingPoint(MPS,rndm);
                sig.setMC( eventsMC );
                sig.prepare();
                //sig.normaliseAmps(excludeNorm);
            }
            FitResult* fr = doFit(ll, events, eventsMC, MPS );
            
            if(fr->LL()>min_LL || TMath::IsNaN(fr->LL())){
                INFO("Fit did not improve: LL = " << fr->LL() << " ; min_LL = " << min_LL);
                continue;
            }
            else {
                INFO("Fit did improve: LL = " << fr->LL() << " ; min_LL = " << min_LL);
                min_LL=fr->LL();
            }
            auto fitFractions = sig.fitFractions( fr->getErrorPropagator() );   
            fr->addFractions( fitFractions );
            vector<double> thresholds{0.1,0.5,1,2,5};
            vector<double> numFracAboveThresholds = sig.numFracAboveThreshold(thresholds);
            
            // Plot spline 
            string head = NamedParameter<std::string>("Head","");
            if(head!="")fr->plotSpline(head,outDir);
          
            // Estimate the chi2 
            auto evaluator = sig.evaluator();
            auto MinEventsChi2 = NamedParameter<Int_t>("MinEventsChi2", 15, "MinEventsChi2" );
            if(pNames.size()==6){
                vector<double> chi2Hyper = getChi2MuMu(events, eventsMC, evaluator, 7, MinEventsChi2, 0. );
                //chi2Hyper = getChi2MuMu(events, eventsMC, evaluator, 5, MinEventsChi2, 0. );
                fr->addChi2( chi2Hyper[0], chi2Hyper[1] );
            }
            else{
                Chi2Estimator chi2( events, eventsMC, evaluator, MinEvents(MinEventsChi2), Dim(3*(pNames.size()-1)-7));
                //chi2Hyper = getChi2(events, eventsMC, evaluator, Dim(3*(pNames.size()-1)-7), MinEventsChi2, 0., 1 );
                //chi2.writeBinningToFile("chi2_binning.txt");
                fr->addChi2( chi2.chi2(), chi2.nBins() );
            }
            
            TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
            fr->print();
            fr->writeToFile(logFile);
            fr->printToLatexTable(tableFile);
            fr->writeToOptionsFile(modelFile);
            fr->writeToRootFile( output, seed, 0, sig.numAmps(), nSig, thresholds, numFracAboveThresholds );
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
                sig.prepare();
                makePlotWeightFile(sig,eventsPlotMC);  
                //Chi2Estimator chi2Plot( events, eventsPlotMC, evaluator, MinEvents(MinEventsChi2), Dim(3*(pNames.size()-1)-7));
                //fr->addChi2( chi2Plot.chi2(), chi2Plot.nBins() );
                //fr->print();
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
