#include <Rtypes.h>
#include <TH1.h>
#include <dlfcn.h>
#include <memory>
#include <string>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#include "AmpGen/DynamicFCN.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Particle.h"
#include "AmpGen/RecursivePhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"

using namespace AmpGen;

struct FixedLibPDF {
  void* lib;
  AmpGen::DynamicFCN<double( const double*, const int& )> PDF;
  void prepare(){};
  void setEvents( AmpGen::EventList& evts ){};
  double prob_unnormalised( const AmpGen::Event& evt ) const { return PDF( evt, 1 ); }
  FixedLibPDF( const std::string& lib )
  {
    void* handle = dlopen( lib.c_str(), RTLD_NOW );
    if ( handle == nullptr ) ERROR( dlerror() );
    PDF = AmpGen::DynamicFCN<double( const double*, const int& )>( handle, "FCN" );
  }
  size_t size() { return 0; }
  void reset( const bool& flag = false ){};
};

template <class PDF_TYPE, class PRIOR_TYPE> 
  void GenerateEvents( EventList& events
                       , PDF_TYPE& pdf 
                       , PRIOR_TYPE& prior
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm )
{
  Generator<PRIOR_TYPE> signalGenerator( prior );
  signalGenerator.setRandom( rndm);
  signalGenerator.setBlockSize( blockSize );
  signalGenerator.fillEventList( pdf, events, nEvents );
}

template <class T> void GenerateSource(T& pdf, const std::string& output, MinuitParameterSet& mps)
{
  bool normalise      = NamedParameter<bool>("Normalise",true);
  std::string type    = NamedParameter<std::string>( "Type", "CoherentSum" );
  int seed            = NamedParameter<int>("Seed", 0);
  double safetyFactor = NamedParameter<double>( "SafefyFactor", 3 );
  size_t NormEvents   = NamedParameter<size_t>( "NormEvents", 1000000 );
  auto oEventType     = NamedParameter<std::string>("EventType").getVector();
  EventType eventType( oEventType ); 
  TRandom3 rnd(seed); 
  Generator<PhaseSpace> phsp(eventType, &rnd);
  EventList normEvents = phsp.generate(NormEvents);
  double norm = 1; 
  if( normalise ){
    double pMax = 0 ;
    pdf.setEvents( normEvents );
    pdf.prepare();
    normEvents[0].print();
    normEvents[0].printCache();
    INFO( pdf.prob_unnormalised( normEvents[0] ) );
    pdf.debug( normEvents[0] ); 
    for ( auto& evt : normEvents ) {
      if( type == "PolarisedSum" ){ 
        double px, py, pz; 
        gRandom->Sphere(px,py,pz, gRandom->Uniform(0,1));
        mps["Px"]->setCurrentFitVal(px);
        mps["Py"]->setCurrentFitVal(py);
        mps["Pz"]->setCurrentFitVal(pz);
        pdf.transferParameters();
      }
      double n = pdf.prob_unnormalised( evt );
      if ( n > pMax ) pMax = n;
    }
    norm = pMax * safetyFactor; 
    INFO( "Making source with " << pMax << " x safety factor = " << safetyFactor );
  }
  mps.resetToInit(); 
  pdf.generateSourceCode( output, norm, true );
}


int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 1, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 100000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  bool sourceOnly     = NamedParameter<bool>       ("SourceOnly",false,"Flag to generate source code only"); 
  std::string outfile = NamedParameter<std::string>("Output"   , sourceOnly ? "Generate_Output.root" : "amp.cpp" , "Name of output file" );
  std::string gen_type = NamedParameter<std::string>( "Type", "CoherentSum", optionalHelpString("Generator configuration to use:", 
    { {"CoherentSum", "Full phase-space generator with (pseudo)scalar amplitude"}
    , {"PolarisedSum", "Full phase-space generator with particles carrying spin in the initial/final states"}
    , {"FixedLib", "Full phase-space generator with an amplitude from a precompiled library"}
    , {"RGenerator", "Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons"} } ) );
  
  std::string lib     = NamedParameter<std::string>("Library","","Name of library to use for a fixed library generation");
  size_t nBins        = NamedParameter<size_t>     ("nBins"     ,100, "Number of bins for monitoring plots." );

  #ifdef _OPENMP
    unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
    unsigned int nCores                    = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported, "Number of cores to use (OpenMP only)" );
    omp_set_num_threads( nCores );
    omp_set_dynamic( 0 );
  #endif

  TRandom3 rand;
  rand.SetSeed( seed + 934534 );
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  
  Particle p;
  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );

  EventList accepted( eventType );

  INFO("Generating events with type = " << eventType );

  if ( gen_type == "CoherentSum" ) {
    CoherentSum sig( eventType, MPS );
    PhaseSpace phsp(eventType,&rand);
    if( !sourceOnly ) GenerateEvents( accepted, sig, phsp , nEvents, blockSize, &rand );
    else GenerateSource(sig, outfile, MPS); 
  } 
  else if ( gen_type == "PolarisedSum" ){
    PolarisedSum sig( eventType, MPS );    
    RecursivePhaseSpace phsp( sig.matrixElements()[0].decayTree.quasiStableTree() , eventType, &rand );
    if( !sourceOnly ) GenerateEvents( accepted, sig, phsp, nEvents, blockSize, &rand );
    else GenerateSource(sig,outfile,MPS);
  }
  else if ( gen_type == "RGenerator" ) {
    CoherentSum sig( eventType, MPS, "" );
    Generator<RecursivePhaseSpace> signalGenerator( sig[0].decayTree.quasiStableTree(), eventType );
    signalGenerator.setRandom( &rand );
    signalGenerator.fillEventList( sig, accepted, nEvents );
  }
  else if ( gen_type == "FixedLib" ) {
    Generator<> signalGenerator( eventType );
    signalGenerator.setRandom( &rand );
    signalGenerator.setBlockSize( blockSize );
    signalGenerator.setNormFlag( false );
    FixedLibPDF pdf( lib  );
    signalGenerator.fillEventList( pdf, accepted, nEvents );
  } 
  else ERROR("Did not recognise configuration: " << gen_type );
  if( sourceOnly ) return true; 

  if( accepted.size() == 0 ) return -1;
  TFile* f = TFile::Open( outfile.c_str(), "RECREATE" );
  accepted.tree( "DalitzEventList" )->Write();
  auto plots = accepted.makeDefaultProjections(Bins(nBins), LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
  if( NamedParameter<bool>("plots_2d",true) == true ){
    auto proj = eventType.defaultProjections(nBins);
    INFO("Making 2D projections...");
    for( size_t i = 0 ; i < proj.size(); ++i ){
      for( size_t j = i+1 ; j < proj.size(); ++j ){
      
        accepted.makeProjection( Projection2D(proj[i],proj[j] ) )->Write(); 
      }
    }
  }  
  INFO( "Writing output file " );
  f->Close();
}
