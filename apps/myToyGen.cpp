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
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/enum.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/AddCPConjugate.h"

#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_t = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_t = AmpGen::EventList; 
#endif

using namespace std;
using namespace AmpGen;

namespace AmpGen { 
  make_enum(pdfTypes, CoherentSum, IncoherentSum, PolarisedSum, FixedLib, Phsp)
  make_enum(phspTypes, PhaseSpace, RecursivePhaseSpace, TreePhaseSpace)
}  

struct FixedLibPDF 
{
  void* lib = {nullptr};
  DynamicFCN<double( const double*, int )> PDF;
  void debug( const Event& event) {};
  void prepare(){};
  void setEvents( AmpGen::EventList& evts ){};
  void setEvents( AmpGen::EventListSIMD& evts ){};
  double operator()( const AmpGen::Event& evt ) const { return PDF( evt, 1 ); }
  double operator()( const double* evt, const unsigned& index )
  {
    return PDF(evt, 1 ); 
  } 
  FixedLibPDF( const std::string& lib )
  {
    void* handle = dlopen( lib.c_str(), RTLD_NOW );
    if ( handle == nullptr ) ERROR( dlerror() );
    PDF = DynamicFCN<double( const double*, int )>( handle, "FCN" );
  }
  size_t size() { return 0; }
  void reset( const bool& flag = false ){};
};

template <class T> void generateSource(T& pdf, const std::string& sourceFile, MinuitParameterSet& mps)
{
  bool normalise      = NamedParameter<bool>("Normalise",true);
  double safetyFactor = NamedParameter<double>( "SafetyFactor", 3 );
  int seed            = NamedParameter<int>("Seed", 1);
  size_t nEvents      = NamedParameter<size_t>( "NormEvents", 1000000 );
  
  TRandom3 rnd(seed);

  Generator<PhaseSpace> phsp(pdf.eventType());
  phsp.setRandom(&rnd);
  EventList normEvents = phsp.generate(nEvents);
  //if constexpr( std::is_same<T, CoherentSum>::value )
  pdf.prepare();

  double norm = 1; 
  if( normalise ){
    double pMax = 0;
    for ( auto& evt : normEvents ) 
    {
      if constexpr ( std::is_same<T, PolarisedSum>::value )
      {
        double px, py, pz; 
        rnd.Sphere(px,py,pz, rnd.Uniform(0,1));
        mps["Px"]->setCurrentFitVal(px);
        mps["Py"]->setCurrentFitVal(py);
        mps["Pz"]->setCurrentFitVal(pz);
        pdf.transferParameters();
      }
      double n = 0;
      if constexpr ( std::is_same<T, CoherentSum>::value ) n = std::norm( pdf.getValNoCache(evt) );
      if constexpr ( std::is_same<T, IncoherentSum>::value ) n =  pdf.prob_unnormalisedNoCache(evt);
      if constexpr ( std::is_same<T, PolarisedSum>::value ) n = pdf.getValNoCache(evt);
      if ( n > pMax ) pMax = n;
    }
    norm = pMax * safetyFactor ; 
    INFO( "Making binary with " << pMax << " x safety factor = " << safetyFactor );
  }
  mps.resetToInit(); 
  pdf.generateSourceCode( sourceFile, norm, true );
}


template <typename pdf_t> Particle getTopology(const pdf_t& pdf)
{
  if constexpr( std::is_same<pdf_t, FixedLibPDF>::value )
  {
    FATAL("Cannot deduce decay topology from a compiled library, check generator options");
  }
  else return pdf.matrixElements()[0].decayTree.quasiStableTree();
}

template <typename pdf_t> std::vector<Particle> getDecayChains( const pdf_t& pdf )
{
  if constexpr( std::is_same<pdf_t, FixedLibPDF>::value )
  {
    FATAL("Cannot deduce decay topology from a compiled library, check generator options");
  }
  else {
    std::vector<Particle> channels; 
    for( auto& chain : pdf.matrixElements() ) channels.push_back( chain.decayTree );
    return channels;
  }
}


template <typename pdf_t> void generateEvents( EventList& events
                       , pdf_t& pdf 
                       , const phspTypes& phsp_type
                       , const size_t& nEvents
                       , const size_t& blockSize
                       , TRandom* rndm
                       , const bool& normalise = true )
{
  if constexpr( std::is_same<pdf_t, FixedLibPDF>::value )
  {
    Generator<PhaseSpace, EventList> signalGenerator(events.eventType(), rndm);
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents );
  }
  else { 
  if( phsp_type == phspTypes::PhaseSpace )
  {
    Generator<PhaseSpace, EventList_t> signalGenerator(events.eventType(), rndm);
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents );
  }
  else if( phsp_type == phspTypes::RecursivePhaseSpace )
  {
    Generator<RecursivePhaseSpace, EventList_t> signalGenerator( getTopology(pdf), events.eventType(), rndm );
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents);
  }
  else if( phsp_type == phspTypes::TreePhaseSpace )
  {
    Generator<TreePhaseSpace, EventList_t> signalGenerator(getDecayChains(pdf), events.eventType(), rndm);
    signalGenerator.setBlockSize(blockSize);
    signalGenerator.setNormFlag(normalise);
    signalGenerator.fillEventList(pdf, events, nEvents );
  }
  else {
    FATAL("Phase space configuration: " << phsp_type << " is not supported");
  }
  }
}


void sanityChecks(MinuitParameterSet& mps){

   for(int i=0;i<mps.size();++i)if(mps[i]->name().find( "cut_dim" ) != std::string::npos){ mps.unregister( mps.at(i)); i=0; }

   // check
   for(int i=0;i<mps.size();++i)if(mps[i]->name().find( "::Spline" ) != std::string::npos){ mps.unregister( mps.at(i)); i=0; }
    
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
            if(sig.findAmp(name_split[0])==0){
                INFO("Removing " << name);
                mps.unregister(mps[i]);
                i=0;
            }
        }
    }
}


int main( int argc, char** argv )
{
  time_t startTime = time(0);
  OptionsParser::setArgs( argc, argv );

  size_t nEvents      = NamedParameter<size_t>     ("nEvents"  , 1000, "Total number of events to generate" );
  size_t blockSize    = NamedParameter<size_t>     ("BlockSize", 5000000, "Number of events to generate per block" );
  int seed            = NamedParameter<int>        ("Seed"     , 0, "Random seed used in event Generation" );
  std::string outfile = NamedParameter<std::string>("Output"   , "Generate_Output.root" , "Name of output file" ); 
  auto phspType        = NamedParameter<phspTypes>( "PhaseSpace", phspTypes::PhaseSpace, optionalHelpString("Phase-space generator to use:",
      std::make_pair(phspTypes::PhaseSpace         , "Phase space generation based on Raubold-Lynch algorithm (recommended).\0")
    , std::make_pair(phspTypes::TreePhaseSpace     , "Divides the phase-space into a series of quasi two-body phase-spaces for efficiently generating narrow states.\0")
    , std::make_pair(phspTypes::RecursivePhaseSpace, "Includes possible quasi-stable particles and the phase spaces of their decay products, such as Λ baryons.\0") ) ); 
  size_t nBins        = NamedParameter<size_t>     ("nBins"     ,100, "Number of bins for monitoring plots." );
 
  auto normAmps = NamedParameter<bool>("normAmps", 1);
  auto excludeNorm = NamedParameter<std::string>("excludeNorm", std::vector<std::string>()).getVector();
  auto combineNorm = NamedParameter<std::string>("combineNorm", std::vector<std::string>()).getVector();
  std::string phspFile = NamedParameter<std::string>("phspFile","","phspFile");

  #ifdef _OPENMP
    unsigned int concurentThreadsSupported = std::thread::hardware_concurrency();
    unsigned int nCores                    = NamedParameter<unsigned int>( "nCores", concurentThreadsSupported, "Number of cores to use (OpenMP only)" );
    omp_set_num_threads( nCores );
    omp_set_dynamic( 0 );
  #endif

  TRandom3 rand;
  seed =  atoi(argv[1]);
  cout << "Using random seed = " << seed << endl;
  rand.SetSeed( seed + 934534 );
  gRandom = &rand;

  MinuitParameterSet MPS;
  MPS.loadFromStream();
  
  EventType eventType; 
  std::string decay   = NamedParameter<std::string>("Decay","","Single decay written on the command line"); 
  if( decay != "" )
  {
    Particle p(decay);
    eventType = p.eventType();
    MPS.add(p.decayDescriptor()+"_Re", Flag::Fix, 1., 0);
    MPS.add(p.decayDescriptor()+"_Im", Flag::Fix, 0., 0);
  } 
  else eventType = EventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(), false);
 
  INFO("Writing output: " << outfile );
  #ifdef _OPENMP
  INFO("Using: " << nCores  << " / " << concurentThreadsSupported  << " threads" );
  #endif

  INFO("Generating events with type = " << eventType );
  EventList events( eventType );
  EventList events_bkg( eventType );

  double f_sig = MPS["f_sig"]==nullptr ? 1.0 : MPS["f_sig"]->mean();
  double f_bkg = MPS["f_bkg"]==nullptr ? 0.0 : MPS["f_bkg"]->mean();
  INFO("Generating " << nEvents << " events with f_sig = " << f_sig << " and f_bkg = " << f_bkg << " ..." );
  TFile* f;
    
  // Generate signal events
  if(f_sig>0){
        {
            PolarisedSum pdf(eventType, MPS);
            checkAmps(pdf, MPS);
            sanityChecks(MPS);
            cout << "Number of amplitudes = " << pdf.numAmps() << endl;
            if(normAmps){
                auto bNamesPhsp = NamedParameter<std::string>("BranchesPhsp", std::vector<std::string>()).getVector();
                std::string weightPhsp = NamedParameter<std::string>("weightPhsp", "weight");
                EventList eventsPhspMC = EventList(phspFile, eventType, GetGenPdf(true),Branches(bNamesPhsp), WeightBranch(weightPhsp));
                EventList dummy( eventType );
                dummy.push_back(eventsPhspMC[0]);
                
                pdf.setEvents( dummy );
                pdf.setMC( eventsPhspMC );
                pdf.prepare();
                INFO("Normalizing amps using the file " << phspFile);
                pdf.normaliseAmps(excludeNorm);
                pdf.reset(true);
            }
            INFO("Generating " << f_sig * nEvents << " signal events ...");
            generateEvents( events, pdf, phspType, f_sig * nEvents, blockSize, &rand );
        }
        TString outfile_sig = TString(outfile.c_str()).ReplaceAll(".root","_sig.root");
        f = new TFile( outfile_sig, "RECREATE" );
        INFO( "Writing signal output file " << outfile_sig );
        events.tree( "DalitzEventList" )->Write();
        f->Close();
  }
  
  // Generate bkg events
  if(f_bkg>0){
        {
            IncoherentSum bkg( eventType, MPS, "Inco");
            INFO("Generating " << f_bkg * nEvents << " bkg events ...");
            generateEvents(events_bkg, bkg, phspType , f_bkg * nEvents, blockSize, &rand );
        }
        TString outfile_bkg = TString(outfile.c_str()).ReplaceAll(".root","_bkg.root");
        f = TFile::Open( outfile_bkg, "RECREATE" );
        INFO( "Writing bkg output file " << outfile_bkg );
        events_bkg.tree( "DalitzEventList" )->Write();
        f->Close();
  }
    
  // Merge sig + bkg events
  events.add(events_bkg);
  f = TFile::Open( outfile.c_str(), "RECREATE" );
  events.tree( "DalitzEventList" )->Write();
  auto plots = events.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack));
  for ( auto& plot : plots ) plot->Write();
  INFO( "Writing output file " << outfile );
  f->Close();
    
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
}
