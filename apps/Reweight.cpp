#include "AmpGen/BinDT.h"
#include "TRandom3.h"
#include "TPad.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TDecompChol.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "AmpGen/BoostedDecisionTree.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
//#include "Publish/Publish.h"
#if ENABLE_AVX2
#include "AmpGen/EventListSIMD.h"
using EventList_type = AmpGen::EventListSIMD;
#else
#include "AmpGen/EventList.h"
using EventList_type = AmpGen::EventList; 
#endif
using namespace AmpGen;
using namespace std;

struct CorrelatedGaussian {
  CorrelatedGaussian( const std::vector<double>& means, 
      const TMatrixD& reducedCovariance, 
      TRandom3* rnd ) : 
    m_parameters( means ), 
    m_decomposedCholesky( means.size(),means.size() ),
    m_rand(rnd) { 

      TDecompChol decomposed( reducedCovariance );
      decomposed.Decompose();
      m_decomposedCholesky = decomposed.GetU();
      transpose();
    }
  std::vector<double> get() {
    const unsigned int N = m_decomposedCholesky.GetNrows();
    TVectorD e(N);
    for( unsigned int i = 0 ; i < N ; ++i ) e[i] = m_rand->Gaus(0,1);

    TVectorD p = m_decomposedCholesky*e; // perturbation in the usual parameter basis
    std::vector<double> rt( N );
    for(int j = 0 ; j < p.GetNrows(); ++j){
      rt[j] = m_parameters[j] + p[j] ;
    }
    return rt;
  }

  void transpose(){
    for( int i = 0 ; i < m_decomposedCholesky.GetNrows(); ++i){
      for( int j = i+1;j<m_decomposedCholesky.GetNrows();++j){
        double tmp = m_decomposedCholesky(j,i);
        m_decomposedCholesky(j,i) = m_decomposedCholesky(i,j);
        m_decomposedCholesky(i,j) = tmp;
      }
    }
  }
  std::vector<double> m_parameters; 
  TMatrixD            m_decomposedCholesky;
  TRandom3*           m_rand;
};

struct TestEvent {
  double x; 
  double y; 
  double w;
  TestEvent( double x, double y, double w ) : x(x),y(y),w(w) {};
};
struct BinningWithWeights {
  BinDT m_binning; 
  std::vector<double> m_leafWeights;
  BinningWithWeights() {};
  BinningWithWeights( BinDT binning ) : m_binning(binning) {} 
  void generateLeafWeights( int dim, std::vector<double*> sourcePointers, std::vector<double*> targetPointers ){
    m_binning.makeNodes( sourcePointers, targetPointers );
    std::vector<double> weight_source(m_binning.const_nodes().size(),0);
    std::vector<double> weight_target(m_binning.const_nodes().size(),0);
    for( auto& evt : targetPointers ) weight_target[m_binning.getBinNumber(evt)] += evt[dim];
    for( auto& evt : sourcePointers ) weight_source[m_binning.getBinNumber(evt)] += evt[dim];
    m_leafWeights.resize( m_binning.const_nodes().size() );
    for( unsigned int i = 0 ; i < m_leafWeights.size(); ++i ){
      m_leafWeights[i] = weight_target[i] / weight_source[i] ;
      INFO("Leaf["<<i<<"] weight = " << m_leafWeights[i] << " " <<
         weight_target[i] << " / " << weight_source[i] 
          ); 
    }
  }
  double getEventWeight( double* event ){
    return m_leafWeights[ m_binning.getBinNumber( event ) ] ; 
  }
};

int test(int argc, char** argv){
  //setupStyle();

  std::vector<TestEvent> source; 
  std::vector<TestEvent> target;

  TRandom3* rnd = new TRandom3(7);
  TMatrixD covMatrix(2,2);
  covMatrix(0,0) = covMatrix(1,1) = 1;
  covMatrix(1,0) = covMatrix(0,1) = 0.5;
  CorrelatedGaussian FG({0.2,0.0},covMatrix, rnd );

  unsigned int nEvents=500000;
  unsigned int nBins2D=40;
  TH1D* histTarget = new TH1D("histTarget","",100,-4.0,4.0);
  TH1D* histSource = new TH1D("histSource","",100,-4.0,4.0);
  TH1D* histReweighted = new TH1D("histReweighted","",100,-4.0,4.0);

  TH2D* histTarget2D = new TH2D("histTarget2D","",nBins2D,-4.0,4.0,nBins2D,-4.0,4.0);
  TH2D* histSource2D = new TH2D("histSource2D","",nBins2D,-4.0,4.0,nBins2D,-4.0,4.0);
  for( unsigned int i = 0 ; i < nEvents ; ++i ){
    TestEvent targetEvent( rnd->Gaus(), rnd->Gaus(), 1);
    auto cm = FG.get();
    TestEvent sourceEvent( cm[0], cm[1], 1);
    target.push_back( targetEvent );
    source.push_back( sourceEvent );
    histTarget->Fill( targetEvent.x );
    histSource->Fill( sourceEvent.x );
    histTarget2D->Fill( targetEvent.x, targetEvent.y );
    histSource2D->Fill( sourceEvent.x, sourceEvent.y );
  }
  histTarget2D->Draw("COLZ");
  gPad->SaveAs("HistTarget2D.pdf");
  histSource2D->Draw("COLZ");
  gPad->SaveAs("HistSource2D.pdf");

  std::vector<double*> sourcePointers; 
  std::vector<double*> targetPointers;

  for( unsigned int i = 0 ; i < nEvents ; ++i ){
    sourcePointers.push_back( &( source[i].x ) );
    targetPointers.push_back( &( target[i].x ) );
  }
  int dim = 2;
  std::vector<BinningWithWeights> binnings;

  for( unsigned int i = 0 ; i < 10; ++i ){ 
    binnings.emplace_back( BinDT( MaxDepth(4),Dim(dim),MinEvents(50) ) );
    binnings.rbegin()->generateLeafWeights(dim, sourcePointers, targetPointers );
    for( auto& evt : source ) evt.w *= binnings.rbegin()->getEventWeight( &(evt.x) );
  }

  std::vector<TestEvent> newEvents;
  for( unsigned int i = 0 ; i < nEvents ; ++i ){
    auto cm = FG.get();
    TestEvent te( cm[0], cm[1],1);
    for( auto& binning : binnings ) te.w *= binning.getEventWeight( &(te.x) );
    newEvents.push_back( te );
  } 
  
  TH2D* histSourceReweighted2D = new TH2D("histSourceReweighted2D","",40,-4.0,4.0,40,-4.0,4.0);
  for( auto& evt : newEvents ){
    histReweighted->Fill( evt.x, evt.w );
    histSourceReweighted2D->Fill( evt.x, evt.y, evt.w );
  }
  TH2D* pull = new TH2D("pull","",40,-4.0,4.0,40,-4.0,4.0);
  TH1D* integrated_pull = new TH1D("integrated_pull","",50,-6,6);
  for( int x = 1 ; x <= 40; ++x){
    for( int y = 1 ; y <= 40; ++y){
      double delta = histSourceReweighted2D->GetBinContent(x,y) - histTarget2D->GetBinContent(x,y);
      if( histSourceReweighted2D->GetBinContent(x,y) < 10 ) continue ; 
      double error = sqrt( 
            pow( histSourceReweighted2D->GetBinError(x,y) , 2 )
          + pow( histTarget2D->GetBinError(x,y)           , 2 ) );
      if( error != 0 ){ 
        pull->SetBinContent(x,y,fabs(delta)/error);
        integrated_pull->Fill( delta/error );
      INFO( x << ", " << y << ", " << delta/error << " "
         << histSourceReweighted2D->GetBinContent(x,y) << " " <<  histTarget2D->GetBinContent(x,y) );
      }
    }
  }

  histTarget->Draw();
  histSource->SetLineColor(kRed);
  histReweighted->SetLineColor(kBlue);
  histReweighted->SetMarkerSize(0);
  histSource->Draw("same");
  histReweighted->Draw("E same");

  //TLegend leg; 
  //leg.addEntry( histTarget ,"Data");
  //leg.addEntry( histSource , "MC");
  //leg.addEntry( histReweighted, "Reweighted MC");

  //leg.Draw(0.175,0.9);
  gPad->SaveAs("hist.pdf");
  histSourceReweighted2D->Draw("COLZ");
  gPad->SaveAs("hist2D_reweighted.pdf");
  pull->Draw("COLZ");
  gPad->SaveAs("pull.pdf");

  integrated_pull->Draw("HIST");

  gPad->SaveAs("ipull.pdf");
  TFile* output = TFile::Open("ReweightPlots.root","RECREATE");
  integrated_pull->Write();
  histSourceReweighted2D->Write();
  pull->Write();
  output->Close();
} 

std::tuple<EventList, EventList> split( EventList& events ){
    std::tuple<EventList, EventList> splitEventList = std::make_tuple( events.eventType(), events.eventType() );;
    for( auto& event : events ){
        ( gRandom->Uniform() > 0.5 ? std::get<0>(splitEventList ) : std::get<1>( splitEventList)).push_back( event );
    }
    INFO("Split eventlist of size: " << events.size() << " into " << std::get<0>(splitEventList).size() << " and " << std::get<1>(splitEventList).size() << " event lists");
    events.clear();
    return splitEventList;
}

int main( int argc, char** argv ){
    OptionsParser::setArgs( argc, argv );

    EventType evtType(NamedParameter<std::string>("EventType").getVector());    

    std::string output  = NamedParameter<std::string>("Output","bdt.txt");
    size_t nTrees       = NamedParameter<size_t>("nTrees",300);
    size_t maxDepth     = NamedParameter<size_t>("MaxDepth",7);
    double learningRate = NamedParameter<double>("LearningRate",0.06);
    size_t minEvents    = NamedParameter<size_t>("MinEvents",25);
    
    std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
    std::string weightData = NamedParameter<std::string>("weightData", "weight");  
    std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
    std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
    
    auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>() ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    auto bNamesMC = NamedParameter<std::string>("BranchesMC", std::vector<std::string>() ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    if(bNamesMC.size()==0)bNamesMC=bNames;
    
    EventList_type events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));
    EventList_type eventsMC(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(false));
    
    auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
    if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing data units from MeV -> GeV");
        events.transform( scale_transform );
    }
    if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing MC units from MeV -> GeV");
        eventsMC.transform( scale_transform );
    }
    
    const std::string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName","Fit_weights.root");  
    TFile* weight_file = TFile::Open(FitWeightFileName.c_str(),"OPEN");
    weight_file->cd();
    auto weight_tree = (TTree*) weight_file->Get("DalitzEventList");
    if(weight_tree->GetEntries() != eventsMC.size()){
        cout << "ERROR inconsistent number of events" << endl;
        throw "ERROR";
    }
    
    double w;
    weight_tree->SetBranchAddress("weight",&w);
    for(int i=0; i< eventsMC.size(); i++ )
    {
        weight_tree->GetEntry(i);
        eventsMC.setWeight(i,w);
    }
    
    auto projections = evtType.defaultProjections(50); 
    
    auto sourceEvts = eventsMC;
    auto targetEvts = events;
    
    auto split_sourceEvts = split( sourceEvts );
    //auto sourceTrain      = sourceEvts; 
    //auto sourceTest       = sourceEvts; 
    auto sourceTrain      = std::get<0>(split_sourceEvts); 
    auto sourceTest       = std::get<1>(split_sourceEvts); 

    //auto split_targetEvts = split( targetEvts );
    auto targetTrain     = targetEvts; // std::get<0>(split_targetEvts); 
    auto targetTest      = targetEvts; //std::get<1>(split_targetEvts); 
    
    std::function< std::vector<double>( const Event& )> functors = 
    [=](const Event& event ) -> std::vector<double> { return {  
        sqrt(event.s(1,2,3)),
        sqrt(event.s(1,3)),
        sqrt(event.s(2,3)),
        sqrt(event.s(0,2,3)),
        sqrt(event.s(0,2)),
        event.weight() }; };
    
    size_t dim = 5;
    
    auto targetTestPlots  = targetTest.makeProjections(projections , PlotOptions::Prefix("TargetTest"));
    auto targetTrainPlots = targetTrain.makeProjections(projections, PlotOptions::Prefix("TargetTrain"));
    auto sourceTestPlots  = sourceTest.makeProjections(projections , PlotOptions::Prefix("SourceTest"));
    auto sourceTrainPlots = sourceTrain.makeProjections(projections, PlotOptions::Prefix("SourceTrain"));
    
    BoostedDecisionTree bdt(Store(sourceTrain, functors, dim), 
                            Store(targetTrain, functors, dim),
                            nTrees, maxDepth, learningRate, minEvents );
    
    bdt.writeToFile(output); 
    
    //auto reweight =  [&bdt,&functors](auto& event){ event.setWeight(event.weight()*bdt.weight( functors(event).data())); };
    auto reweight =  [&bdt,&functors](auto& event){ event.setWeight(bdt.weight( functors(event).data())); };
    auto reweightRes =  [&bdt,&functors](auto& event){ event.setWeight(event.weight()/(bdt.weight( functors(event).data()) ) ) ; };

    EventList_type events_copy(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));    
    if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") events_copy.transform( scale_transform );
    events_copy.transform(reweightRes);
    TFile* outputData = TFile::Open("OutputData.root","RECREATE");
    outputData->cd(0);
    events_copy.tree("DalitzEventList")->Write();
    outputData->Close();

    EventList_type eventsMC_copy(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
    if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") eventsMC_copy.transform( scale_transform );    
    eventsMC_copy.transform(reweightRes);
    TFile* outputMC = TFile::Open("OutputMC.root","RECREATE");
    outputMC->cd(0);
    eventsMC_copy.tree("DalitzEventList")->Write();
    outputMC->Close();

    targetTrain.transform( reweightRes );
    sourceTest.transform( reweight );
    sourceTrain.transform( reweight );
    
    //INFO( "Original" << sourceEvts_copy.size() );
    INFO( "Source test has: " << sourceTest.size() << " events; train has: " << sourceTrain.size() << " events");
    
    TFile* outputHists = TFile::Open("Output.root","RECREATE");
    auto sourceTestPlots_reweighted  =  sourceTest.makeProjections(projections, PlotOptions::Prefix("SourceReweighted")     );
    auto sourceTrainPlots_reweighted = sourceTrain.makeProjections(projections, PlotOptions::Prefix("SourceTrainReweighted"));
    auto targetTrainPlots_reweighted = targetTrain.makeProjections(projections, PlotOptions ::Prefix("TargetTrainReweighted"));

    
    TCanvas* c = new TCanvas();
    for( unsigned int i = 0 ; i < sourceTestPlots.size(); ++i){
        THStack* hs = new THStack();  
        sourceTestPlots[i]->Scale( 1. / sourceTestPlots[i]->Integral() );
        sourceTestPlots_reweighted[i]->Scale( 1. / sourceTestPlots_reweighted[i]->Integral() );
        targetTrainPlots_reweighted[i]->Scale( 1. / targetTrainPlots_reweighted[i]->Integral() );
        targetTestPlots[i]->Scale( 1. / targetTestPlots[i]->Integral() );

        hs->Add( sourceTestPlots[i], "HIST");
        hs->Add( sourceTestPlots_reweighted[i], "HIST");
        hs->Add( targetTestPlots[i], "E"); 
        hs->Add( targetTrainPlots_reweighted[i], "HIST"); 

        sourceTestPlots[i]->SetLineColor(kBlue); 
        sourceTestPlots_reweighted[i]->SetLineColor(kRed);
        targetTrainPlots_reweighted[i]->SetLineColor(kGreen+3); 
 
        TLegend leg(0.6,0.8,0.9,0.9,"");
        leg.AddEntry(targetTrainPlots[i]          , "Target", "E");
        leg.AddEntry(sourceTestPlots[i]           , "Source" );
        leg.AddEntry(sourceTestPlots_reweighted[i], "Source [weighted]"   );
        leg.SetTextSize(0.05);
        hs->Draw("nostack");
        leg.Draw();
        
        hs->GetXaxis()->SetTitle( sourceTestPlots[i]->GetXaxis()->GetTitle() );
        hs->GetYaxis()->SetTitle( sourceTestPlots[i]->GetYaxis()->GetTitle() );
        hs->GetYaxis()->SetTitleOffset(1.15);
        c->Print((sourceTestPlots[i]->GetName() + std::string(".pdf")).c_str() );
        delete hs;
    } 
    outputHists->Close();
}

