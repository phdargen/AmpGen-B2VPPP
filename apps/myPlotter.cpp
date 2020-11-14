#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/LHCbStyle.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/Kinematics.h"

#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLorentzVector.h>

using namespace AmpGen;
using namespace std;

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

vector<TH1D*> createHistos(vector<unsigned int> dim,string name, string title, int nBins, vector<double> limits, vector<string> weights){

  vector<TH1D*> histos;
  TH1D* histo = new TH1D(name.c_str(),"",dim.size()==1 ? nBins/2 : nBins,limits[0],limits[1]);
  histo->SetMinimum(0.);
  histo->GetXaxis()->SetTitle(title.c_str());
  histo->GetYaxis()->SetTitle("Yield (norm.)");
  histo->SetMarkerSize(1.);
	//histo->SetMarkerStyle(21);
  histos.push_back(histo);

  for(int i=1; i<weights.size();i++){
      TH1D* h = (TH1D*) histo->Clone((name+weights[i]).c_str());

      if(i==1){
            h->SetLineColor(kBlue);
            h->SetLineWidth(3);
      }else if(i==2){
            h->SetLineColor(kRed+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kRed+1);
            //h->SetFillStyle(3353);
      }else if(i==7){
            h->SetLineColor(kGreen+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kGreen+3);
            //h->SetFillStyle(3353);
      }else if(i==4){
            h->SetLineColor(kMagenta+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kMagenta+3);
            //h->SetFillStyle(3353);
      }else if(i==5){
            h->SetLineColor(kBlack);
            h->SetLineWidth(3);
            //h->SetLineStyle(kDashed);
      }else if(i==6){
            h->SetLineColor(kCyan+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kGray+3);
            //h->SetFillStyle(1001);
      }else if(i==3){
            h->SetLineColor(kGray+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kGray+3);
            //h->SetFillStyle(1001);
      }else if(i==8){
          h->SetLineColor(kOrange+1);
          h->SetLineWidth(3);
          //h->SetFillColor(kGray+3);
          //h->SetFillStyle(1001);
      }else {
          h->SetLineColor(kBlack);
          h->SetLineWidth(3);
          h->SetLineStyle(kDashed);
      }
      histos.push_back(h);
  }
  return histos;
}

void plotHistos(vector<TH1D*>histos,bool plotComponents = true, int style = 0){

  if(style == 1){
    histos[0]->SetLineWidth(1);
    histos[0]->SetMarkerSize(.5);
    histos[0]->GetYaxis()->SetTitle("");
  }
  histos[0]->SetMinimum(1);
  if(plotComponents == true && style == 0)histos[0]->SetMaximum(histos[0]->GetMaximum()*1.3);
  histos[0]->DrawNormalized("",1);

  //for (int i = (plotComponents == true ? histos.size()-1 : 1); i > 0; i--)
  for (int i = 1; i < (plotComponents == true ? histos.size() : 2); i++)
  {
      histos[i]->SetMinimum(1);
      if(style == 1)histos[i]->SetLineWidth(2);
      double norm = histos[i]->Integral()/histos[1]->Integral();
      histos[i]->DrawNormalized("histcsame",norm);
  }
  histos[0]->DrawNormalized("same",1);
}

vector<TH2D*> createHistos2D(vector<unsigned int> dim1, vector<unsigned int> dim2, string name, string title, int nBins, vector<double> limits1, vector<double> limits2, vector<string> weights){

  vector<TH2D*> histos;
  TH2D* histo = new TH2D(name.c_str(),title.c_str(),nBins*2,limits1[0],limits1[1],nBins*2,limits2[0],limits2[1]);
  histo->SetMinimum(0.);
  histo->SetMarkerSize(0.1);
	//histo->SetMarkerStyle(21);
  histos.push_back(histo);

  for(int i=1; i<weights.size();i++){
      TH2D* h = (TH2D*) histo->Clone((name+weights[i]).c_str());

      if(i==1){
            h->SetLineColor(kBlue);
            h->SetLineWidth(3);
      }else if(i==2){
            h->SetLineColor(kRed+1);
            h->SetLineWidth(2);
            //h->SetFillColor(kRed+1);
            //h->SetFillStyle(3353);
      }else if(i==7){
            h->SetLineColor(kGreen+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kGreen+3);
            //h->SetFillStyle(3353);
      }else if(i==4){
            h->SetLineColor(kMagenta+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kMagenta+3);
            //h->SetFillStyle(3353);
      }else if(i==5){
            h->SetLineColor(kBlack);
            h->SetLineWidth(3);
            //h->SetLineStyle(kDashed);
      }else if(i==6){
            h->SetLineColor(kCyan+1);
            h->SetLineWidth(3);
            //h->SetFillColor(kGray+3);
            //h->SetFillStyle(1001);
      }else if(i==3){
            h->SetLineColor(kGray+1);
            h->SetLineWidth(2);
            //h->SetFillColor(kGray+3);
            //h->SetFillStyle(1001);
      }else h->SetLineColor(i+1);

      histos.push_back(h);
  }
  return histos;
}

void makePlots(){

  std::string outDir = NamedParameter<std::string>("outDir", ".");
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string weightData = NamedParameter<std::string>("weightData", "weight");  
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto bNamesMC = NamedParameter<std::string>("BranchesMC", std::vector<std::string>() ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  if(bNamesMC.size()==0)bNamesMC=bNames;
  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  EventType evtType(pNames);
  EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData) );
  EventList eventsMC(intFile, evtType, Branches(bNamesMC), GetGenPdf(true), WeightBranch(weightMC));
  auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
    if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing data units from MeV -> GeV");
        events.transform( scale_transform );
  }
  if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing MC units from MeV -> GeV");
        eventsMC.transform( scale_transform );
  }

  auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
  auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
  auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
  auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();

  phsp_cut filter(cut_dim,cut_limits,invertCut);
  if(useFilter==1)events.filter(filter);
  if(useFilter==1)eventsMC.filter(filter);

    auto invertCut_plot1 = NamedParameter<bool>("invertCut_plot1", 0,"Invert cut logic");
    auto cut_dim_plot1 = NamedParameter<unsigned int>("cut_dim_plot1", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot1 = NamedParameter<double>("cut_limits_plot1", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot1(cut_dim_plot1,cut_limits_plot1,invertCut_plot1);    

    auto invertCut_plot2 = NamedParameter<bool>("invertCut_plot2", 0,"Invert cut logic");
    auto cut_dim_plot2 = NamedParameter<unsigned int>("cut_dim_plot2", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot2 = NamedParameter<double>("cut_limits_plot2", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot2(cut_dim_plot2,cut_limits_plot2,invertCut_plot2);    

    auto invertCut_plot3 = NamedParameter<bool>("invertCut_plot3", 0,"Invert cut logic");
    auto cut_dim_plot3 = NamedParameter<unsigned int>("cut_dim_plot3", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot3 = NamedParameter<double>("cut_limits_plot3", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot3(cut_dim_plot3,cut_limits_plot3,invertCut_plot3);    

    auto invertCut_plot4 = NamedParameter<bool>("invertCut_plot4", 0,"Invert cut logic");
    auto cut_dim_plot4 = NamedParameter<unsigned int>("cut_dim_plot4", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot4 = NamedParameter<double>("cut_limits_plot4", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot4(cut_dim_plot4,cut_limits_plot4,invertCut_plot4);    

    auto invertCut_plot5 = NamedParameter<bool>("invertCut_plot5", 0,"Invert cut logic");
    auto cut_dim_plot5 = NamedParameter<unsigned int>("cut_dim_plot5", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot5 = NamedParameter<double>("cut_limits_plot5", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot5(cut_dim_plot5,cut_limits_plot5,invertCut_plot5);    

    auto invertCut_plot6 = NamedParameter<bool>("invertCut_plot6", 0,"Invert cut logic");
    auto cut_dim_plot6 = NamedParameter<unsigned int>("cut_dim_plot6", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot6 = NamedParameter<double>("cut_limits_plot6", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot6(cut_dim_plot6,cut_limits_plot6,invertCut_plot6);    
    
    auto invertCut_plot7 = NamedParameter<bool>("invertCut_plot7", 0,"Invert cut logic");
    auto cut_dim_plot7 = NamedParameter<unsigned int>("cut_dim_plot7", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot7 = NamedParameter<double>("cut_limits_plot7", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot7(cut_dim_plot7,cut_limits_plot7,invertCut_plot7);    
    
    auto invertCut_plot8 = NamedParameter<bool>("invertCut_plot8", 0,"Invert cut logic");
    auto cut_dim_plot8 = NamedParameter<unsigned int>("cut_dim_plot8", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot8 = NamedParameter<double>("cut_limits_plot8", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot8(cut_dim_plot8,cut_limits_plot8,invertCut_plot8);    
    
    auto invertCut_plot9 = NamedParameter<bool>("invertCut_plot9", 0,"Invert cut logic");
    auto cut_dim_plot9 = NamedParameter<unsigned int>("cut_dim_plot9", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot9 = NamedParameter<double>("cut_limits_plot9", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot9(cut_dim_plot9,cut_limits_plot9,invertCut_plot9);    
    
    auto invertCut_plot10 = NamedParameter<bool>("invertCut_plot10", 0,"Invert cut logic");
    auto cut_dim_plot10 = NamedParameter<unsigned int>("cut_dim_plot10", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot10 = NamedParameter<double>("cut_limits_plot10", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot10(cut_dim_plot10,cut_limits_plot10,invertCut_plot10);    
    
    auto invertCut_plot11 = NamedParameter<bool>("invertCut_plot11", 0,"Invert cut logic");
    auto cut_dim_plot11 = NamedParameter<unsigned int>("cut_dim_plot11", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot11 = NamedParameter<double>("cut_limits_plot11", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot11(cut_dim_plot11,cut_limits_plot11,invertCut_plot11);    
    
    auto invertCut_plot12 = NamedParameter<bool>("invertCut_plot12", 0,"Invert cut logic");
    auto cut_dim_plot12 = NamedParameter<unsigned int>("cut_dim_plot12", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot12 = NamedParameter<double>("cut_limits_plot12", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot12(cut_dim_plot12,cut_limits_plot12,invertCut_plot12);    

    
  const std::string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName","Fit_weights.root");  
  TFile* weight_file = TFile::Open((outDir+"/"+FitWeightFileName).c_str(),"OPEN");
  weight_file->cd();
  auto weight_tree = (TTree*) weight_file->Get("DalitzEventList");
  if(weight_tree->GetEntries() != eventsMC.size()){
      cout << "ERROR inconsistent number of events" << endl;
      throw "ERROR";
  }
    
  cout << "Using data file: " << dataFile << endl;
  cout << "Using MC file: " << intFile << endl;
  cout << "Using weight file: " << FitWeightFileName << endl;

  //Dimensions to plot
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
  vector<unsigned int> m12{1,2};

  vector<vector<unsigned int>> dims{m123,m13,m23,m023,m02,m03,m01,m013,{1},{2},m012,m12};
  vector<string> labels{"m_Kpipi","m_Kpi","m_pipi","m_psipipi","m_psipi","m_psipi2","m_psiK","m_psiKpi","cosTheta","chi","m_psiKpi2","m_Kpi2"};
  
  vector<string> titles{"m(K#pi#pi) [GeV]","m(K#pi) [GeV]","m(#pi#pi) [GeV]","m(#psi(2S)#pi#pi) [GeV]","m(#psi(2S)#pi^{+}) [GeV]","m(#psi(2S)#pi^{-}) [GeV]", "m(#psi(2S)K) [GeV]","m(#psi(2S)K#pi) [GeV]"};
  titles = {"m(K^{#plus}#pi^{#plus}#pi^{#minus}) [GeV]","m(K^{#plus}#pi^{#minus}) [GeV]","m(#pi^{#plus}#pi^{#minus}) [GeV]","m(#psi(2S)#pi^{#plus}#pi^{#minus}) [GeV]","m(#psi(2S)#pi^{+}) [GeV]","m(#psi(2S)#pi^{-}) [GeV]", "m(#psi(2S)K^{#plus}) [GeV]","m(#psi(2S)K^{#plus}#pi^{#minus}) [GeV]"};

  vector<double> lim123{0.9,1.65};
  vector<double> lim13{0.6,1.45};
  vector<double> lim23{0.2,1.15};
  vector<double> lim023{4.05,4.85};
  vector<double> lim02{3.8,4.65};
  vector<double> lim03{3.8,4.65};
  vector<double> lim01{4.15,4.9};
  vector<double> lim013{4.3,5.2};
  vector<double> lim012{4.3,5.2};
  vector<double> lim12{0.6,1.45};

  if(pNames[1] != "psi(2S)0"){
      titles = {"m(K^{#plus}#pi^{#plus}#pi^{#minus}) [GeV]","m(K^{#plus}#pi^{#minus}) [GeV]","m(#pi^{#plus}#pi^{#minus}) [GeV]","m(J/#psi#pi^{#plus}#pi^{#minus}) [GeV]","m(J/#psi#pi^{+}) [GeV]","m(J/#psi#pi^{-}) [GeV]", "m(J/#psiK^{#plus}) [GeV]","m(J/#psiK^{#plus}#pi^{#minus}) [GeV]"};
      lim123 = {0.9,2.3};
      lim13={0.6,2};
      lim23={0.2,1.8};
      lim023={3.5,4.85};
      lim02={3.,4.65};
      lim03={3.,4.65};
      lim01={3.5,4.9};
      lim013={3.7,5.2};
      lim012={3.7,5.2};
      lim12={0.6,2};
  }
    
  vector<vector<double>> limits{lim123,lim13,lim23,lim023,lim02,lim03,lim01,lim013};
    
  titles.push_back("cos(#theta)");
  titles.push_back("#chi");
  limits.push_back({-1,1});
  limits.push_back({-3.141,3.141});

  titles.push_back("m(J/#psiK^{#plus}#pi^{#minus}) [GeV]");
  titles.push_back("m(K^{#plus}#pi^{#minus}) [GeV]");
  limits.push_back(lim012);
  limits.push_back(lim12);

  //Amps to plot
  auto legend = NamedParameter<string>("plot_legend", std::vector<string>() ).getVector();
  auto plot_weights = NamedParameter<string>("plot_weights", std::vector<string>(),"plot weight names" ).getVector();

  vector<string> weights{"data","weight"};
  for (int i = 0; i < plot_weights.size(); i++)weights.push_back(plot_weights[i]);
  
	vector<double> w(weights.size());	
  for(int i=1; i<weights.size();i++)weight_tree->SetBranchAddress(weights[i].c_str(),&w[i]);
 
  //Create histograms
  auto nBins = NamedParameter<Int_t>("nBins", 50, "Number of bins");
  vector<vector<TH1D*>> histo_set,histo_set_cut1,histo_set_cut2,histo_set_cut3,histo_set_cut4,histo_set_cut5,histo_set_cut6;
  vector<vector<TH1D*>> histo_set_cut7,histo_set_cut8,histo_set_cut9,histo_set_cut10,histo_set_cut11,histo_set_cut12;

  for(int i=0;i<dims.size();i++){
      histo_set.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut1.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut2.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut3.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut4.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut5.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut6.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));    
      histo_set_cut7.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut8.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut9.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut10.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut11.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut12.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));    
  }
  //Fill data
  for( auto& evt : events ){
	  for(int j=0;j<dims.size();j++){
          double val = 0;
          if(dims[j].size()>1) val = sqrt(evt.s(dims[j]));
          else if(dims[j][0] == 1) val = cosTheta(evt);
          else if(dims[j][0] == 2) val = chi(evt);
          
          histo_set[j][0]->Fill(val,evt.weight());
          if(filter_plot1(evt))histo_set_cut1[j][0]->Fill(val,evt.weight());
          if(filter_plot2(evt))histo_set_cut2[j][0]->Fill(val,evt.weight());
          if(filter_plot3(evt))histo_set_cut3[j][0]->Fill(val,evt.weight());
          if(filter_plot4(evt))histo_set_cut4[j][0]->Fill(val,evt.weight());
          if(filter_plot5(evt))histo_set_cut5[j][0]->Fill(val,evt.weight());
          if(filter_plot6(evt))histo_set_cut6[j][0]->Fill(val,evt.weight());        
          if(filter_plot7(evt))histo_set_cut7[j][0]->Fill(val,evt.weight());
          if(filter_plot8(evt))histo_set_cut8[j][0]->Fill(val,evt.weight());
          if(filter_plot9(evt))histo_set_cut9[j][0]->Fill(val,evt.weight());
          if(filter_plot10(evt)) histo_set_cut10[j][0]->Fill(val,evt.weight());
          if(filter_plot11(evt))histo_set_cut11[j][0]->Fill(val,evt.weight());
          if(filter_plot12(evt))histo_set_cut12[j][0]->Fill(val,evt.weight());        
      }
  }

  //Fill fit projections
  for(int i=0; i< eventsMC.size(); i++ ){
    Event evt(eventsMC[i]);
    weight_tree->GetEntry(i);

    for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++){
        double val = 0;
        if(dims[j].size()>1) val = sqrt(evt.s(dims[j]));
        else if(dims[j][0] == 1) val = cosTheta(evt);
        else if(dims[j][0] == 2) val = chi(evt);

        histo_set[j][k]->Fill(val,w[k]);
        
        if(filter_plot1(evt))histo_set_cut1[j][k]->Fill(val,w[k]);
        if(filter_plot2(evt))histo_set_cut2[j][k]->Fill(val,w[k]);
        if(filter_plot3(evt))histo_set_cut3[j][k]->Fill(val,w[k]);
        if(filter_plot4(evt))histo_set_cut4[j][k]->Fill(val,w[k]);
        if(filter_plot5(evt))histo_set_cut5[j][k]->Fill(val,w[k]);
        if(filter_plot6(evt))histo_set_cut6[j][k]->Fill(val,w[k]);
        if(filter_plot7(evt))histo_set_cut7[j][k]->Fill(val,w[k]);
        if(filter_plot8(evt))histo_set_cut8[j][k]->Fill(val,w[k]);
        if(filter_plot9(evt))histo_set_cut9[j][k]->Fill(val,w[k]);
        if(filter_plot10(evt))histo_set_cut10[j][k]->Fill(val,w[k]);
        if(filter_plot11(evt))histo_set_cut11[j][k]->Fill(val,w[k]);
        if(filter_plot12(evt))histo_set_cut12[j][k]->Fill(val,w[k]);
     }
  }
  //Plot
  TCanvas* c = new TCanvas();  
  TLegend leg(0.,0.,1,1,"");
	leg.SetLineStyle(0);
	leg.SetLineColor(0);
	leg.SetFillColor(0);
	leg.SetTextFont(22);
	leg.SetTextColor(1);
	leg.SetTextSize(0.075);
	leg.SetTextAlign(12);
	for(int k=2; k<weights.size();k++)leg.AddEntry(histo_set[0][k],legend[k].c_str(),"l");

  for(int j=0;j<dims.size();j++){ 
      plotHistos(histo_set[j]);
      c->Print((labels[j]+".eps").c_str());
  }
  
  c->Clear();
  c->Divide(3,2);
  for(int j=0;j<6;j++){ 
      c->cd(j+1);
      plotHistos(histo_set[j],false,1);
  }
  c->Print("plots.eps");

  c->Clear();
  c->Divide(3,2);
  for(int j=0;j<5;j++){ 
      c->cd(j+1);
      plotHistos(histo_set[j],true,1);
  }
  c->cd(6);
  leg.Draw();
  c->Print("amp_plots.eps");
  c->Print("amp_plots.C");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
    }
    c->Print("amp_plots2.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
        gPad->SetLogy(1);
    }
    c->Print("amp_plots2_log.eps");

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
    }
    c->Print("amp_plots3.eps");

    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1);
    }
    c->Print("amp_plots_cut1.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1);
    }
    c->Print("amp_plots_cut2.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1);
    }
    c->Print("amp_plots_cut3.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1);
    }
    c->Print("amp_plots_cut4.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1);
    }
    c->Print("amp_plots_cut5.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1);
    }
    c->Print("amp_plots_cut6.eps");
    

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut7[j],true,1);
    }
    c->Print("amp_plots_cut7.eps");
    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut8[j],true,1);
    }
    c->Print("amp_plots_cut8.eps");
    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut9[j],true,1);
    }
    c->Print("amp_plots_cut9.eps");
    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut10[j],true,1);
    }
    c->Print("amp_plots_cut10.eps");
    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut11[j],true,1);
    }
    c->Print("amp_plots_cut11.eps");
    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut12[j],true,1);
    }
    c->Print("amp_plots_cut12.eps");

    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1);
    }
    c->Print("amp_plots3_cut1.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1);
    }
    c->Print("amp_plots3_cut2.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1);
    }
    c->Print("amp_plots3_cut3.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1);
    }
    c->Print("amp_plots3_cut4.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1);
    }
    c->Print("amp_plots3_cut5.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1);
    }
    c->Print("amp_plots3_cut6.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut7[j],true,1);
    }
    c->Print("amp_plots3_cut7.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut8[j],true,1);
    }
    c->Print("amp_plots3_cut8.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut9[j],true,1);
    }
    c->Print("amp_plots3_cut9.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut10[j],true,1);
    }
    c->Print("amp_plots3_cut10.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut11[j],true,1);
    }
    c->Print("amp_plots3_cut11.eps");
    
    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut12[j],true,1);
    }
    c->Print("amp_plot3_cut12.eps");
    
	c->Clear();
	leg.Draw();
	c->Print("leg.eps");
    
    leg.SetNColumns(2);
    c->Clear();
    leg.Draw();
    c->Print("leg2.eps");
    
}

void makePlots3body(){
    
    std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
    std::string weightData = NamedParameter<std::string>("weightData", "weight");  
    std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
    std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
    auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
                                              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    auto bNamesMC = NamedParameter<std::string>("BranchesMC", std::vector<std::string>() ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    if(bNamesMC.size()==0)bNamesMC=bNames;
    auto pNames = NamedParameter<std::string>("EventType" , ""    
                                              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
    
    EventType evtType(pNames);
    EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData) );
    EventList eventsMC(intFile, evtType, Branches(bNamesMC), GetGenPdf(true), WeightBranch(weightMC));
    auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
    if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing data units from MeV -> GeV");
        events.transform( scale_transform );
    }
    if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing MC units from MeV -> GeV");
        eventsMC.transform( scale_transform );
    }
    
    auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
    auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
    auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();
    
    phsp_cut filter(cut_dim,cut_limits,invertCut);
    if(useFilter==1)events.filter(filter);
    if(useFilter==1)eventsMC.filter(filter);
    
    auto invertCut_plot1 = NamedParameter<bool>("invertCut_plot1", 0,"Invert cut logic");
    auto cut_dim_plot1 = NamedParameter<unsigned int>("cut_dim_plot1", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot1 = NamedParameter<double>("cut_limits_plot1", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot1(cut_dim_plot1,cut_limits_plot1,invertCut_plot1);    
    
    auto invertCut_plot2 = NamedParameter<bool>("invertCut_plot2", 0,"Invert cut logic");
    auto cut_dim_plot2 = NamedParameter<unsigned int>("cut_dim_plot2", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot2 = NamedParameter<double>("cut_limits_plot2", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot2(cut_dim_plot2,cut_limits_plot2,invertCut_plot2);    
    
    auto invertCut_plot3 = NamedParameter<bool>("invertCut_plot3", 0,"Invert cut logic");
    auto cut_dim_plot3 = NamedParameter<unsigned int>("cut_dim_plot3", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot3 = NamedParameter<double>("cut_limits_plot3", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot3(cut_dim_plot3,cut_limits_plot3,invertCut_plot3);    
    
    auto invertCut_plot4 = NamedParameter<bool>("invertCut_plot4", 0,"Invert cut logic");
    auto cut_dim_plot4 = NamedParameter<unsigned int>("cut_dim_plot4", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot4 = NamedParameter<double>("cut_limits_plot4", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot4(cut_dim_plot4,cut_limits_plot4,invertCut_plot4);    
    
    auto invertCut_plot5 = NamedParameter<bool>("invertCut_plot5", 0,"Invert cut logic");
    auto cut_dim_plot5 = NamedParameter<unsigned int>("cut_dim_plot5", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot5 = NamedParameter<double>("cut_limits_plot5", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot5(cut_dim_plot5,cut_limits_plot5,invertCut_plot5);    
    
    auto invertCut_plot6 = NamedParameter<bool>("invertCut_plot6", 0,"Invert cut logic");
    auto cut_dim_plot6 = NamedParameter<unsigned int>("cut_dim_plot6", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits_plot6 = NamedParameter<double>("cut_limits_plot6", std::vector<double>(),"cut window" ).getVector();    
    phsp_cut filter_plot6(cut_dim_plot6,cut_limits_plot6,invertCut_plot6);    
    
    const std::string FitWeightFileName = NamedParameter<std::string>("FitWeightFileName","Fit_weights.root");  
    TFile* weight_file = TFile::Open(FitWeightFileName.c_str(),"OPEN");
    weight_file->cd();
    auto weight_tree = (TTree*) weight_file->Get("DalitzEventList");
    if(weight_tree->GetEntries() != eventsMC.size()){
        cout << "ERROR inconsistent number of events" << endl;
        throw "ERROR";
    }
    
    cout << "Using data file: " << dataFile << endl;
    cout << "Using MC file: " << intFile << endl;
    cout << "Using weight file: " << FitWeightFileName << endl;
    
    //Dimensions to plot
    //EventType B0 D- K0S0 pi+ 
    vector<unsigned int> m12{1,2};
    vector<unsigned int> m02{0,2};
    vector<unsigned int> m01{0,1};
    
    vector<vector<unsigned int>> dims{m12,m02,m01};
    vector<string> labels{"m_Kspi","m_Dpi","m_DKs"};
    vector<string> titles{"m(K_{s}#pi) [GeV]","m(D#pi) [GeV]","m(DK_{s}) [GeV]"};
    
    //Limits
    vector<double> lim12{0.5,3.5};
    vector<double> lim02{1.9,5};
    vector<double> lim01{2.1,5.25};
    vector<vector<double>> limits{lim12,lim02,lim01};
    
    //Amps to plot
    auto legend = NamedParameter<string>("plot_legend", std::vector<string>() ).getVector();
    auto plot_weights = NamedParameter<string>("plot_weights", std::vector<string>(),"plot weight names" ).getVector();
    
    vector<string> weights{"data","weight"};
    for (int i = 0; i < plot_weights.size(); i++)weights.push_back(plot_weights[i]);
    
    vector<double> w(weights.size());    
    for(int i=1; i<weights.size();i++)weight_tree->SetBranchAddress(weights[i].c_str(),&w[i]);
    
    //Create histograms
    auto nBins = NamedParameter<Int_t>("nBins", 50, "Number of bins");
    vector<vector<TH1D*>> histo_set,histo_set_cut1,histo_set_cut2,histo_set_cut3,histo_set_cut4,histo_set_cut5,histo_set_cut6;
    vector<vector<TH2D*>> histo2D_set;
    for(int i=0;i<dims.size();i++){
        histo_set.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
        histo_set_cut1.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
        histo_set_cut2.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
        histo_set_cut3.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
        histo_set_cut4.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
        histo_set_cut5.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
        histo_set_cut6.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));      
        for(int j=i+1;j<dims.size();j++)histo2D_set.push_back(createHistos2D(dims[i],dims[j],labels[i],";"+titles[i]+";"+titles[j],nBins,limits[i],limits[j],weights));

    }
    //Fill data
    for( auto& evt : events )
    {
        for(int j=0;j<dims.size();j++) histo_set[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());
        if(filter_plot1(evt))for(int j=0;j<dims.size();j++) histo_set_cut1[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());
        if(filter_plot2(evt))for(int j=0;j<dims.size();j++) histo_set_cut2[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());
        if(filter_plot3(evt))for(int j=0;j<dims.size();j++) histo_set_cut3[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());
        if(filter_plot4(evt))for(int j=0;j<dims.size();j++) histo_set_cut4[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());
        if(filter_plot5(evt))for(int j=0;j<dims.size();j++) histo_set_cut5[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());
        if(filter_plot6(evt))for(int j=0;j<dims.size();j++) histo_set_cut6[j][0]->Fill(sqrt(evt.s(dims[j])),evt.weight());        

        int n = 0;
        for(int i=0;i<dims.size();i++)for(int j=i+1;j<dims.size();j++){ histo2D_set[n][0]->Fill(sqrt(evt.s(dims[i])),sqrt(evt.s(dims[j])),evt.weight()); n++;}
    }
    
    //Fill fit projections
    for(int i=0; i< eventsMC.size(); i++ )
    {
        Event evt(eventsMC[i]);
        weight_tree->GetEntry(i);
        
        for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);
        if(filter_plot1(evt))for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set_cut1[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);
        if(filter_plot2(evt))for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set_cut2[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);
        if(filter_plot3(evt))for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set_cut3[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);
        if(filter_plot4(evt))for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set_cut4[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);
        if(filter_plot5(evt))for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set_cut5[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);
        if(filter_plot6(evt))for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++) histo_set_cut6[j][k]->Fill(sqrt(evt.s(dims[j])),w[k]);

        int n = 0;
        for(int i=0;i<dims.size();i++)for(int j=i+1;j<dims.size();j++){
          for(int k=1; k<weights.size();k++){
            histo2D_set[n][k]->Fill(sqrt(evt.s(dims[i])),sqrt(evt.s(dims[j])),w[k]);
          }
          n++;
          }
    }
    
    //Plot
    TCanvas* c = new TCanvas();  
    TLegend leg(0.,0.,1,1.,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(22);
    leg.SetTextColor(1);
    leg.SetTextSize(0.1);
    leg.SetTextAlign(12);
    for(int k=2; k<weights.size();k++)leg.AddEntry(histo_set[0][k],legend[k].c_str(),"l");
    
    for(int j=0;j<dims.size();j++){ 
        plotHistos(histo_set[j]);
        c->Print((labels[j]+".eps").c_str());
    }
    
    c->Clear();
    c->Divide(3);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],false,1);
    }
    c->Print("plots.eps");
    
    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots.eps");

    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(4);
    gPad->SetLogy(0);
    leg.Draw();
    c->Print("amp_plots_log.eps");
    

    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots_cut1.eps");

    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots_cut2.eps");


    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots_cut3.eps");

    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots_cut4.eps");

    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots_cut5.eps");

    c->Clear();
    c->Divide(2,2);
    for(int j=0;j<3;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1);
    }
    c->cd(4);
    leg.Draw();
    c->Print("amp_plots_cut6.eps");

    
    c->Clear();
    leg.Draw();
    c->Print("leg.eps");

    c->Clear();
    int n = 0;
    for(int i=0;i<dims.size();i++)for(int j=i+1;j<dims.size();j++){ 
      histo2D_set[n][0]->Draw("colz");
      c->Print(("dalitz_"+to_string(i)+"_"+to_string(j)+".eps").c_str());
      histo2D_set[n][0]->Draw();
      c->Print(("dalitzScatter_"+to_string(i)+"_"+to_string(j)+".eps").c_str());

      histo2D_set[n][1]->Draw("colz");
      c->Print(("dalitzFit_"+to_string(i)+"_"+to_string(j)+".eps").c_str());
      histo2D_set[n][1]->Draw();
      c->Print(("dalitzFitScatter_"+to_string(i)+"_"+to_string(j)+".eps").c_str());

      for(int k=2; k<weights.size();k++){
        histo2D_set[n][k]->Draw("colz");
        c->Print(("dalitzFit_"+to_string(i)+"_"+to_string(j)+"_amp"+to_string(k)+".eps").c_str());
        histo2D_set[n][k]->Draw();
        c->Print(("dalitzFitScatter_"+to_string(i)+"_"+to_string(j)+"_amp"+to_string(k)+".eps").c_str());
      }
      n++; 
    }
}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();

  auto pNames = NamedParameter<std::string>("EventType" , "").getVector(); 
  std::string outDir = NamedParameter<std::string>("outDir", ".");

  if(pNames.size()==4)makePlots3body();
  else makePlots();
    
  int status = 0;
  if(outDir != ".") status &= system( ("mkdir -p " + outDir).c_str() );
  if( status != 0 ) ERROR("Building OutDir directory failed");  
  else if(outDir != "."){
        system(("mv *.eps " + outDir + "/").c_str());
  }
    
  return 0;
}
