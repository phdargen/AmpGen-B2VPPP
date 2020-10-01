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
#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif
#include "AmpGen/LHCbStyle.h"
#include "AmpGen/PolarisedSum.h"
#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

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

void readPlotFile(){

  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  
  TFile* f_plots = TFile::Open( plotFile.c_str(), "OPEN" ); f_plots->cd();
  auto s_Kpipi = (TH1D*) f_plots->Get("Data_s123");
  auto s_pipi = (TH1D*) f_plots->Get("Data_s23");
  auto s_psipipi = (TH1D*) f_plots->Get("Data_s023");
  auto s_psipi = (TH1D*) f_plots->Get("Data_s02");
  auto s_psipi2 = (TH1D*) f_plots->Get("Data_s03");
  auto s_Kpi = (TH1D*) f_plots->Get("Data_s13");

  auto s_Kpipi_fit = (TH1D*) f_plots->Get("Model_s123");
  auto s_pipi_fit = (TH1D*) f_plots->Get("Model_s23");  
  auto s_psipipi_fit = (TH1D*) f_plots->Get("Model_s023");
  auto s_psipi_fit = (TH1D*) f_plots->Get("Model_s02");
  auto s_psipi2_fit = (TH1D*) f_plots->Get("Model_s03");
  auto s_Kpi_fit = (TH1D*) f_plots->Get("Model_s13");

  s_Kpipi_fit->SetLineColor(kRed);
  s_pipi_fit->SetLineColor(kRed); 
  s_psipipi_fit->SetLineColor(kRed);
  s_psipi_fit->SetLineColor(kRed);
  s_psipi2_fit->SetLineColor(kRed);
  s_Kpi_fit->SetLineColor(kRed);

  TCanvas* c = new TCanvas("c");
  c->Divide(3,2);
  
  c->cd(1);
  s_Kpipi->DrawNormalized("e",1);
  s_Kpipi_fit->DrawNormalized("histcsame",1);
 
  c->cd(2);
  s_Kpi->DrawNormalized("e",1);
  s_Kpi_fit->DrawNormalized("histcsame",1);

  c->cd(3);
  s_pipi->DrawNormalized("e",1);
  s_pipi_fit->DrawNormalized("histcsame",1);
 
  c->cd(4);
  s_psipipi->DrawNormalized("e",1);
  s_psipipi_fit->DrawNormalized("histcsame",1);

  c->cd(5);
  s_psipi->DrawNormalized("e",1);
  s_psipi_fit->DrawNormalized("histcsame",1);
 
  c->cd(6);
  s_psipi2->DrawNormalized("e",1);
  s_psipi2_fit->DrawNormalized("histcsame",1);

  c->Print("plots.eps");

  f_plots->Close();

}

vector<TH1D*> createHistos(vector<unsigned int> dim,string name, string title, int nBins, vector<double> limits, vector<string> weights){

  vector<TH1D*> histos;
  TH1D* histo = new TH1D(name.c_str(),"",nBins,limits[0],limits[1]);
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

void plotHistos(vector<TH1D*>histos,bool plotComponents = true, int style = 0){

  if(style == 1){
    histos[0]->SetMarkerSize(.1);
    histos[0]->GetYaxis()->SetTitle("");
  }
  if(plotComponents == true && style == 0)histos[0]->SetMaximum(histos[0]->GetMaximum()*1.2);
  histos[0]->DrawNormalized("",1);

  //for (int i = (plotComponents == true ? histos.size()-1 : 1); i > 0; i--)
  for (int i = 1; i < (plotComponents == true ? histos.size() : 2); i++)
  {
      if(style == 1)histos[i]->SetLineWidth(1);
      double norm = histos[i]->Integral()/histos[1]->Integral();
      histos[i]->DrawNormalized("histcsame",norm);
  }
  histos[0]->DrawNormalized("same",1);
}

void makePlots(){

  std::string dataFile = NamedParameter<std::string>("DataSample", "", "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample","", "Name of file containing events to use for MC integration.");
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  EventType evtType(pNames);
  EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch("weight") );
  EventList eventsMC(intFile, evtType, Branches(bNames), GetGenPdf(true), WeightBranch("weight"));

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
  //EventType B+ psi(2S)0 K+ pi+ pi- 
  vector<unsigned int> m123{1,2,3};
  vector<unsigned int> m13{1,3};
  vector<unsigned int> m23{2,3};
  vector<unsigned int> m023{0,2,3};
  vector<unsigned int> m02{0,2};
  vector<unsigned int> m03{0,3};
  vector<unsigned int> m01{0,1};
  vector<unsigned int> m013{0,1,3};

  vector<vector<unsigned int>> dims{m123,m13,m23,m023,m02,m03,m01,m013};
  vector<string> labels{"m_Kpipi","m_Kpi","m_pipi","m_psipipi","m_psipi","m_psipi2","m_psiK","m_psiKpi"};
  vector<string> titles{"m(K#pi#pi) [GeV]","m(K#pi) [GeV]","m(#pi#pi) [GeV]","m(#psi(2S)#pi#pi) [GeV]","m(#psi(2S)#pi^{+}) [GeV]","m(#psi(2S)#pi^{-}) [GeV]", "m(#psi(2S)K) [GeV]","m(#psi(2S)K#pi) [GeV]"};

  //Limits
  vector<double> lim123{0.9,1.7};
  vector<double> lim13{0.5,1.5};
  vector<double> lim23{0.2,1.2};
  vector<double> lim023{4,4.85};
  vector<double> lim02{3.8,4.7};
  vector<double> lim03{3.8,4.7};
  vector<double> lim01{4.1,4.9};
  vector<double> lim013{4.25,5.25};

  vector<vector<double>> limits{lim123,lim13,lim23,lim023,lim02,lim03,lim01,lim013};

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
  for(int i=0;i<dims.size();i++){
      histo_set.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut1.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut2.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut3.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut4.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut5.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));
      histo_set_cut6.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));      
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
	leg.Draw();
	c->Print("leg.eps");
}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  makePlots();

  return 0;
}
