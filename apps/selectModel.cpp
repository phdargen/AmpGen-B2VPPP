#include <chrono>
#include <ctime>
#include <iostream>
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
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
using namespace std;
using namespace AmpGen;

void selectModel(){

  std::string resultsFile = NamedParameter<std::string>("ResultsFile"     , "results.root", "Name of the output plot file");
  
  TChain* chain = new TChain("Result");
  chain->Add((replaceAll(resultsFile,".root","_*.root")).c_str());  
 
  double nll,chi2,sumFractions;
  unsigned status, nPar,nAmps, seed;
    
  chain->SetBranchAddress("nll",&nll);
  chain->SetBranchAddress("chi2",&chi2);
  chain->SetBranchAddress("sumFractions",&sumFractions);
  chain->SetBranchAddress("status",&status);
  chain->SetBranchAddress("nPar",&nPar);
  chain->SetBranchAddress("nAmps",&nAmps);
  chain->SetBranchAddress("seed",&seed);
        
  chain->GetEntry(0);
  double min_nll = nll; 
  double max_nll = nll; 
  double min_chi2 = chi2;   
  int min_i = 0;  
  int min_nPar = nPar;  
    
  const int n = chain->GetEntries();
  for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
    
        if(TMath::IsNaN(nll))continue;
      
        if(nll>max_nll)max_nll=nll;
      
        //if(status>0)continue;
        if(nll<min_nll){ 
            min_nll=nll;
            min_i=i;
            min_nPar = nPar;
        }
        if(chi2<min_chi2)min_chi2=chi2;      
  }
    
  TH1D * h_nll = new TH1D("h_nll","h_nll",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll0 = new TH1D("h_nll1","h_nll1",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll1 = new TH1D("h_nll1","h_nll1",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll2 = new TH1D("h_nll2","h_nll2",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll3 = new TH1D("h_nll3","h_nll3",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll4 = new TH1D("h_nll4","h_nll4",40,-0.5,min(max_nll-min_nll,3000.5));  
    
  h_nll0->SetLineColor(kRed);  
  h_nll1->SetLineColor(kBlue);  
  h_nll2->SetLineColor(kGreen);  
  h_nll3->SetLineColor(kYellow);  
  h_nll4->SetLineColor(kOrange);  

  TH1D * h_chi2 = new TH1D("h_chi2","h_chi2",50,0,10);  
  TH1D * h_chi20 = new TH1D("h_chi20","h_chi20",50,0,10);  
  TH1D * h_chi21 = new TH1D("h_chi21","h_chi21",50,0,10);  
  TH1D * h_chi22 = new TH1D("h_chi22","h_chi22",50,0,10);  
  TH1D * h_chi23 = new TH1D("h_chi23","h_chi23",50,0,10);  
  TH1D * h_chi24 = new TH1D("h_chi24","h_chi24",50,0,10);  

  h_chi20->SetLineColor(kRed);  
  h_chi21->SetLineColor(kBlue);  
  h_chi22->SetLineColor(kGreen);  
  h_chi23->SetLineColor(kYellow);  
  h_chi24->SetLineColor(kOrange);  
  
  for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
        if(TMath::IsNaN(nll))continue;

        h_nll->Fill(nll-min_nll);
        h_chi2->Fill(chi2);
        
        double AIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) ;
        double BIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) * log(30000) ;
      
        cout << endl << "model " << i << endl;
        cout << "n2ll = " << nll-min_nll << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "status = " << status << endl;
        cout << "nAmps = " << nAmps << endl;
        cout << "nPar = " << nPar << endl;
        cout << "AIC = " << AIC << " ( " << sqrt( AIC ) << " sigma) " << endl;
        cout << "BIC = " << BIC << " ( " << sqrt( BIC ) << " sigma) " << endl;
        //cout << "sigma = " << TMath::ErfcInverse(TMath::Prob(nll-min_nll,max(1.,abs((double)min_nPar-(double)nPar)))) * sqrt(2) << endl << endl;       
      
      if(status==0){
            h_nll0->Fill(nll-min_nll);
            h_chi20->Fill(chi2);            
        }
      if(status==1){
          h_nll1->Fill(nll-min_nll);
          h_chi21->Fill(chi2);            
      }
      if(status==2){
          h_nll2->Fill(nll-min_nll);
          h_chi22->Fill(chi2);            
      }
      if(status==3){
          h_nll3->Fill(nll-min_nll);
          h_chi23->Fill(chi2);            
      }
      if(status==4){
          h_nll4->Fill(nll-min_nll);
          h_chi24->Fill(chi2);            
      }
  }
    
  TCanvas* c = new TCanvas("c");
  //c->Divide(3,2);

  h_nll->Draw("h");
  h_nll0->Draw("hsame");
  h_nll1->Draw("hsame");
  h_nll2->Draw("hsame");
  h_nll3->Draw("hsame");
  h_nll4->Draw("hsame");
  c->Print("n2ll.eps");  

  h_chi2->Draw("h");
  h_chi20->Draw("hsame");
  h_chi21->Draw("hsame");
  h_chi22->Draw("hsame");
  h_chi23->Draw("hsame");
  h_chi24->Draw("hsame");
  c->Print("chi2.eps");
    
  chain->GetEntry(min_i);
  cout << endl << "Best model: " << min_i << endl;
  cout << "with seed: " << seed << endl;
  cout << "n2ll = " << nll << endl;
  cout << "chi2 = " << chi2 << endl;
  cout << "status = " << status << endl;
  cout << "nAmps = " << nAmps << endl;
  cout << "nPar = " << nPar << endl;

}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  selectModel();

  return 0;
}
