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
#include <TGraph.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
using namespace std;
using namespace AmpGen;

void bubble_sort(double* arr, double* arr2, double* arr3, double* arr4, size_t len)
{
    double temp, temp2, temp3, temp4;
    for (int i = 0; i < len; i++)
        {
            for (int j = 0; j + 1 < len - i; j++)
            {
                // Swaping the elements if first one
                // is greater than second one.
                if (arr[j] > arr[j + 1])
                {
                    temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;

                    temp2 = arr2[j];
                    arr2[j] = arr2[j + 1];
                    arr2[j + 1] = temp2;
                    
                    temp3 = arr3[j];
                    arr3[j] = arr3[j + 1];
                    arr3[j + 1] = temp3;

                    temp4 = arr4[j];
                    arr4[j] = arr4[j + 1];
                    arr4[j + 1] = temp4;
                }
            }
        }
}

void selectModel(){

  std::string resultsFile = NamedParameter<std::string>("ResultsFile"     , "results.root", "Name of the output plot file");
  auto scanParam = NamedParameter<std::string>("scanParam", std::vector<std::string>()).getVector(); // name, min, max, nSteps

  string outDir = resultsFile;
  outDir = replaceAll(outDir,"result.root","");
  
  TChain* chain = new TChain("Result");
  chain->Add((replaceAll(resultsFile,".root","_*.root")).c_str());  
 
  double nll,chi2,sumFractions, nSig;
  unsigned status, nPar,nAmps, seed;
    
  chain->SetBranchAddress("nll",&nll);
  chain->SetBranchAddress("chi2",&chi2);
  chain->SetBranchAddress("Sum_Bp",&sumFractions);
  chain->SetBranchAddress("status",&status);
  chain->SetBranchAddress("nPar",&nPar);
  chain->SetBranchAddress("nAmps",&nAmps);
  chain->SetBranchAddress("nSig",&nSig);
  chain->SetBranchAddress("seed",&seed);
        
  chain->GetEntry(0);
  double min_nll = nll; 
  double max_nll = nll; 
  double min_chi2 = chi2;   
  double max_chi2 = chi2;   
  int min_i = 0;  
  int min_nPar = nPar;  
  int min_chi2_i = 0;  
    
  const int n = chain->GetEntries();
  for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
    
        if(TMath::IsNaN(nll))continue;
      
        if(nll>max_nll)max_nll=nll;
        if(chi2>max_chi2)max_chi2=chi2;

        //if(status>0)continue;
        if(nll<min_nll){ 
            min_nll=nll;
            min_i=i;
            min_nPar = nPar;
        }
        if(chi2<min_chi2){
            min_chi2=chi2;    
            min_chi2_i=i;
        }
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

  TH1D * h_chi2 = new TH1D("h_chi2","h_chi2",50,0,max_chi2*1.1);  
  TH1D * h_chi20 = new TH1D("h_chi20","h_chi20",50,0,max_chi2*1.1);  
  TH1D * h_chi21 = new TH1D("h_chi21","h_chi21",50,0,max_chi2*1.1);  
  TH1D * h_chi22 = new TH1D("h_chi22","h_chi22",50,0,max_chi2*1.1);  
  TH1D * h_chi23 = new TH1D("h_chi23","h_chi23",50,0,max_chi2*1.1);  
  TH1D * h_chi24 = new TH1D("h_chi24","h_chi24",50,0,max_chi2*1.1);  

  h_chi20->SetLineColor(kRed);  
  h_chi21->SetLineColor(kBlue);  
  h_chi22->SetLineColor(kGreen);  
  h_chi23->SetLineColor(kYellow);  
  h_chi24->SetLineColor(kOrange);  
    
  double nll_vals[n];
  double sig_vals[n];
  double chi2_vals[n]; 
  double seed_vals[n]; 
  
  for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
        if(TMath::IsNaN(nll))continue;

        h_nll->Fill(nll-min_nll);
        h_chi2->Fill(chi2);
      
        nll_vals[i]=nll-min_nll;
        sig_vals[i]= TMath::ErfcInverse(TMath::Prob(nll-min_nll,max(1.,abs((double)min_nPar-(double)nPar)))) * sqrt(2);
        chi2_vals[i]=chi2;
        if(scanParam.size()==4)seed_vals[i]= std::stod(scanParam[1]) + ( std::stod(scanParam[2])-std::stod(scanParam[1]) ) * (double)seed/(double)std::stoi(scanParam[3]);
        else seed_vals[i]= seed;

        double AIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) ;
        double BIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) * log(nSig) ;
      
        cout << "seed " << seed << endl;
        if(scanParam.size()==4) cout << "seed val " << seed_vals[i] << endl;
        cout << "n2ll = " << nll-min_nll << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "status = " << status << endl;
        cout << "nAmps = " << nAmps << endl;
        cout << "nPar = " << nPar << endl;
        cout << "Sum of fractions = " << sumFractions << endl;
        cout << "AIC = " << AIC << " ( " << sqrt( AIC ) << " sigma) " << endl;
        cout << "BIC = " << BIC << " ( " << sqrt( BIC ) << " sigma) " << endl;
        cout << "sigma = " << TMath::ErfcInverse(TMath::Prob(nll-min_nll,max(1.,abs((double)min_nPar-(double)nPar)))) * sqrt(2) << endl << endl;       
      
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
  c->Print((outDir+"n2ll.pdf").c_str());  

  h_chi2->Draw("h");
  h_chi20->Draw("hsame");
  h_chi21->Draw("hsame");
  h_chi22->Draw("hsame");
  h_chi23->Draw("hsame");
  h_chi24->Draw("hsame");
  c->Print((outDir+"chi2.pdf").c_str());
    
  chain->GetEntry(min_i);
  cout << endl << "Best nll model: " << endl;
  cout << "with seed: " << seed << endl;
  if(scanParam.size()==4) cout << "seed val " << seed_vals[min_i] << endl;
  cout << "n2ll = " << nll << endl;
  cout << "chi2 = " << chi2 << endl;
  cout << "status = " << status << endl;
  cout << "nAmps = " << nAmps << endl;
  cout << "nPar = " << nPar << endl;
    
  chain->GetEntry(min_chi2_i);
  cout << endl << "Best chi2 model: " << endl;
  cout << "with seed: " << seed << endl;
  if(scanParam.size()==4) cout << "seed val " << seed_vals[min_chi2_i] << endl;
  cout << "n2ll = " << nll << endl;
  cout << "chi2 = " << chi2 << endl;
  cout << "status = " << status << endl;
  cout << "nAmps = " << nAmps << endl;
  cout << "nPar = " << nPar << endl;
    
  bubble_sort(seed_vals,nll_vals,sig_vals,chi2_vals,n);

  TGraph* g_nll = new TGraph(n,seed_vals,nll_vals);
  TGraph* g_chi2 = new TGraph(n,seed_vals,chi2_vals);
  TGraph* g_sig = new TGraph(n,seed_vals,sig_vals);

  if(scanParam.size()==4)g_nll->GetXaxis()->SetTitle(scanParam[0].c_str());  
  else g_nll->GetXaxis()->SetTitle("seed"); 
  g_nll->GetYaxis()->SetTitle("n2ll"); 
  g_nll->SetMarkerColor(kRed);
  g_nll->SetMarkerSize(1);
    
  if(scanParam.size()==4)g_sig->GetXaxis()->SetTitle(scanParam[0].c_str());  
  else g_sig->GetXaxis()->SetTitle("seed"); 
  g_sig->GetYaxis()->SetTitle("#sigma_{n2ll}"); 
  g_sig->SetMarkerColor(kRed);
  g_sig->SetMarkerSize(1);

  if(scanParam.size()==4)g_chi2->GetXaxis()->SetTitle(scanParam[0].c_str());  
  else g_chi2->GetXaxis()->SetTitle("seed"); 
  g_chi2->GetYaxis()->SetTitle("chi2"); 
  g_chi2->SetMarkerColor(kRed);
  g_chi2->SetMarkerSize(1.);
    
  g_nll->Draw("AC*");
  c->Print((outDir+"n2ll_seed.pdf").c_str());  
    
  g_sig->Draw("AC*");
  c->Print((outDir+"sigma_seed.pdf").c_str());  

  g_chi2->Draw("AC*");
  c->Print((outDir+"chi2_seed.pdf").c_str());  
    
}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  selectModel();

  return 0;
}
