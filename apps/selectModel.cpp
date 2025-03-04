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
#include "AmpGen/FitResult.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
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

  std::string resultsFile = NamedParameter<std::string>("ResultsFile"     , "result.root", "Name of the output plot file");
  auto scanParam = NamedParameter<std::string>("scanParam", std::vector<std::string>()).getVector(); // name, min, max, nSteps
  bool isScan = (scanParam.size()==4);
  bool isParam = (scanParam.size()==3);

  string outDir = resultsFile;
  outDir = replaceAll(outDir,"result.root","");
  
  TChain* chain = new TChain("Result");
  chain->Add((replaceAll(resultsFile,".root","_*.root")).c_str());  
 
  double nll, chi2, sumFractions, nSig, par;
  unsigned status, nPar,nAmps, seed;
    
  chain->SetBranchAddress("nll",&nll);
  chain->SetBranchAddress("chi2",&chi2);
  chain->SetBranchAddress("Sum_Bp",&sumFractions);
  chain->SetBranchAddress("status",&status);
  chain->SetBranchAddress("nPar",&nPar);
  chain->SetBranchAddress("nAmps",&nAmps);
  chain->SetBranchAddress("nSig",&nSig);
  chain->SetBranchAddress("seed",&seed);
  if(isParam)chain->SetBranchAddress(scanParam[0].c_str(),&par);

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
        nll*=2;
    
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
    
  TH1D * h_nll = new TH1D("h_nll","; #Delta(-2log L); # Fits",60,-0.5,min(max_nll-min_nll,120.5));
  TH1D * h_nll0 = new TH1D("h_nll1","h_nll1",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll1 = new TH1D("h_nll1","h_nll1",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll2 = new TH1D("h_nll2","h_nll2",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll3 = new TH1D("h_nll3","h_nll3",40,-0.5,min(max_nll-min_nll,3000.5));  
  TH1D * h_nll4 = new TH1D("h_nll4","h_nll4",40,-0.5,min(max_nll-min_nll,3000.5));  
   
  h_nll->SetLineColor(kBlue);
  h_nll0->SetLineColor(kRed);  
  //h_nll1->SetLineColor(kBlue);
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
    
  TH1D * h_par = new TH1D("h_par",("; " + scanParam[0] + "; # Fits").c_str(),60, std::stod(scanParam[1]), std::stod(scanParam[2])  );
    
  double nll_vals[n];
  double sig_vals[n];
  double chi2_vals[n]; 
  double seed_vals[n]; 
  
  for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
        nll*=2;
        if(TMath::IsNaN(nll))continue;

        h_nll->Fill(nll-min_nll);
        h_chi2->Fill(chi2);
      
        nll_vals[i]=nll-min_nll;
        sig_vals[i]= TMath::ErfcInverse(TMath::Prob(nll-min_nll,max(1.,abs((double)min_nPar-(double)nPar)))) * sqrt(2);
        chi2_vals[i]=chi2;
      
        if(isScan)seed_vals[i]= std::stod(scanParam[1]) + ( std::stod(scanParam[2])-std::stod(scanParam[1]) ) * (double)seed/(double)std::stoi(scanParam[3]);
        else if(isParam)seed_vals[i]= par;
        else seed_vals[i]= seed;
      
        h_par->Fill(seed_vals[i]);

        double AIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) ;
        double BIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) * log(nSig) ;
      
        cout << "seed " << seed << endl;
        if(isScan || isParam) cout << "seed val " << seed_vals[i] << endl;
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
  //h_nll0->Draw("hsame");
  //h_nll1->Draw("hsame");
  //h_nll2->Draw("hsame");
  //h_nll3->Draw("hsame");
  //h_nll4->Draw("hsame");
  c->Print((outDir+"n2ll.pdf").c_str());  

  h_chi2->Draw("h");
  h_chi20->Draw("hsame");
  h_chi21->Draw("hsame");
  h_chi22->Draw("hsame");
  h_chi23->Draw("hsame");
  h_chi24->Draw("hsame");
  c->Print((outDir+"chi2.pdf").c_str());
    
  h_par->Draw("h");
  c->Print((outDir+"h_par_" + scanParam[0] + ".pdf").c_str());
    
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
  g_nll->SetMaximum(120);
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
  c->Print((outDir+"scan_n2ll_seed.pdf").c_str());
    
  g_sig->Draw("AC*");
  c->Print((outDir+"scan_sigma_seed.pdf").c_str());

  g_chi2->Draw("AC*");
  c->Print((outDir+"scan_chi2_seed.pdf").c_str());
    
  cout << "Not converged = " << h_nll1->GetEntries() + h_nll2->GetEntries() + h_nll3->GetEntries() + h_nll4->GetEntries() <<  endl;
    cout << "Overflow = " << h_nll->GetBinContent(61) << endl;
    
}

void scan(){

  string resultsFile = NamedParameter<std::string>("ResultsFile", "result.root", "");
  string resultsBaseFile = NamedParameter<std::string>("ResultsBaseFile", "", "");
  if(resultsBaseFile=="")resultsBaseFile=resultsFile;
    
  auto scanParam = NamedParameter<std::string>("scanParam", std::vector<std::string>()).getVector(); // name, min, max, nSteps
  bool isScan = (scanParam.size()==4);
   
  auto nSteps = std::stoi(scanParam[3]);
  auto nScan = NamedParameter<int>("nScan", nSteps);

  TCanvas* c = new TCanvas("c");
  //c->Divide(3,2);
  string outDir = resultsFile;
  outDir = replaceAll(outDir,"result.root","");
    
  double min_nll_all = 0;
  vector<TGraph*> vec_g_nll;
  vector<TGraph*> vec_g_chi2;
  vector<TGraph*> vec_g_sig;
  
  for (int k=0; k<nScan/nSteps; k++) {

        int n = 0;
        TChain* chain = new TChain("Result");

        if(std::ifstream(((string)resultsBaseFile).c_str()).good()){
            chain->Add(resultsBaseFile.c_str());
            n++;
        }
        else cout << "ERROR::Baseline file not found" << endl;
      
        for(int j = nSteps * k ; j < nSteps * (k+1); j++){
            string fileName = replaceAll(resultsFile,".root","_"+ to_string(j) +".root");
            if(std::ifstream((fileName).c_str()).good()){
                chain->Add((fileName).c_str());
                n++;
                cout << "added file = " << fileName << " ; k = " << k << " ; j = " << j << endl;
            }
        }
      
        double nll, chi2, sumFractions, nSig;
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
        
        double nll_baseline = nll;
        double chi2_baseline = chi2;
        int nPar_baseline = nPar;
        
        //int n = chain->GetEntries();
        for (int i = 1; i<n; i++) {
            chain->GetEntry(i);
            
            if(TMath::IsNaN(nll))continue;
            if(nll>nll_baseline)continue;
            if(status>0)continue;
            
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
        TH1D * h_chi2 = new TH1D("h_chi2","h_chi2",50,0,max_chi2*1.1);
        
        double nll_vals[n];
        double sig_vals[n];
        double chi2_vals[n];
        double seed_vals[n];
        double seed_vals2[n];

        for (int i = 1; i<n; i++) {
            chain->GetEntry(i);
            
            if(TMath::IsNaN(nll) || status>0 || nll>nll_baseline){
                chain->GetEntry(i+1);
                double tmp_nll = nll;
                double tmp_chi2 = chi2;
                
                chain->GetEntry(i-1);
                tmp_nll = (nll + tmp_nll) / 2.;
                tmp_chi2 = (chi2 + tmp_chi2) / 2.;
                
                chain->GetEntry(i);
                nll  = tmp_nll;
                chi2 = tmp_chi2;
            }
                        
            min_nll_all = nll-nll_baseline > min_nll_all ? min_nll_all : nll-nll_baseline;
            h_nll->Fill(nll-nll_baseline);
            h_chi2->Fill(chi2);
            
            nll_vals[i]=nll-nll_baseline;
            sig_vals[i]= TMath::ErfcInverse(TMath::Prob(abs(nll-nll_baseline),max(1.,abs((double)nPar_baseline-(double)nPar)))) * sqrt(2);
            chi2_vals[i]=chi2;
            
            if(isScan)seed_vals[i]= std::stod(scanParam[1]) + ( std::stod(scanParam[2])-std::stod(scanParam[1]) ) * (double)(seed%nSteps)/(double)nSteps;
            else seed_vals[i]= seed;
            seed_vals2[i]= seed;
            
            double AIC = (nll-nll_baseline) + 2. * ((double) nPar - (double) nPar_baseline) ;
            double BIC = (nll-nll_baseline) + 2. * ((double) nPar - (double) nPar_baseline) * log(nSig) ;
            
            /*
            cout << "seed " << seed << endl;
            if(scanParam.size()==4) cout << "seed val " << seed_vals[i] << endl;
            cout << "n2ll = " << nll-nll_baseline << endl;
            cout << "chi2 = " << chi2 << endl;
            cout << "status = " << status << endl;
            cout << "nAmps = " << nAmps << endl;
            cout << "nPar = " << nPar << endl;
            cout << "Sum of fractions = " << sumFractions << endl;
            cout << "AIC = " << AIC << " ( " << sqrt( AIC ) << " sigma) " << endl;
            cout << "BIC = " << BIC << " ( " << sqrt( BIC ) << " sigma) " << endl;
            cout << "sigma = " << TMath::ErfcInverse(TMath::Prob(abs(nll-nll_baseline),max(1.,abs((double)nPar_baseline-(double)nPar)))) * sqrt(2) << endl << endl;
            */
        }
               
        for (int i = 2; i<n; i++){
            
            if( abs( nll_vals[i] - nll_vals[i-1] ) > 20 ){
            //if( abs( nll_vals[i] - nll_vals[i-1] ) > 50   ||  abs( nll_vals[i] - (nll_vals[i-1]+nll_vals[i+1])/2. ) > 50 ){
                cout << "Outlier:: nll= " << nll_vals[i] << " ; seed = " << seed_vals[i] << " ; seed2 = " << seed_vals2[i] << endl;
            }
            
            if( nll_vals[i] > -1 ){
            //if( abs( nll_vals[i] - nll_vals[i-1] ) > 50   ||  abs( nll_vals[i] - (nll_vals[i-1]+nll_vals[i+1])/2. ) > 50 ){
                cout << "Small:: nll= " << nll_vals[i] << " ; seed = " << seed_vals[i] << " ; seed2 = " << seed_vals2[i] << endl;
            }
            
            //if( abs( nll_vals[i] - (nll_vals[i-1]+nll_vals[i+1])/2. ) > 50 ) nll_vals[i] = (nll_vals[i-1]+nll_vals[i+1])/2.;
            
            //if( i > 3 && abs( nll_vals[i] - (nll_vals[i-1]+nll_vals[i+1] + nll_vals[i-2]+nll_vals[i+2] )/4. ) > 30 ) nll_vals[i] = (nll_vals[i-1]+nll_vals[i+1] + nll_vals[i-2]+nll_vals[i+2] )/4.;

        }
      
        h_nll->Draw("h");
        c->Print((outDir+"n2ll.pdf").c_str());
        
        h_chi2->Draw("h");
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
        
        TString title(scanParam[0]);
        title.ReplaceAll("Xs","X_{s}");
        title.ReplaceAll("Zs","Z_{s}");
        title.ReplaceAll("_mass"," mass [MeV]");
        title.ReplaceAll(")0",")^{0}");

        if(scanParam.size()==4)g_nll->GetXaxis()->SetTitle(title);
        else g_nll->GetXaxis()->SetTitle("seed");
        if(scanParam.size()==4)g_nll->GetXaxis()->SetLimits(std::stod(scanParam[1]),std::stod(scanParam[2]));
        g_nll->GetYaxis()->SetTitle("#Delta(-2log L)");
        if(k==0)g_nll->SetMarkerColor(kRed);
        if(k==2)g_nll->SetMarkerColor(kBlue);
        if(k==0)g_nll->SetLineColor(kRed);
        if(k==2)g_nll->SetLineColor(kBlue);
        if(k==0)g_nll->SetMarkerStyle(25);
        if(k==2)g_nll->SetMarkerStyle(26);
        g_nll->SetMarkerSize(1);

        if(scanParam.size()==4)g_sig->GetXaxis()->SetTitle(title);
        else g_sig->GetXaxis()->SetTitle("seed");
        if(scanParam.size()==4)g_sig->GetXaxis()->SetLimits(std::stod(scanParam[1]),std::stod(scanParam[2]));
        g_sig->GetYaxis()->SetTitle("#sigma_{n2ll}");
        g_sig->SetMarkerColor(kRed);
        g_sig->SetMarkerSize(1);
        
        if(scanParam.size()==4)g_chi2->GetXaxis()->SetTitle(title);
        else g_chi2->GetXaxis()->SetTitle("seed");
        if(scanParam.size()==4)g_chi2->GetXaxis()->SetLimits(std::stod(scanParam[1]),std::stod(scanParam[2]));
        g_chi2->GetYaxis()->SetTitle("chi2");
        g_chi2->SetMarkerColor(kRed);
        g_chi2->SetMarkerSize(1.);
        
        g_nll->Draw("APC");
        c->Print((outDir+"scan_n2ll_seed_"+to_string(k)+".pdf").c_str());
        
        g_sig->Draw("AC*");
        c->Print((outDir+"scan_sigma_seed_"+to_string(k)+".pdf").c_str());
        
        g_chi2->Draw("AC*");
        c->Print((outDir+"scan_chi2_seed_"+to_string(k)+".pdf").c_str());
      
        //TGraphSmooth *gs = new TGraphSmooth();
        //auto g_nll_s = gs->SmoothKern(g_nll);
        vec_g_nll.push_back(g_nll);
    }
    
    vec_g_nll[0]->SetMinimum(min_nll_all);
    vec_g_nll[0]->Draw("APL");
    for(auto& g : vec_g_nll)g->Draw("PLsame");
    c->Print((outDir+"scan_n2ll_seed.pdf").c_str());
    
}

void analyzeFits(){
    INFO("Analyze fit results ... ");

    string baseLineFit = NamedParameter<string>("baseLineFit", "log.txt");
    string outDir = NamedParameter<string>("outDir", "out/");

    FitResult fr_base(baseLineFit);
    fr_base.print();
    
    auto fitFiles = NamedParameter<string>("fitFiles", vector<string>()).getVector(); // basename, first file number, last file number

    vector<FitResult> fits;
    for(unsigned int i = stoi(fitFiles[1]) ; i <= stoi(fitFiles[2]) ; ++i){
      if(! std::ifstream((fitFiles[0] + to_string(i)+".txt").c_str()).good()){
          ERROR("File " << fitFiles[0] + to_string(i)+".txt" << " not found");
          continue;
      }
      FitResult fr(fitFiles[0] + to_string(i)+".txt");
      fits.push_back(fr);
      
      //fr.print();
      double dLL = fr.LL() - fr_base.LL();
      //double LL_sig = sqrt(2) * TMath::ErfcInverse(TMath::Prob(abs(dLL),1));
      //if(abs(dLL)>0 && TMath::Prob(abs(dLL),1) == 0 ) LL_sig = 100;
      //if(dLL>0) LL_sig *= -1.;
      int dParams = fr.floating().size() - fr_base.floating().size();
      double dChi2 = fr.chi2()/((double)fr.nBins()-1.-fr.nParam()) - fr_base.chi2()/((double)fr_base.nBins()-1.-fr_base.nParam()) ;
      double chi2 = fr.chi2()/((double)fr.nBins()-1.);
      INFO("Fit " << i << ":: dLL = " << dLL << " chi2/bin = " << chi2 << " dChi2 = " << dChi2 << " dParams =  " << dParams << endl );
    }
    
}


void fitChiSquareDistributionForSignificanceTests() {
    
    auto nFits = NamedParameter<int>("fitChi2SignTests::nFits", 1);
    auto nBins = NamedParameter<int>("fitChi2SignTests::nBins", 100);
    auto histMin = NamedParameter<double>("fitChi2SignTests::histMin", -10);
    auto histMax = NamedParameter<double>("fitChi2SignTests::histMax", 100);
    auto fitMin = NamedParameter<double>("fitChi2SignTests::fitMin", -10);
    auto fitMax = NamedParameter<double>("fitChi2SignTests::fitMax", 100);
    auto sign = NamedParameter<int>("fitChi2SignTests::sign", 1);
    auto fitGauss = NamedParameter<int>("fitChi2SignTests::fitGauss", 0);

    string outDir = NamedParameter<string>("outDir", "");
    TH1D* hist = new TH1D("hist",";#Delta(2lnL); # Fits",nBins,histMin,histMax);
    
    for(unsigned int i = 0 ; i < nFits ; ++i){
        if(! std::ifstream((outDir+"log_" + to_string(i)+".txt").c_str()).good()){
            ERROR("File " << outDir+"log_" + to_string(i)+".txt" << " not found");
            continue;
        }
        FitResult fr0(outDir+"log_H0_" + to_string(i)+".txt");
        FitResult fr(outDir+"log_" + to_string(i)+".txt");

        double dLL = sign * ( fr0.LL() - fr.LL() );
        
        hist->Fill(dLL);
  }
    
  // Create a TF1 function for the chi-square distribution
  TF1* func;
  if(fitGauss) func =  new TF1("func", "[0]*ROOT::Math::normal_pdf(x, [2], [1])", fitMin, fitMax);
  else func = new TF1("func", "[0]*ROOT::Math::chisquared_pdf(x, [1])", fitMin, fitMax);
    
  // Set some initial parameters for the fit
  if(fitGauss){
        func->SetParameters(hist->GetEntries(), hist->GetMean(), hist->GetRMS());
        func->SetParName(0, "Normalization");
        func->SetParName(1, "Mean");
        func->SetParName(2, "Sigma");
        func->SetParLimits(2, 0.1, hist->GetRMS() * 10);

  }
  else{
        func->SetParameters(hist->GetEntries(), 1.0);
        func->SetParName(0, "Normalization");
        func->SetParName(1, "Degrees of Freedom");
  }
      
  // Fit the histogram with the chi-square distribution
  hist->Fit(func, "R"); // "R" option for adjusting the fit range
  
  // Print the fit result
  TF1* fitFunction = hist->GetFunction("func");
  if (fitFunction) {
    double normalization = fitFunction->GetParameter(0);
    printf("Fit Result:\n");
    printf("Normalization: %f\n", normalization);
    
    if(fitGauss){
        double meanFit = fitFunction->GetParameter(1);
        double sigmaFit = fitFunction->GetParameter(2);
        printf("Mean: %f\n", meanFit);
        printf("Sigma: %f\n", sigmaFit);
    }
    else{
        double degreesOfFreedom = fitFunction->GetParameter(1);
        printf("Degrees of Freedom: %f\n", degreesOfFreedom);
    }
  }
  
  // Create a canvas to display the histogram and the fit function
  TCanvas* canvas = new TCanvas("canvas", "Chi-Square Distribution Fit", 800, 600);
  
  // Draw the histogram
  hist->SetLineColor(kBlue);
  hist->Draw("hist");
  
  // Draw the fit function
  func->SetLineColor(kRed);
  func->Draw("same");
  
  // Update the canvas
  canvas->Update();
  
  // Save the canvas as a PDF file
  canvas->Print(((string)outDir+"dLL_fit.pdf").c_str());
}


int main( int argc, char* argv[] ){

  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
    
  auto mode = NamedParameter<int>("selectModel::mode", 1);

  if(mode==0)selectModel();
  if(mode==1)scan();
  if(mode==2)analyzeFits();
  if(mode==3)fitChiSquareDistributionForSignificanceTests();

  return 0;
}
