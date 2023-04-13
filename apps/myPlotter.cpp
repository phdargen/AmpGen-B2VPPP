#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <sstream>      

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
#include <TStatistic.h>

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

inline double cosHel(const Event& evt, int i, int j){

   TLorentzVector particle = pFromEvent( evt, i );
   TLorentzVector parent = pFromEvent( evt, i ) + pFromEvent( evt, j );
   TLorentzVector grandparent = pFromEvent( evt, 0 ) + pFromEvent( evt, 1 ) + pFromEvent( evt, 2 ) + pFromEvent( evt, 3 ) + pFromEvent( evt, 4 );  
    
   TVector3 boosttoparent = -(parent.BoostVector());

   particle.Boost(boosttoparent);
   grandparent.Boost(boosttoparent);

   TVector3 particle3 = particle.Vect();
   TVector3 grandparent3 = grandparent.Vect();
   double numerator = particle3.Dot(grandparent3);
   double denominator = (particle3.Mag())*(grandparent3.Mag());
   double temp = numerator/denominator;

   return temp;
}


inline double helicityCos(const Event& evt, int i, int j, vector<int> kl){
    
    TLorentzVector pi = pFromEvent( evt, i ); 
    TLorentzVector pj = pFromEvent( evt, j ); 
    TLorentzVector PR = pFromEvent( evt, kl[0] ) + pFromEvent( evt, kl[1] );
    
    return dotProduct(pi, pj, PR) / sqrt( dotProduct( pi, pi, PR ) * dotProduct( pj, pj, PR ) );    
}

inline double cosTheta( const Event& evt ){

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

inline double chi( const Event& evt ){
    
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
        
    //return cosChi;
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


vector<TH1D*> createHistos(vector<unsigned int> dim,string name, string title, int nBins, vector<double> limits, vector<string> weights){

  vector<TH1D*> histos;
  TH1D* histo = new TH1D(name.c_str(),"",dim.size()==1 ? nBins : nBins*2,limits[0],limits[1]);
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
            h->SetLineWidth(2);
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
            h->SetLineColor(kOrange-7);
            h->SetLineWidth(3);
            //h->SetFillColor(kGray+3);
            //h->SetFillStyle(1001);
      }else if(i==8){
          h->SetLineColor(kOrange+1);
          h->SetLineWidth(3);
          //h->SetFillColor(kGray+3);
          //h->SetFillStyle(1001);
      }else if(i==weights.size()-1){
          h->SetLineColor(kGray);
          h->SetLineWidth(3);
          h->SetFillColor(kGray);
          //h->SetFillStyle(1001);
      }else {
          h->SetLineColor(kBlack);
          h->SetLineWidth(3);
          h->SetLineStyle(kDashed);
      }
      histos.push_back(h);
  }
    
  unsigned int nPermErrorBands   = NamedParameter<unsigned int>("nPermErrorBands",0);  
  for(int i = 0; i < nPermErrorBands; i++)  {
      TH1D* h = (TH1D*) histo->Clone((name+"Err"+to_string(i)).c_str());
      histos.push_back(h);
  }
    
  unsigned int addSysErrBand   = NamedParameter<unsigned int>("addSysErrBand",0);    
  if(addSysErrBand)for(int i = 0; i < nPermErrorBands; i++)  {
        TH1D* h = (TH1D*) histo->Clone((name+"ErrTot"+to_string(i)).c_str());
        histos.push_back(h);
  }  
    
  auto FitWeightFileNames = NamedParameter<std::string>("AltFitWeightFileNames", std::vector<std::string>(),"AltFitWeightFileNames" ).getVector();
  for(int i = 0; i < FitWeightFileNames.size(); i++)  {
        TH1D* h = (TH1D*) histo->Clone((name+"Alt"+to_string(i)).c_str());
        histos.push_back(h);
  }
    
  return histos;
}

void plotHistos(vector<TH1D*>histos, bool plotComponents = true, int style = 0, bool computeErrorBands = false){
  
  if(style == 0){
      histos[0]->SetLineWidth(1);
      histos[0]->SetMarkerSize(1.25);
      histos[0]->GetYaxis()->SetTitle("");
  }
  if(style == 1){
    histos[0]->SetLineWidth(1);
    histos[0]->SetMarkerSize(.5);
    histos[0]->GetYaxis()->SetTitle("");
  }
  histos[0]->SetMinimum(1);
  if(plotComponents == true && style == 0)histos[0]->SetMaximum(histos[0]->GetMaximum()*1.3);
  histos[0]->DrawNormalized("",1);
  
  //Compute err bands for fit
  unsigned int nPermErrorBands   = NamedParameter<unsigned int>("nPermErrorBands",0);  
  auto FitWeightFileNames = NamedParameter<std::string>("AltFitWeightFileNames", std::vector<std::string>(),"AltFitWeightFileNames" ).getVector();
  unsigned int nAltModels = FitWeightFileNames.size();
  unsigned int addSysErrBand   = NamedParameter<unsigned int>("addSysErrBand",0);  
  unsigned int plotAltModels   = NamedParameter<unsigned int>("plotAltModels",0);

//  INFO("nPermErrorBands = " << nPermErrorBands);  
//  INFO("nAltModels = " << nAltModels);  
//  INFO("addSysErrBand = " << addSysErrBand);  

  // Plot alt models
  if(plotAltModels){
        histos[1]->SetLineColor(kRed);
        histos[1]->SetMarkerColor(kRed);
        histos[histos.size()-1]->SetLineColor(kAzure - 4);
        histos[histos.size()-1]->SetMarkerColor(kAzure - 4);
        histos[histos.size()-1]->SetFillColorAlpha(kAzure - 4, 1);
        histos[histos.size()-1]->SetMarkerSize(0.);
        histos[histos.size()-1]->DrawNormalized("histsame",1);
//        for (unsigned int i=histos.size()-1 ; i >= histos.size()-nAltModels; i--) {
//            histos[i]->SetLineColor(i);
//            histos[i]->SetMarkerColor(i);
//            histos[i]->DrawNormalized("e1same",1);
//        }
        computeErrorBands = false;
  }
    
  // Only compute once for every histo set !  
  if(computeErrorBands){
      
      auto h_norm = (TH1D*) histos[1]->Clone(((string)histos[1]->GetName()+"norm").c_str());
      h_norm->Scale(1./h_norm->Integral());
      
      for (unsigned int i=histos.size()-1 ; i >= histos.size()-nPermErrorBands*(1+addSysErrBand)-nAltModels; i--){
                histos[i]->Scale(1./histos[i]->Integral());
      }  
        
      // Stat error bars
      if(nPermErrorBands>3){
          
//          for (unsigned int i=histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand ; i >= histos.size()-nPermErrorBands*(1+addSysErrBand)-nAltModels; i--) {
//              histos[i]->SetLineColor(i);                    
//              histos[i]->SetMarkerColor(i);
//              histos[i]->DrawNormalized("e1same",1);
//          }
          
            for(unsigned int b = 1; b <= histos[0]->GetXaxis()->GetNbins(); b++){
                
                TStatistic stat;
                double min = TMath::Limits<Double_t>::Max();
                double max = 0;
                
                for (unsigned int i=histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand ; i >= histos.size()-nPermErrorBands*(1+addSysErrBand)-nAltModels; i--) {
                    stat.Fill(histos[i]->GetBinContent(b)) ;
                    min = (histos[i]->GetBinContent(b) < min) ? histos[i]->GetBinContent(b) : min;
                    max = (histos[i]->GetBinContent(b) > max) ? histos[i]->GetBinContent(b) : max;

                    //histos[i]->SetLineColor(i);                    
                    //histos[i]->SetMarkerColor(i);
                    //histos[i]->DrawNormalized("e1same",1);
                }
                //if(b==20)stat.Print();
                histos[histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand]->SetBinContent(b, stat.GetMean());
                histos[histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand]->SetBinError(b,  (stat.GetMean() > 0 && stat.GetRMS() > 0) > 0 ? stat.GetRMS() : 0.);
                
                histos[histos.size()-2-nAltModels-nPermErrorBands*addSysErrBand]->SetBinContent(b, stat.GetMean());
                histos[histos.size()-2-nAltModels-nPermErrorBands*addSysErrBand]->SetBinError(b,  stat.GetMean() > 0 ? abs(max-min)/2. : 0.);
                
                histos[histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand]->SetMarkerSize(0.);
                histos[histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand]->SetFillColor( kAzure - 4);

                histos[histos.size()-2-nAltModels-nPermErrorBands*addSysErrBand]->SetMarkerSize(0.);
                histos[histos.size()-2-nAltModels-nPermErrorBands*addSysErrBand]->SetFillColor( kAzure - 4);

                if(addSysErrBand){
                    TStatistic stat_tot;
                    min = TMath::Limits<Double_t>::Max();
                    max = 0;
                    
//                    for (unsigned int i=histos.size()-1-nAltModels ; i >= histos.size()-nPermErrorBands-nAltModels; i--) {
//                        histos[i]->SetLineColor(i);                    
//                        histos[i]->SetMarkerColor(i);
//                        histos[i]->DrawNormalized("e1same",1);
//                    }
                    
                    for (unsigned int i=histos.size()-1-nAltModels ; i >= histos.size()-nPermErrorBands-nAltModels; i--) {
                        stat_tot.Fill(histos[i]->GetBinContent(b)) ;
                        min = (histos[i]->GetBinContent(b) < min) ? histos[i]->GetBinContent(b) : min;
                        max = (histos[i]->GetBinContent(b) > max) ? histos[i]->GetBinContent(b) : max;
                    }
                    //if(b==20)stat_tot.Print();

                    double stat_err = nPermErrorBands>3 ? histos[histos.size()-1-nAltModels-nPermErrorBands*addSysErrBand]->GetBinError(b) : 0.;
                    
                    histos[histos.size()-1-nAltModels]->SetBinContent(b, stat_tot.GetMean());
                    histos[histos.size()-1-nAltModels]->SetBinError(b,  (stat_tot.GetMean() > 0 && stat_tot.GetRMS() > 0) > 0 ? stat_tot.GetRMS() : 0);
                    
                    histos[histos.size()-2-nAltModels]->SetBinContent(b, stat_tot.GetMean());
                    histos[histos.size()-2-nAltModels]->SetBinError(b,  stat_tot.GetMean() > 0 ? abs(max-min)/2. : stat_err);
                    
                    histos[histos.size()-1-nAltModels]->SetMarkerSize(0.);
                    histos[histos.size()-1-nAltModels]->SetFillColor(kAzure - 4);

                    histos[histos.size()-2-nAltModels]->SetMarkerSize(0.);
                    histos[histos.size()-2-nAltModels]->SetFillColor(kAzure - 4);
                }
                
            }
        

      }  
        
      //Compute err bands from alternative fit
      if(nAltModels>0){
          
              for(unsigned int b = 1; b <= histos[0]->GetXaxis()->GetNbins(); b++){
                  
                  TStatistic stat;
                  double min = TMath::Limits<Double_t>::Max();
                  double max = 0;
                  
                  for (unsigned int i=histos.size()-1 ; i >= histos.size()-nAltModels; i--) {
                      stat.Fill(histos[i]->GetBinContent(b)) ;
                      min = (histos[i]->GetBinContent(b) < min) ? histos[i]->GetBinContent(b) : min;
                      max = (histos[i]->GetBinContent(b) > max) ? histos[i]->GetBinContent(b) : max;
                  }
                  
                  double stat_err = nPermErrorBands>3 ? histos[histos.size()-2-nAltModels]->GetBinError(b) : 0.;
                  
                  histos[histos.size()-1]->SetBinContent(b, stat.GetMean());
                  histos[histos.size()-1]->SetBinError(b,  (stat.GetMean() > 0 && stat.GetRMS() > 0) ? sqrt( pow(stat.GetRMS(),2) + pow(stat_err,2) ) : stat_err);
                  
                  histos[histos.size()-2]->SetBinContent(b, stat.GetMean());
                  //histos[histos.size()-2]->SetBinError(b,  stat.GetMean() > 0 ? sqrt( pow(abs(max-min)/2.,2) + pow(stat_err,2) ) : stat_err);
                  histos[histos.size()-2]->SetBinError(b,  stat.GetMean() > 0 ? abs(max-min)/2. + stat_err  : stat_err);

              }
          
              histos[histos.size()-1]->SetMarkerSize(0.);
              histos[histos.size()-1]->SetFillColor(862);
              histos[histos.size()-1]->SetFillColorAlpha(kAzure - 4, 1);

              histos[histos.size()-2]->SetMarkerSize(0.);
              histos[histos.size()-2]->SetFillColor(862);
              histos[histos.size()-2]->SetFillColorAlpha(kAzure - 4, 1);        
      } 
      
      // Plot error bands
      if(nAltModels>0){
          histos[histos.size()-2]->DrawNormalized("e5same",1);
      }
      if(nPermErrorBands>3){
            if(addSysErrBand)histos[histos.size()-2-nAltModels]->DrawNormalized("e5same",1);
            histos[histos.size()-2-nAltModels-nPermErrorBands*addSysErrBand]->DrawNormalized("e5same",1);
      }
          
  }

  // Plot fit projections
  for (unsigned int i = (plotComponents == true ? histos.size()-nPermErrorBands-nAltModels-1 : 1); i >= 1 ; i--)
  {
      histos[i]->SetMinimum(1);
      histos[i]->SetLineWidth(2);
      if(style == 0 && i==1)histos[i]->SetLineWidth(4);
      if(nPermErrorBands>3)histos[1]->SetLineWidth(1);
      double norm = histos[i]->Integral()/histos[1]->Integral();
      histos[i]->DrawNormalized("histcsame",norm);
  }
  histos[0]->DrawNormalized("same",1);
  gPad->RedrawAxis();
}

void plotData(vector<TH1D*>histos,bool plotComponents = false, int style = 0){
    if(style == 1){
        histos[0]->SetLineWidth(1);
        histos[0]->SetMarkerSize(.5);
        histos[0]->GetYaxis()->SetTitle("");
    }
    histos[0]->SetMinimum(1);
    histos[0]->DrawNormalized("",1);
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
    std::string weightData = NamedParameter<std::string>("weightData", "");  
    if(weightData==" " || weightData=="noWeight")  weightData="";
    std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
    std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
    auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    auto bNamesMC = NamedParameter<std::string>("BranchesMC", std::vector<std::string>() ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    if(bNamesMC.size()==0)bNamesMC=bNames;
    auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

    auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
    auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
    auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
    auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();
    phsp_cut filter(cut_dim,cut_limits,invertCut);

    EventType evtType(pNames);
    vector<size_t> entryList,entryListMC;
    if(useFilter==1){
        EventList dummyEvents(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));
        EventList dummyEventsMC(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
        
        if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") dummyEvents.transform( scale_transform );      
        if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") dummyEventsMC.transform( scale_transform );
        
        for(int i=0; i< dummyEvents.size(); i++ )if(filter(dummyEvents[i]))entryList.push_back(i);
        for(int i=0; i< dummyEventsMC.size(); i++ )if(filter(dummyEventsMC[i]))entryListMC.push_back(i);
    }

    EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData),EntryList(entryList) );
    EventList eventsMC(intFile, evtType, Branches(bNamesMC), GetGenPdf(true), WeightBranch(weightMC),EntryList(entryListMC));
    if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing data units from MeV -> GeV");
        events.transform( scale_transform );
    }
    if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing MC units from MeV -> GeV");
        eventsMC.transform( scale_transform );
    }

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
    vector<string> titles{"#it{m(K^{#plus}#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(K^{#plus}#pi^{#minus})} [GeV]","#it{m(#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(#psi(2S)#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(#psi(2S)#pi^{+})} [GeV]","#it{m(#psi(2S)#pi^{#minus})} [GeV]", "#it{m(#psi(2S)K^{#plus})} [GeV]","#it{m(#psi(2S)K^{#plus}#pi^{#minus})} [GeV]"};

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
      titles = {"#it{m(K^{#plus}#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(K^{#plus}#pi^{#minus})} [GeV]","#it{m(#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(J/#psi#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(J/#psi#pi^{+})} [GeV]","#it{m(J/#psi#pi^{#minus})} [GeV]", "#it{m(J/#psiK^{#plus})} [GeV]","#it{m(J/#psiK^{#plus}#pi^{#minus})} [GeV]"};
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

    if(pNames[1] != "psi(2S)0")titles.push_back("#it{m(J/#psiK^{#plus}#pi^{#minus})} [GeV]");
    else titles.push_back("#it{m(#psi}(2S)it{K^{#plus}#pi^{#minus})} [GeV]");
    titles.push_back("#it{m(K^{#plus}#pi^{#minus})} [GeV]");
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
    weight_tree->GetEntry(i);

    for(int j=0;j<dims.size();j++)for(int k=1; k<weights.size();k++){
        double val = 0;
        if(dims[j].size()>1) val = sqrt(eventsMC[i].s(dims[j]));
        else if(dims[j][0] == 1) val = cosTheta(eventsMC[i]);
        else if(dims[j][0] == 2) val = chi(eventsMC[i]);

        histo_set[j][k]->Fill(val,w[k]);
        
        if(filter_plot1(eventsMC[i]))histo_set_cut1[j][k]->Fill(val,w[k]);
        if(filter_plot2(eventsMC[i]))histo_set_cut2[j][k]->Fill(val,w[k]);
        if(filter_plot3(eventsMC[i]))histo_set_cut3[j][k]->Fill(val,w[k]);
        if(filter_plot4(eventsMC[i]))histo_set_cut4[j][k]->Fill(val,w[k]);
        if(filter_plot5(eventsMC[i]))histo_set_cut5[j][k]->Fill(val,w[k]);
        if(filter_plot6(eventsMC[i]))histo_set_cut6[j][k]->Fill(val,w[k]);
        if(filter_plot7(eventsMC[i]))histo_set_cut7[j][k]->Fill(val,w[k]);
        if(filter_plot8(eventsMC[i]))histo_set_cut8[j][k]->Fill(val,w[k]);
        if(filter_plot9(eventsMC[i]))histo_set_cut9[j][k]->Fill(val,w[k]);
        if(filter_plot10(eventsMC[i]))histo_set_cut10[j][k]->Fill(val,w[k]);
        if(filter_plot11(eventsMC[i]))histo_set_cut11[j][k]->Fill(val,w[k]);
        if(filter_plot12(eventsMC[i]))histo_set_cut12[j][k]->Fill(val,w[k]);
     }
    }

    //Plot
    TCanvas* c = new TCanvas();  
    TLegend leg(0.,0.,1.0,1,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(132);
    leg.SetTextColor(1);
    leg.SetTextSize(0.065);
    leg.SetTextAlign(12);
    for(int k=2; k<weights.size();k++)leg.AddEntry(histo_set[0][k],legend[k].c_str(),"l");
    
    //Get chi2 for legend
    string resultsFileName = NamedParameter<string>("ResultsFile","result.root");
    TFile* results_file = TFile::Open((outDir+"/"+resultsFileName).c_str(),"OPEN");
    results_file->cd();
    auto result_tree = (TTree*) results_file->Get("Result");
    
    double nll,chi2,sumFractions;
    int status, nPar,nAmps, seed;
    result_tree->SetBranchAddress("nll",&nll);
    result_tree->SetBranchAddress("chi2",&chi2);
    result_tree->SetBranchAddress("status",&status);
    result_tree->SetBranchAddress("nPar",&nPar);
    result_tree->SetBranchAddress("nAmps",&nAmps);
    result_tree->SetBranchAddress("seed",&seed);
    result_tree->GetEntry(0);
    stringstream ss_chi;
    ss_chi << setprecision(3) << chi2;
    leg.AddEntry((TObject *)0,("#chi^{2}/ndf = " + ss_chi.str()).c_str() ,"");
    leg.AddEntry((TObject *)0,"" ,"");
    leg.AddEntry((TObject *)0,("nPar = " + to_string(nPar)).c_str() ,"");
    leg.AddEntry((TObject *)0,("nAmps = " + to_string(nAmps)).c_str() ,"");
    leg.SetNColumns(2);

//    for(int j=0;j<dims.size();j++){ 
//        plotData(histo_set[j]);
//        c->Print(("data_"+labels[j]+".eps").c_str());
//    }
    
    for(int j=0;j<dims.size();j++){ 
      plotHistos(histo_set[j]);
      c->Print((outDir+"/"+labels[j]+".eps").c_str());
    }

    c->Clear();
    c->Divide(3,2);
    for(int j=0;j<6;j++){ 
      c->cd(j+1);
      plotHistos(histo_set[j],false,1);
    }
    c->Print((outDir+"/"+"plots.eps").c_str());

    c->Clear();
    c->Divide(3,2);
    for(int j=0;j<5;j++){ 
      c->cd(j+1);
      plotHistos(histo_set[j],true,1);
    }
    c->cd(6);
    leg.Draw();
    c->Print((outDir+"/"+"amp_plots.eps").c_str());
    c->Print((outDir+"/"+"amp_plots.C").c_str());
    c->Print((outDir+"/"+"amp_plots.pdf").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots2.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
        gPad->SetLogy(1);
    }
    c->Print((outDir+"/"+"amp_plots2_log.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3.eps").c_str());
    c->Print((outDir+"/"+"amp_plots3.pdf").c_str());


    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut1.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut2.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut3.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut4.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut5.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut6.eps").c_str());


    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut7[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut7.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut8[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut8.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut9[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut9.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut10[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut10.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut11[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut11.eps").c_str());

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut12[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots_cut12.eps").c_str());


    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut1.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut2.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut3.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut4.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut5.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut6.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut7[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut7.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut8[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut8.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut9[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut9.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut10[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut10.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut11[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plots3_cut11.eps").c_str());

    c->Clear();
    c->Divide(5,2);
    for(int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut12[j],true,1);
    }
    c->Print((outDir+"/"+"amp_plot3_cut12.eps").c_str());

    c->Clear();
    leg.Draw();
    c->Print((outDir+"/"+"leg.eps").c_str());

    leg.SetNColumns(2);
    c->Clear();
    leg.Draw();
    c->Print((outDir+"/"+"leg2.eps").c_str());

}

void makePlotsMuMu(){
    
    std::string decayMuMu = NamedParameter<std::string>("decayMuMu", "psi(2S)0");
    std::string outDir = NamedParameter<std::string>("outDir", ".");
    std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
    std::string weightData = NamedParameter<std::string>("weightData", "weight");  
    if(weightData==" " || weightData=="noWeight")  weightData="";
    std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
    std::string weightMC = NamedParameter<std::string>("weightMC", "weight");
    auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>(),"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
    auto bNamesMC = NamedParameter<std::string>("BranchesMC", std::vector<std::string>() ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
   // if(bNamesMC.size()==0)bNamesMC=bNames;
    auto pNames = NamedParameter<std::string>("EventType" , "", "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
    
    auto scale_transform = [](auto& event){ for( size_t x = 0 ; x < event.size(); ++x ) event[x] /= 1000.; };
    auto useFilter = NamedParameter<Int_t>("useFilter", 0,"Apply phsp cut");
    auto invertCut = NamedParameter<bool>("invertCut", 0,"Invert cut logic");
    auto cut_dim = NamedParameter<unsigned int>("cut_dim", std::vector<unsigned int>(),"dimension to cut on" ).getVector();
    auto cut_limits = NamedParameter<double>("cut_limits", std::vector<double>(),"cut window" ).getVector();
    phsp_cut filter(cut_dim,cut_limits,invertCut);
    
    EventType evtType(pNames);
    vector<size_t> entryList,entryListMC;
    if(useFilter==1){
        EventList dummyEvents(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData));
        EventList dummyEventsMC(intFile, evtType, Branches(bNamesMC), WeightBranch(weightMC), GetGenPdf(true));
        
        if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") dummyEvents.transform( scale_transform );      
        if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") dummyEventsMC.transform( scale_transform );
        
        for(unsigned int i=0; i< dummyEvents.size(); i++ )if(filter(dummyEvents[i]))entryList.push_back(i);
        for(unsigned int i=0; i< dummyEventsMC.size(); i++ )if(filter(dummyEventsMC[i]))entryListMC.push_back(i);
    }
    
    EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weightData),EntryList(entryList) );
    EventList eventsMC(intFile, evtType, Branches(bNamesMC), GetGenPdf(true), WeightBranch(weightMC),EntryList(entryListMC));
    if( NamedParameter<std::string>("DataUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing data units from MeV -> GeV");
        events.transform( scale_transform );
    }
    if( NamedParameter<std::string>("MCUnits", "GeV").getVal()  == "MeV") {
        INFO("Changing MC units from MeV -> GeV");
        eventsMC.transform( scale_transform );
    }
        
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
    auto FitWeightFileNames = NamedParameter<std::string>("AltFitWeightFileNames", std::vector<std::string>(),"AltFitWeightFileNames" ).getVector();
    unsigned int nAltModels = FitWeightFileNames.size();
    FitWeightFileNames.insert(FitWeightFileNames.begin(),outDir+"/"+FitWeightFileName);
    
    vector<TTree*> weight_trees;
    for( const auto& f : FitWeightFileNames ){
        TFile* weight_file = TFile::Open(f.c_str(),"OPEN");
        weight_file->cd();
        auto weight_tree = (TTree*) weight_file->Get("DalitzEventList");
        if(weight_tree->GetEntries() != eventsMC.size()){
            cout << "ERROR inconsistent number of events" << endl;
            throw "ERROR";
        }
        weight_trees.push_back(weight_tree);
    }
    
    cout << "Using data file: " << dataFile << endl;
    cout << "Using MC file: " << intFile << endl;
    cout << "Using weight file: " << FitWeightFileNames[0] << endl;
    
    //Dimensions to plot
    //EventType B+ K+ pi+ pi- mu+ mu-
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

    vector<vector<unsigned int>> dims{m012,m02,m12,m3412,m341,m342,m340,m3402,{1},{2},m3401,m01,m34,{3}};
    vector<string> labels{"m_Kpipi","m_Kpi","m_pipi","m_psipipi","m_psipi","m_psipi2","m_psiK","m_psiKpi","cosTheta","chi","m_psiKpi2","m_Kpi2","m_mumu","cosThetaKs"};
    vector<string> titles{"#it{m(K^{#plus}#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(K^{#plus}#pi^{#minus})} [GeV]","#it{m(#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(#psi(2S)#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(#psi(2S)#pi^{+})} [GeV]","#it{m(#psi(2S)#pi^{#minus})} [GeV]", "#it{m(#psi(2S)K^{#plus})} [GeV]","#it{m(#psi(2S)K^{#plus}#pi^{#minus})} [GeV]"};
    
    vector<double> lim012{0.9,1.65};
    vector<double> lim02{0.6,1.45};
    vector<double> lim12{0.2,1.15};
    vector<double> lim3412{4.05,4.85};
    vector<double> lim341{3.8,4.65};
    vector<double> lim342{3.8,4.65};
    vector<double> lim340{4.15,4.9};
    vector<double> lim3402{4.3,5.2};
    vector<double> lim3401{4.3,5.2};
    vector<double> lim01{0.6,1.45};
    vector<double> lim34{3.68609-0.002*10,3.68609+0.002*10};
    
    if(decayMuMu != "psi(2S)0"){
        titles = {"#it{m(K^{#plus}#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(K^{#plus}#pi^{#minus})} [GeV]","#it{m(#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(J/#psi#pi^{#plus}#pi^{#minus})} [GeV]","#it{m(J/#psi#pi^{+})} [GeV]","#it{m(J/#psi#pi^{#minus})} [GeV]", "#it{m(J/#psiK^{#plus})} [GeV]","#it{m(J/#psiK^{#plus}#pi^{#minus})} [GeV]"};
        lim012 = {0.9,2.3};
        lim02={0.6,2};
        lim12={0.2,1.8};
        lim3412={3.5,4.85};
        lim341={3.,4.65};
        lim342={3.,4.65};
        lim340={3.5,4.9};
        lim3402={3.7,5.2};
        lim3401={3.7,5.2};
        lim01={0.6,2};
    }
    
    vector<vector<double>> limits{lim012,lim02,lim12,lim3412,lim341,lim342,lim340,lim3402};
    
    titles.push_back("cos(#theta)");
    titles.push_back("#chi");
    limits.push_back({-1,1});
    limits.push_back({-3.141,3.141});
    
    if(decayMuMu != "psi(2S)0")titles.push_back("#it{m(J/#psiK^{#plus}#pi^{#minus})} [GeV]");
    else titles.push_back("#it{m(#psi}(2S)#it{K^{#plus}#pi^{+})} [GeV]");
    titles.push_back("#it{m(K^{#plus}#pi^{#minus})} [GeV]");
    titles.push_back("#it{m(#psi}(2S))} [GeV]");

    limits.push_back(lim3401);
    limits.push_back(lim01);
    limits.push_back(lim34);
    
    titles.push_back("cos(#theta_K*)");
    limits.push_back({-1,1});
    
    // Plot fit
    auto legend = NamedParameter<string>("plot_legend", std::vector<string>() ).getVector();
    vector<string> weights{"data","weight"};

    // Amps to plot    
    auto plot_weights = NamedParameter<string>("plot_weights", std::vector<string>(),"plot weight names" ).getVector();
    for (unsigned int i = 0; i < plot_weights.size(); i++)weights.push_back(plot_weights[i]);

    // Plot bkg?
    auto f_bkg = NamedParameter<double>("f_bkg", std::vector<double>(),"f_bkg" ).getVector();
    if(f_bkg.size()>1 )if(f_bkg[1] != 0.0){
        //weights.push_back("weight_sig");
        weights.push_back("weight_bkg");
        //legend.push_back("Signal");
        legend.push_back("Background");
    }
    
    vector<double> w(weights.size());    
    for(unsigned int i=1; i<weights.size();i++)weight_trees[0]->SetBranchAddress(weights[i].c_str(),&w[i]);
    
    // Err bands 
    unsigned int nPermErrorBands   = NamedParameter<unsigned int>("nPermErrorBands",0);  
    vector<double> wErr(nPermErrorBands);    
    for(unsigned int i = 0; i < nPermErrorBands; i++) weight_trees[0]->SetBranchAddress( ("weightErr_"+to_string(i)).c_str(),&wErr[i]); 

    unsigned int addSysErrBand   = NamedParameter<unsigned int>("addSysErrBand",0);  
    vector<double> wErrTot(nPermErrorBands);    
    if(addSysErrBand)for(unsigned int i = 0; i < nPermErrorBands; i++) weight_trees[0]->SetBranchAddress( ("weightErrTot_"+to_string(i)).c_str(),&wErrTot[i]); 
    
    // Alternative models
    vector<double> wAlt(nAltModels);    
    for(unsigned int i = 0; i < nAltModels; i++) weight_trees[i+1]->SetBranchAddress( "weight",&wAlt[i]); 

    //Create histograms
    auto nBins = NamedParameter<Int_t>("nBins", 50, "Number of bins");
    vector<vector<TH1D*>> histo_set,histo_set_cut1,histo_set_cut2,histo_set_cut3,histo_set_cut4,histo_set_cut5,histo_set_cut6;
    vector<vector<TH1D*>> histo_set_cut7,histo_set_cut8,histo_set_cut9,histo_set_cut10,histo_set_cut11,histo_set_cut12;
    vector<vector<TH1D*>> histo_set_cutCombo1,histo_set_cutCombo2,histo_set_cutCombo3,histo_set_cutCombo4;

    for(unsigned int i=0;i<dims.size();i++){
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
        histo_set_cutCombo1.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));  
        histo_set_cutCombo2.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));  
        histo_set_cutCombo3.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));  
        histo_set_cutCombo4.push_back(createHistos(dims[i],labels[i],titles[i],nBins,limits[i],weights));  
    }
    //Fill data
    auto kstar_hcos = HelicityCosine({3},{0,1,2},{3,4}); 
    
    cout << "Fill data hists ..." << endl;
    for( auto& evt : events ){
        for(unsigned int j=0;j<dims.size();j++){
            double val = 0;
            if(dims[j].size()>1) val = sqrt(evt.s(dims[j]));
            else if(dims[j][0] == 1) val = cosThetaMuAngle(evt);
            else if(dims[j][0] == 2) val = chiMuAngle(evt);
            else if(dims[j][0] == 3) val =  cosHel(evt,2,0); // -kstar_hcos(evt);

            //evt.print();
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
            
            if(filter_plot1(evt) && filter_plot6(evt))  histo_set_cutCombo1[j][0]->Fill(val,evt.weight());   
            if(filter_plot2(evt) && filter_plot5(evt))  histo_set_cutCombo2[j][0]->Fill(val,evt.weight());   
            if(filter_plot3(evt) && filter_plot5(evt))  histo_set_cutCombo3[j][0]->Fill(val,evt.weight());   
            if(filter_plot4(evt) && filter_plot5(evt))  histo_set_cutCombo4[j][0]->Fill(val,evt.weight());   
        }
    }
    
    //Fill fit projections
    cout << "Fill fits hists ..." << endl;
    auto nHists = weights.size();
    
    for(unsigned int i=0; i< eventsMC.size(); i++ ){
        
        for( const auto& weight_tree : weight_trees )weight_tree->GetEntry(i);
        
        for(unsigned int j=0;j<dims.size();j++){

            double val = 0;
            if(dims[j].size()>1) val = sqrt(eventsMC[i].s(dims[j]));
            else if(dims[j][0] == 1) val = cosThetaMuAngle(eventsMC[i]);
            else if(dims[j][0] == 2) val = chiMuAngle(eventsMC[i]);
            else if(dims[j][0] == 3) val = cosHel(eventsMC[i],2,0); // -kstar_hcos(eventsMC[i]);
            
            for(unsigned int k=1; k < nHists; k++){

                histo_set[j][k]->Fill(val,w[k]);
                if(filter_plot1(eventsMC[i]))histo_set_cut1[j][k]->Fill(val,w[k]);
                if(filter_plot2(eventsMC[i]))histo_set_cut2[j][k]->Fill(val,w[k]);
                if(filter_plot3(eventsMC[i]))histo_set_cut3[j][k]->Fill(val,w[k]);
                if(filter_plot4(eventsMC[i]))histo_set_cut4[j][k]->Fill(val,w[k]);
                if(filter_plot5(eventsMC[i]))histo_set_cut5[j][k]->Fill(val,w[k]);
                if(filter_plot6(eventsMC[i]))histo_set_cut6[j][k]->Fill(val,w[k]);
                if(filter_plot7(eventsMC[i]))histo_set_cut7[j][k]->Fill(val,w[k]);
                if(filter_plot8(eventsMC[i]))histo_set_cut8[j][k]->Fill(val,w[k]);
                if(filter_plot9(eventsMC[i]))histo_set_cut9[j][k]->Fill(val,w[k]);
                if(filter_plot10(eventsMC[i]))histo_set_cut10[j][k]->Fill(val,w[k]);
                if(filter_plot11(eventsMC[i]))histo_set_cut11[j][k]->Fill(val,w[k]);
                if(filter_plot12(eventsMC[i]))histo_set_cut12[j][k]->Fill(val,w[k]);
                
                if(filter_plot1(eventsMC[i]) && filter_plot6(eventsMC[i]))histo_set_cutCombo1[j][k]->Fill(val,w[k]);   
                if(filter_plot2(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo2[j][k]->Fill(val,w[k]);   
                if(filter_plot3(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo3[j][k]->Fill(val,w[k]);   
                if(filter_plot4(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo4[j][k]->Fill(val,w[k]);   
            }
            
            for(unsigned int k = 0; k < nPermErrorBands; k++){
                histo_set[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot1(eventsMC[i]))histo_set_cut1[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot2(eventsMC[i]))histo_set_cut2[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot3(eventsMC[i]))histo_set_cut3[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot4(eventsMC[i]))histo_set_cut4[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot5(eventsMC[i]))histo_set_cut5[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot6(eventsMC[i]))histo_set_cut6[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot7(eventsMC[i]))histo_set_cut7[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot8(eventsMC[i]))histo_set_cut8[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot9(eventsMC[i]))histo_set_cut9[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot10(eventsMC[i]))histo_set_cut10[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot11(eventsMC[i]))histo_set_cut11[j][k+nHists]->Fill(val,wErr[k]);
                if(filter_plot12(eventsMC[i]))histo_set_cut12[j][k+nHists]->Fill(val,wErr[k]);
                
                if(filter_plot1(eventsMC[i]) && filter_plot6(eventsMC[i]))histo_set_cutCombo1[j][k+nHists]->Fill(val,wErr[k]);   
                if(filter_plot2(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo2[j][k+nHists]->Fill(val,wErr[k]);   
                if(filter_plot3(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo3[j][k+nHists]->Fill(val,wErr[k]);   
                if(filter_plot4(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo4[j][k+nHists]->Fill(val,wErr[k]);   
            }
            
            if(addSysErrBand)for(unsigned int k = 0; k < nPermErrorBands; k++){
                histo_set[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot1(eventsMC[i]))histo_set_cut1[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot2(eventsMC[i]))histo_set_cut2[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot3(eventsMC[i]))histo_set_cut3[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot4(eventsMC[i]))histo_set_cut4[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot5(eventsMC[i]))histo_set_cut5[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot6(eventsMC[i]))histo_set_cut6[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot7(eventsMC[i]))histo_set_cut7[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot8(eventsMC[i]))histo_set_cut8[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot9(eventsMC[i]))histo_set_cut9[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot10(eventsMC[i]))histo_set_cut10[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot11(eventsMC[i]))histo_set_cut11[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                if(filter_plot12(eventsMC[i]))histo_set_cut12[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);
                
                if(filter_plot1(eventsMC[i]) && filter_plot6(eventsMC[i]))histo_set_cutCombo1[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);   
                if(filter_plot2(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo2[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);   
                if(filter_plot3(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo3[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);   
                if(filter_plot4(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo4[j][k+nHists+nPermErrorBands]->Fill(val,wErrTot[k]);   
            }
            
            for(unsigned int k = 0; k < nAltModels; k++){
                histo_set[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot1(eventsMC[i]))histo_set_cut1[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot2(eventsMC[i]))histo_set_cut2[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot3(eventsMC[i]))histo_set_cut3[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot4(eventsMC[i]))histo_set_cut4[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot5(eventsMC[i]))histo_set_cut5[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot6(eventsMC[i]))histo_set_cut6[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot7(eventsMC[i]))histo_set_cut7[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot8(eventsMC[i]))histo_set_cut8[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot9(eventsMC[i]))histo_set_cut9[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot10(eventsMC[i]))histo_set_cut10[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot11(eventsMC[i]))histo_set_cut11[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                if(filter_plot12(eventsMC[i]))histo_set_cut12[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);
                
                if(filter_plot1(eventsMC[i]) && filter_plot6(eventsMC[i]))histo_set_cutCombo1[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);   
                if(filter_plot2(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo2[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);   
                if(filter_plot3(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo3[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);   
                if(filter_plot4(eventsMC[i]) && filter_plot5(eventsMC[i]))histo_set_cutCombo4[j][k+nHists+nPermErrorBands*(1+addSysErrBand)]->Fill(val,wAlt[k]);   
            }
            
        }
    }
    
    //Plot
    cout << "Create plots ..." << endl;
    TCanvas* c = new TCanvas();      
    TLegend leg(0.,0.,1.0,1,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(132);
    leg.SetTextColor(1);
    leg.SetTextSize(0.065);
    leg.SetTextAlign(12);
    for(unsigned int k=2; k<nHists;k++)leg.AddEntry(histo_set[0][k],legend[k].c_str(),"l");
    
    //Get chi2 for legend
    string resultsFileName = NamedParameter<string>("ResultsFile","result.root");
    TFile* results_file = TFile::Open((outDir+"/"+resultsFileName).c_str(),"OPEN");
    results_file->cd();
    auto result_tree = (TTree*) results_file->Get("Result");
    
    double nll,chi2,sumFractions;
    int status, nPar,nAmps, seed;
    result_tree->SetBranchAddress("nll",&nll);
    result_tree->SetBranchAddress("chi2",&chi2);
    result_tree->SetBranchAddress("status",&status);
    result_tree->SetBranchAddress("nPar",&nPar);
    result_tree->SetBranchAddress("nAmps",&nAmps);
    result_tree->SetBranchAddress("seed",&seed);
    result_tree->GetEntry(0);
    stringstream ss_chi;
    ss_chi << setprecision(3) << chi2;
    leg.AddEntry((TObject *)0,("#chi^{2}/ndf = " + ss_chi.str()).c_str() ,"");
    //leg.AddEntry((TObject *)0,"" ,"");
    leg.AddEntry((TObject *)0,("nPar = " + to_string(nPar)).c_str() ,"");
    leg.AddEntry((TObject *)0,("nAmps = " + to_string(nAmps)).c_str() ,"");
    leg.SetNColumns(2);
    
    //    for(int j=0;j<dims.size();j++){ 
    //        plotData(histo_set[j]);
    //        c->Print(("data_"+labels[j]+".eps").c_str());
    //    }
    
    // Apply smoothing?
    unsigned int smooth   = NamedParameter<unsigned int>("smoothFitProj",0);  
    for (unsigned int i=0 ; i < histo_set[0].size(); i++){
                for(int j=0;j<dims.size();j++){
                    if(i>0 && smooth> 0){
                        histo_set[j][i]->Smooth(smooth);
                        histo_set_cut1[j][i]->Smooth(smooth);
                        histo_set_cut2[j][i]->Smooth(smooth);
                        histo_set_cut3[j][i]->Smooth(smooth);
                        histo_set_cut4[j][i]->Smooth(smooth);
                        histo_set_cut5[j][i]->Smooth(smooth);
                        histo_set_cut6[j][i]->Smooth(smooth);
                        histo_set_cut7[j][i]->Smooth(smooth);
                        histo_set_cut8[j][i]->Smooth(smooth);
                        histo_set_cut8[j][i]->Smooth(smooth);
                        histo_set_cut10[j][i]->Smooth(smooth);
                        histo_set_cut11[j][i]->Smooth(smooth);
                        histo_set_cut12[j][i]->Smooth(smooth);
                        histo_set_cutCombo1[j][i]->Smooth(smooth);
                        histo_set_cutCombo2[j][i]->Smooth(smooth);
                        histo_set_cutCombo3[j][i]->Smooth(smooth);
                        histo_set_cutCombo4[j][i]->Smooth(smooth);
                    }
                    histo_set[j][i]->Rebin(2);
                    histo_set_cut1[j][i]->Rebin(2);
                    histo_set_cut2[j][i]->Rebin(2);
                    histo_set_cut3[j][i]->Rebin(2);
                    histo_set_cut4[j][i]->Rebin(2);
                    histo_set_cut5[j][i]->Rebin(2);
                    histo_set_cut6[j][i]->Rebin(2);
                    histo_set_cut7[j][i]->Rebin(2);
                    histo_set_cut8[j][i]->Rebin(2);
                    histo_set_cut9[j][i]->Rebin(2);
                    histo_set_cut10[j][i]->Rebin(2);
                    histo_set_cut11[j][i]->Rebin(2);
                    histo_set_cut12[j][i]->Rebin(2);
                    histo_set_cutCombo1[j][i]->Rebin(2);
                    histo_set_cutCombo2[j][i]->Rebin(2);
                    histo_set_cutCombo3[j][i]->Rebin(2);
                    histo_set_cutCombo4[j][i]->Rebin(2);
                }
    }
      
    for(unsigned int j=0;j<dims.size();j++){ 
        plotHistos(histo_set[j], true, 1, true);
        c->Print((outDir+"/"+labels[j]+".pdf").c_str());
        c->Print((outDir+"/"+labels[j]+".png").c_str());
        if(j==0)c->Print((outDir+"/"+labels[j]+".C").c_str());
        if(j==0)c->Print((outDir+"/"+labels[j]+".root").c_str());
    }
    
    cout << "f_bkg = " << histo_set[0][histo_set[0].size()-1]->Integral()/histo_set[0][1]->Integral() << endl;
    cout << "f_bkg = " << histo_set[1][histo_set[1].size()-1]->Integral()/histo_set[1][1]->Integral() << endl;

    c->Clear();
    c->Divide(3,2);
    c->SetCanvasSize(1000, 500);
    for(unsigned int j=0;j<6;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],false,0,false);
    }
    c->Print((outDir+"/"+"plots.pdf").c_str());
    
    c->Clear();
    c->Divide(3,2);
    for(unsigned int j=0;j<5;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1,false);
    }
    c->cd(6);
    leg.Draw();
    //c->Print((outDir+"/"+"amp_plots.C").c_str());
    c->Print((outDir+"/"+"amp_plots.pdf").c_str());
    
    c->Clear();
    c->Divide(4,2);
    c->SetCanvasSize(1400, 600);
    for(unsigned int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1,false);
    }
    c->Print((outDir+"/"+"amp_plots2.pdf").c_str());
    
    c->Clear();
    c->Divide(4,2);
    for(unsigned int j=0;j<8;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1,false);
        gPad->SetLogy(1);
    }
    c->Print((outDir+"/"+"amp_plots2_log.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1,false);
    }
    c->Print((outDir+"/"+"amp_plots3.pdf").c_str());
    
        
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut1.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut2.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut3.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut4.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut5.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut6.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut7[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut7.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut8[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut8.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut9[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut9.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut10[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut10.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut11[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plots3_cut11.pdf").c_str());
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut12[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plot3_cut12.pdf").c_str());
    
    
    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cutCombo1[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plot3_cutCombo1.pdf").c_str());

    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cutCombo2[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plot3_cutCombo2.pdf").c_str());

    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cutCombo3[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plot3_cutCombo3.pdf").c_str());

    c->Clear();
    c->Divide(5,2);
    for(unsigned int j=0;j<10;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cutCombo4[j],true,1,true);
    }
    c->Print((outDir+"/"+"amp_plot3_cutCombo4.pdf").c_str());
    
    
    c->SetCanvasSize(500, 500);
    c->Clear();
    leg.Draw();
    c->Print((outDir+"/"+"leg.pdf").c_str());
    
    leg.SetNColumns(2);
    c->Clear();
    leg.Draw();
    c->Print((outDir+"/"+"leg2.pdf").c_str());
    
}


void makePlots3body(){
    
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
    //EventType B0 D- K0S0 pi+ 
    vector<unsigned int> m12{1,2};
    vector<unsigned int> m02{0,2};
    vector<unsigned int> m01{0,1};
    
    vector<vector<unsigned int>> dims{m12,m02,m01,{1},{2},{3}};
    vector<string> labels{"m_Kspi","m_Dpi","m_DKs","cos_Kspi","cos_Dpi","cos_DKs"};
    vector<string> titles{"m(K^{0}_{s}#pi^{#plus}) [GeV]","m(D^{#minus}#pi^{#plus}) [GeV]","m(D^{#minus}K^{0}_{s}) [GeV]","cos(#theta_{K^{0}_{s}#pi^{#plus}})","cos(#theta_{D^{#minus}#pi^{#plus}})","cos(#theta_{D^{#minus}K^{0}_{s}})"
    };
    
    //Limits
    vector<double> lim12{0.5,3.5};
    vector<double> lim02{1.9,5};
    vector<double> lim01{2.1,5.25};
    vector<vector<double>> limits{lim12,lim02,lim01,{-1,1},{-1,1},{-1,1}};
    
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
        for(int j=0;j<dims.size();j++){
            // EventType B0 D- K0S0 pi+ 
            double val = 0;
            if(dims[j].size()>1) val = sqrt(evt.s(dims[j]));
            else if(dims[j][0] == 1) val = helicityCos(evt,1,0,{1,2}); 
            else if(dims[j][0] == 2) val = helicityCos(evt,2,1,{0,2}); 
            else if(dims[j][0] == 3) val = helicityCos(evt,1,2,{0,1}); 

            histo_set[j][0]->Fill(val,evt.weight());
            if(filter_plot1(evt)) histo_set_cut1[j][0]->Fill(val,evt.weight());
            if(filter_plot2(evt))histo_set_cut2[j][0]->Fill(val,evt.weight());
            if(filter_plot3(evt))histo_set_cut3[j][0]->Fill(val,evt.weight());
            if(filter_plot4(evt))histo_set_cut4[j][0]->Fill(val,evt.weight());
            if(filter_plot5(evt))histo_set_cut5[j][0]->Fill(val,evt.weight());
            if(filter_plot6(evt))histo_set_cut6[j][0]->Fill(val,evt.weight());        

//            int n = 0;
//            for(int i=0;i<dims.size();i++)for(int j=i+1;j<dims.size();j++){ histo2D_set[n][0]->Fill(sqrt(evt.s(dims[i])),sqrt(evt.s(dims[j])),evt.weight()); n++;}
            }
    }
    
    //Fill fit projections
    for(int i=0; i< eventsMC.size(); i++ )
    {
        Event evt(eventsMC[i]);
        weight_tree->GetEntry(i);
        
        for(int j=0;j<dims.size();j++){
            // EventType B0 D- K0S0 pi+ 
            double val = 0;
            if(dims[j].size()>1) val = sqrt(evt.s(dims[j]));
            else if(dims[j][0] == 1) val = helicityCos(evt,1,0,{1,2}); 
            else if(dims[j][0] == 2) val = helicityCos(evt,2,1,{0,2}); 
            else if(dims[j][0] == 3) val = helicityCos(evt,1,2,{0,1}); 

            for(int k=1; k<weights.size();k++) histo_set[j][k]->Fill(val,w[k]);
            if(filter_plot1(evt))for(int k=1; k<weights.size();k++) histo_set_cut1[j][k]->Fill(val,w[k]);
            if(filter_plot2(evt))for(int k=1; k<weights.size();k++) histo_set_cut2[j][k]->Fill(val,w[k]);
            if(filter_plot3(evt))for(int k=1; k<weights.size();k++) histo_set_cut3[j][k]->Fill(val,w[k]);
            if(filter_plot4(evt))for(int k=1; k<weights.size();k++) histo_set_cut4[j][k]->Fill(val,w[k]);
            if(filter_plot5(evt))for(int k=1; k<weights.size();k++) histo_set_cut5[j][k]->Fill(val,w[k]);
            if(filter_plot6(evt))for(int k=1; k<weights.size();k++) histo_set_cut6[j][k]->Fill(val,w[k]);

    //        int n = 0;
    //        for(int i=0;i<dims.size();i++)for(int j=i+1;j<dims.size();j++){
    //          for(int k=1; k<weights.size();k++){
    //            histo2D_set[n][k]->Fill(sqrt(evt.s(dims[i])),sqrt(evt.s(dims[j])),w[k]);
    //          }
    //          n++;
    //          }

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
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
            c->cd(j+1);
            plotHistos(histo_set[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots2.eps");
        
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
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots2.eps");

    
    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut1[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set_cut1[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots_cut1.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut2[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set_cut2[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots_cut2.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut3[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set_cut3[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots_cut3.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut4[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set_cut4[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots_cut4.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut5[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set_cut5[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots_cut5.eps");

    c->Clear();
    c->Divide(4,2);
    for(int j=0;j<4;j++){ 
        c->cd(j+1);
        plotHistos(histo_set_cut6[j],true,1);
    }
    for(int j=0;j<3;j++){ 
        c->cd(j+5);
        plotHistos(histo_set_cut6[j],true,1);
        gPad->SetLogy(1);
    }
    c->cd(8);
    gPad->SetLogy(0);
    leg.Draw();    
    c->Print("amp_plots_cut6.eps");

    
    c->Clear();
    leg.Draw();
    c->Print("leg.eps");

//    c->Clear();
//    int n = 0;
//    for(int i=0;i<dims.size();i++)for(int j=i+1;j<dims.size();j++){ 
//      histo2D_set[n][0]->Draw("colz");
//      c->Print(("dalitz_"+to_string(i)+"_"+to_string(j)+".eps").c_str());
//      histo2D_set[n][0]->Draw();
//      c->Print(("dalitzScatter_"+to_string(i)+"_"+to_string(j)+".eps").c_str());
//
//      histo2D_set[n][1]->Draw("colz");
//      c->Print(("dalitzFit_"+to_string(i)+"_"+to_string(j)+".eps").c_str());
//      histo2D_set[n][1]->Draw();
//      c->Print(("dalitzFitScatter_"+to_string(i)+"_"+to_string(j)+".eps").c_str());
//
//      for(int k=2; k<weights.size();k++){
//        histo2D_set[n][k]->Draw("colz");
//        c->Print(("dalitzFit_"+to_string(i)+"_"+to_string(j)+"_amp"+to_string(k)+".eps").c_str());
//        histo2D_set[n][k]->Draw();
//        c->Print(("dalitzFitScatter_"+to_string(i)+"_"+to_string(j)+"_amp"+to_string(k)+".eps").c_str());
//      }
//      n++; 
//    }
}

int main( int argc, char* argv[] ){

  time_t startTime = time(0);
  OptionsParser::setArgs( argc, argv );

  gStyle->SetOptStat(0);
  LHCbStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  auto pNames = NamedParameter<std::string>("EventType" , "").getVector(); 
  std::string outDir = NamedParameter<std::string>("outDir", ".");
  
  if(pNames.size()==4)makePlots3body();
  if(pNames.size()==6)makePlotsMuMu();
  else makePlots();
  
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;  
    
  return 0;
}
