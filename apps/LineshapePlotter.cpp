//// STL / C ////
#include <iostream>
#include <fstream>
#include <math.h>
#include <dlfcn.h>
#include <complex>

//// ROOT ////
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMultiGraph.h"
#include "TFile.h"

//// AMPGEN ////
#include "AmpGen/Lineshapes.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/Particle.h"
#include "AmpGen/EventType.h"
#include "AmpGen/EventList.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/ThreeBodyCalculators.h"
#include "AmpGen/LHCbStyle.h"

using namespace AmpGen;
using namespace std;

void flexible_multigraph( const std::vector<TGraph*> graphs, const std::string& title="" ){

 // TCanvas* c1 = new TCanvas("g1","",500,400);
  TMultiGraph* mg = new TMultiGraph();
  graphs[0]->SetLineColor(kBlue);
  graphs[1]->SetLineColor(kBlack);
  graphs[1]->SetMarkerStyle(5);
  graphs[1]->SetLineWidth(2);
  graphs[1]->SetMarkerSize(1.4); 
  graphs[1]->SetMarkerColor(kBlack);
  mg->Add( graphs[0],"C");
  mg->Add( graphs[1],"P");
  //  for( auto& g : graphs ) mg->Add(g,"CP");
  
  mg->Draw("ACP");
  mg->GetYaxis()->SetTitle( graphs[0]->GetYaxis()->GetTitle() );
  mg->GetXaxis()->SetTitle( graphs[0]->GetXaxis()->GetTitle() );
  mg->Draw("ACP"); 
 // c1->SaveAs( title.c_str() );
 // delete mg;
 // delete c1;
}

int test(int argc , char* argv[] ){
  //setupStyle();
  EventType evt( NamedParameter<std::string>("EventType").getVector() );
  bool debug                = NamedParameter<unsigned int>("debug",0);
  std::string     lineshape = NamedParameter<std::string>("Lineshape");
  std::string ref_lineshape = NamedParameter<std::string>("RefLineshape","");
  
  unsigned int orbital = NamedParameter<unsigned int>("L",0);
  std::vector<DebugSymbol> db   ;
  Expression BW = Lineshape::Factory::get(lineshape, Parameter("s"),
      evt.mass(0)*evt.mass(0) , evt.mass(1)*evt.mass(1) , evt.mother() , orbital, &db  ) ;
  MinuitParameterSet mps; mps.loadFromStream();
  
  auto props = ParticlePropertiesList::getMe()->get( evt.mother() );
    /*
  if( ref_lineshape != "" ){
    ReferenceAmplitude* refAmp =  Factory<ReferenceAmplitude>::get(ref_lineshape);
    refAmp->set(props->mass(), props->width(), props->radius(), evt.mass(0)*evt.mass(0), evt.mass(1)*evt.mass(1) , orbital , "");
  if( debug ){
    INFO("ENTERING DEBUGGER");
    double s = 0.6*1000*1000;
    //bwc.debug(&s);
   // bwc.print();
    refAmp->debug(s );
  }
  }
     */
  auto limits = NamedParameter<double> ("Limits").getVector();
  double min = limits[0];
  double max = limits[1];
  double step = limits[2];

  INFO("Plotting with limits: " << min << " max = " << max << " step = " << step );
  /*
  CompiledExpression<std::complex<double>,  const double*, const double*  > bwc( BW, "Lineshape", {{"s",0}}, db, &mps );

  double tmp = (1000.*1000.);

  TGraph* amp = new TGraph();
  TGraph* phase = new TGraph();
  TGraph* amp2 = new TGraph();
  TGraph* phase2 = new TGraph();

  TGraph* real = new TGraph();
  TGraph* imag = new TGraph();
  TGraph* real2 = new TGraph();
  TGraph* imag2 = new TGraph();
  
  unsigned int j = 0 ;

  for( double s = min; s < max; s+=step/10. ){
    std::complex<double> value2 = bwc( &s );
    //std::complex<double> lau_value = (*refAmp)( s );
    double m = sqrt(s);
    amp->SetPoint( amp->GetN(),     s, sqrt( std::norm( value2 ) ) );
    phase->SetPoint( phase->GetN(), s, 180 * std::arg( value2 ) / M_PI );
    real->SetPoint( real->GetN(),   s, std::real( value2 ) );
    imag->SetPoint( imag->GetN(),   s, std::imag( value2 ) );
    //INFO( value2 << " , " << lau_value );
  //  if( ++j == 10 ){
  //    imag2->SetPoint( imag2->GetN(),     m, std::imag(lau_value ) );
  //    real2->SetPoint( real2->GetN(),     m, std::real(lau_value ) );
  //    amp2->SetPoint( amp2->GetN(),       m, std::norm( lau_value ) );
  //    phase2->SetPoint( phase2->GetN(),   m, std::arg( lau_value ) );
  //    j=0;
  //  }
  }
//  TCanvas* c1 = new TCanvas("c1","",1000,800);
//  c1->Divide(2,2);
  std::string xAxisTitle = NamedParameter<std::string>("axisTitle");
  std::string cnvName    = NamedParameter<std::string>("cnvName");
  std::string yAxisTitle = NamedParameter<std::string>("yAxisName","\\mathcal{A}(s)"); 
  real->GetYaxis()->SetTitle("Re(A)");
  imag->GetYaxis()->SetTitle("Im(A)");
  phase->GetYaxis()->SetTitle(("\\mathrm{arg}\\left("+yAxisTitle+"\\right) \\left[^\\mathrm{o}\\right]").c_str() );
  amp->GetYaxis()->SetTitle(  ("\\left|" + yAxisTitle+ "\\right|").c_str() );
  
  real->GetXaxis()->SetTitle( xAxisTitle.c_str());
  imag->GetXaxis()->SetTitle( xAxisTitle.c_str());
  phase->GetXaxis()->SetTitle(xAxisTitle.c_str());
  amp->GetXaxis()->SetTitle(  xAxisTitle.c_str());
  
  amp->SetLineWidth(2);
  phase->SetLineWidth(2);
  amp->SetMarkerSize(0);
  phase->SetMarkerSize(0);
  TCanvas* c1 = new TCanvas("c1","",500,400);
  amp->Draw();
  LatexViewer().view( c1 ,"amp.pdf");
  TCanvas* c2 = new TCanvas("c2","",500,400);
  phase->Draw();
  LatexViewer().view(c2,"phase.pdf");
  //  c1->cd(1); flexible_multigraph( {real,real2} );
//  c1->cd(2); flexible_multigraph( {imag,imag2} );
//  c1->cd(3); flexible_multigraph( {amp,amp2} );
//  c1->cd(4); flexible_multigraph( {phase,phase2} );
//  c1->SaveAs(cnvName.c_str());
  */
}

double getX( const TGraph* g, const size_t& i ){
  double x,y;
  g->GetPoint(i,x,y);
  return x; 
}

double getY( const TGraph* g, const size_t& i ){
  double x,y;
  g->GetPoint(i,x,y);
  return y; 
}

void plotRunningWidths(){
    
    MinuitParameterSet mps; 
    mps.loadFromStream();
    
    std::string outDir = NamedParameter<std::string>("outDir", ".");

    std::vector<std::string> threeBodiesToIntegrate = NamedParameter<std::string>( "ThreeBodiesToIntegrate" ).getVector();
    for ( auto& head : threeBodiesToIntegrate ){
    
        double mass0      = mps( head+"_mass");  
        double width0      = mps( head+"_width");    

        auto spline_params = NamedParameter<double>( head + "::Spline").getVector();
        size_t nBins;
        double min, max;
        if( spline_params.size() == 3 ){
            nBins = size_t( spline_params[0] );
            min   =         spline_params[1] ; 
            max   =         spline_params[2];
        }
        else {
            nBins = NamedParameter<double>( head + "::Spline::N"  , 0. );
            min   = NamedParameter<double>( head + "::Spline::Min", 0. );
            max   = NamedParameter<double>( head + "::Spline::Max", 0. );
        }
        
        cout << endl << "Plot running width for " << head << endl; 
        cout << left << head << "::Spline" << "  " << nBins << " " << min << " " << max << endl; 

        double step   = ( max - min ) / double(nBins-1);
        
        TGraph* width_m = new TGraph(nBins);
        TGraph* width_m_norm = new TGraph(nBins);

        for ( size_t c = 0; c < nBins; ++c ) {
            double s = min + double(c) * step;
            double val = mps( head+"::Spline::Gamma::"+to_string(c));
            width_m->SetPoint(c, sqrt(s), val );
        }

        double norm = width_m->Eval(mass0);
        for ( size_t c = 0; c < nBins; ++c ) {
            double s = min + double(c) * step;
            double val = mps( head+"::Spline::Gamma::"+to_string(c));
            width_m_norm->SetPoint(c, sqrt(s), val/norm * width0 );
        }
        
                
        TCanvas* c = new TCanvas();

        TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                            0.87 - gStyle->GetPadTopMargin(),
                                            gStyle->GetPadLeftMargin() + 0.20,
                                            0.95 - gStyle->GetPadTopMargin(),
                                            "BRNDC");
        lhcbName->AddText("LHCb");
        lhcbName->SetFillColor(0);
        lhcbName->SetTextAlign(12);
        lhcbName->SetBorderSize(0);
        lhcbName->SetTextSize(0.08);
        lhcbName->SetTextFont(132);

        TPaveText *text= new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                       0.77 - gStyle->GetPadTopMargin(),
                                       gStyle->GetPadLeftMargin() + 0.20,
                                       0.85 - gStyle->GetPadTopMargin(),
                                       "BRNDC");
        auto fitResult = new FitResult();
        text->AddText(fitResult->latexName(head).c_str());
        text->SetLineColor(kWhite);
        text->SetFillColor(kWhite);
        text->SetShadowColor(0);
        text->SetTextAlign(12);
        text->SetTextSize(0.08);
        text->SetTextFont(132);
        text->SetTextColor(kBlack);     

        width_m_norm->SetLineColor(kBlue);
        width_m_norm->SetTitle("; #sqrt{#it{s}} [GeV]  ; #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
        width_m_norm->Draw("A*C");
        text->Draw();
        c->Print( ( outDir + "/" + head + "_runningWidth.pdf").c_str());

    }    
        
}


void calculateRunningWidths(){
    
    MinuitParameterSet mps; 
    mps.loadFromStream();
    
    std::string outDir = NamedParameter<std::string>("outDir", ".");
    std::string outName = NamedParameter<std::string>("outName", "splineKnots.txt");
    ofstream resultsFile;
    resultsFile.open((outDir + "/" + outName).c_str(),std::ofstream::trunc);

    std::vector<std::string> threeBodiesToIntegrate = NamedParameter<std::string>( "ThreeBodiesToIntegrate" ).getVector();
    for ( auto& head : threeBodiesToIntegrate ){
    
        double mass0      = mps( head+"_mass");  
        double width0      = mps( head+"_width");    

        auto spline_params = NamedParameter<double>( head + "::Spline").getVector();
        size_t nBins;
        double min, max;
        if( spline_params.size() == 3 ){
            nBins = size_t( spline_params[0] );
            min   =         spline_params[1] ; 
            max   =         spline_params[2];
        }
        else {
            nBins = NamedParameter<double>( head + "::Spline::N"  , 0. );
            min   = NamedParameter<double>( head + "::Spline::Min", 0. );
            max   = NamedParameter<double>( head + "::Spline::Max", 0. );
        }
        
        cout << endl << "Running width for " << head << " from AmpGen "<< endl; 
        resultsFile << "# Running width for " << head << " from AmpGen " << endl; 
        cout << left << head << "::Spline" << "  " << nBins << " " << min << " " << max << endl; 
        resultsFile << left << head << "::Spline" << "  " << nBins << " " << min << " " << max << endl;  

        ThreeBodyCalculator calc ( head ,mps, nBins , min, max );
        auto width = calc.widthGraph( mass0 * mass0);
        
        auto width_m = new TGraph(*width);
        auto width_m_norm = new TGraph(*width);

        double norm = calc.getWidth(mass0*mass0);
        
        for (int i=0; i<width_m->GetN(); i++) {
            width_m->SetPoint(i, sqrt(getX(width,i)), getY(width,i)/norm * width0 );
            width_m_norm->SetPoint(i, sqrt(getX(width,i)), getY(width,i)/norm  );            
        }
        
        TCanvas* c = new TCanvas();

        TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                            0.87 - gStyle->GetPadTopMargin(),
                                            gStyle->GetPadLeftMargin() + 0.20,
                                            0.95 - gStyle->GetPadTopMargin(),
                                            "BRNDC");
        lhcbName->AddText("LHCb");
        lhcbName->SetFillColor(0);
        lhcbName->SetTextAlign(12);
        lhcbName->SetBorderSize(0);
        lhcbName->SetTextSize(0.08);
        lhcbName->SetTextFont(132);

        TPaveText *text= new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                       0.77 - gStyle->GetPadTopMargin(),
                                       gStyle->GetPadLeftMargin() + 0.20,
                                       0.85 - gStyle->GetPadTopMargin(),
                                       "BRNDC");
        auto fitResult = new FitResult();
        text->AddText(fitResult->latexName(head).c_str());
        text->SetLineColor(kWhite);
        text->SetFillColor(kWhite);
        text->SetShadowColor(0);
        text->SetTextAlign(12);
        text->SetTextSize(0.08);
        text->SetTextFont(132);
        text->SetTextColor(kBlack);     

        width_m->SetLineColor(kBlue);
        width_m->SetTitle("; #sqrt{#it{s}} [GeV]  ; #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
        width_m->Draw("A*C");
        //lhcbName->Draw();
        text->Draw();
        c->Print( ( outDir + "/" + head + "_runningWidth.pdf").c_str());

        for( size_t i = 0 ; i < nBins ; ++i ){
          cout << left << head << "::Spline::Gamma::" << setw(20) << i << "  2   " << setw(14) << getY(width,i) << " 0 " << endl; 
          resultsFile << left << head << "::Spline::Gamma::" << setw(20) << i << "  2   " << setw(14) << getY(width,i) << " 0 " << endl; 
        }
        resultsFile << endl;

    }    
        
}

void prepareRunningWidthFromFiles(){
        
    MinuitParameterSet mps; 
    mps.loadFromStream();
    
    vector<string> threeBodiesToIntegrate = NamedParameter<string>( "ThreeBodiesToIntegrate" ).getVector();
    string inDir = NamedParameter<string>("RunningWidthFilesDir", ".");
    vector<string> files = NamedParameter<string>( "RunningWidthFiles" ).getVector();
    std::string outDir = NamedParameter<std::string>("outDir", ".");
    std::string outName = NamedParameter<std::string>("outName", "splineKnots.txt");

    ofstream resultsFile;
    resultsFile.open((outDir + "/" + outName).c_str(),std::ofstream::trunc);

    int counter = 0;
    for ( auto& head : threeBodiesToIntegrate ){
        
        cout << endl << "Running width for " << head << " from file " << files[counter] << endl; 
        
        resultsFile << "# Running width for " << head << " from file " << files[counter] << endl; 

        TFile* f=TFile::Open((inDir+files[counter]).c_str());
        TH1D* h=dynamic_cast<TH1D*>(f->Get("RunningWidth"));
        counter++;
        
        auto spline_params = NamedParameter<double>( head + "::Spline").getVector();
        size_t nBins;
        double min, max;
        if( spline_params.size() == 3 ){
            nBins = size_t( spline_params[0] );
            min   =         spline_params[1] ; 
            max   =         spline_params[2];
        }
        else {
            nBins = NamedParameter<double>( head + "::Spline::N"  , 0. );
            min   = NamedParameter<double>( head + "::Spline::Min", 0. );
            max   = NamedParameter<double>( head + "::Spline::Max", 0. );
        }

        double mass0      = mps.find( head+"_mass") != 0 ? mps( head+"_mass") : sqrt( min + (max-min)/2. );  
        double width0      = mps.find( head+"_width") != 0 ? mps( head+"_width") : 0.1;

        cout << "mass0 = " << mass0 << endl; 
        cout << "width0 = " << width0 << endl; 
        cout << left << head << "::Spline" << "  " << nBins << " " << min << " " << max << endl; 
        resultsFile << left << head << "::Spline" << "  " << nBins << " " << min << " " << max << endl;  
        
        double step   = ( max - min ) / double(nBins-1);
        
        TGraph* width_m = new TGraph(nBins);
        TGraph* width_m_norm = new TGraph(nBins);
        double norm = h->Interpolate(mass0*mass0);

        for ( size_t c = 0; c < nBins; ++c ) {
            double s = min + double(c) * step;
            double val = h->Interpolate(s);

            width_m->SetPoint(c, sqrt(s), val/norm * width0 );
            width_m_norm->SetPoint(c, sqrt(s), val/norm  );            

            cout << left << head << "::Spline::Gamma::" << setw(20) << c << "  2   " << setw(14) << val << " 0 " << endl; 
            resultsFile << left << head << "::Spline::Gamma::" << setw(20) << c << "  2   " << setw(14) << val << " 0 " << endl; 
        }
        resultsFile << endl;
        
        TCanvas* c = new TCanvas();
        
        TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                            0.87 - gStyle->GetPadTopMargin(),
                                            gStyle->GetPadLeftMargin() + 0.20,
                                            0.95 - gStyle->GetPadTopMargin(),
                                            "BRNDC");
        lhcbName->AddText("LHCb");
        lhcbName->SetFillColor(0);
        lhcbName->SetTextAlign(12);
        lhcbName->SetBorderSize(0);
        lhcbName->SetTextSize(0.08);
        lhcbName->SetTextFont(132);

        TPaveText *text= new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                       0.77 - gStyle->GetPadTopMargin(),
                                       gStyle->GetPadLeftMargin() + 0.20,
                                       0.85 - gStyle->GetPadTopMargin(),
                                       "BRNDC");
        auto fitResult = new FitResult();
        text->AddText(fitResult->latexName(head).c_str());
        text->SetLineColor(kWhite);
        text->SetFillColor(kWhite);
        text->SetShadowColor(0);
        text->SetTextAlign(12);
        text->SetTextSize(0.08);
        text->SetTextFont(132);
        text->SetTextColor(kBlack);            

        width_m->SetLineColor(kBlue);
        width_m->SetTitle("; #sqrt{#it{s}} [GeV]  ; #sqrt{#it{s}} / #it{m_{0} #it{#Gamma(s)}} [GeV]");
        width_m->Draw("A*C");
        //lhcbName->Draw();
        text->Draw();

        TString n(head);
        n.ReplaceAll("*","");
        n.ReplaceAll("+","");
        n.ReplaceAll("-","");
        n.ReplaceAll(")0","");
        n.ReplaceAll("(","_");
        n.ReplaceAll(")","");

        c->Print( ( outDir + "/" + (string) n + "_runningWidth.pdf").c_str());        
    }
    
}

int main(int argc , char* argv[] ){

      OptionsParser::setArgs(argc,argv);
      LHCbStyle();
      TH1::SetDefaultSumw2();
      gStyle->SetOptStat(0);
      gStyle->SetTitleOffset(0.9,"X");
      gStyle->SetTitleOffset(0.8,"Y");
      gStyle->SetTitleSize(0.08,"x");
      gStyle->SetTitleSize(0.08,"y");
      gStyle->SetLabelOffset(0.005,"X");
      gStyle->SetLabelOffset(0.005,"Y");
      gStyle->SetLabelSize(0.065,"x");
      gStyle->SetLabelSize(0.065,"y");
    
      auto plotRw = NamedParameter<bool>("plotRunningWidths", 0);
      auto calculateRw = NamedParameter<bool>("calculateRunningWidths", 0);
      auto prepareRw = NamedParameter<bool>("prepareRunningWidthFromFiles", 0);
      auto plotSpline = NamedParameter<bool>("plotSpline", 0);

      if(plotRw)plotRunningWidths();
      if(calculateRw)calculateRunningWidths();
      if(prepareRw)prepareRunningWidthFromFiles();
    
      if(plotSpline){
        string logFile = NamedParameter<std::string>("LogFile", "log.txt", "Name of the output log file");
        string head = NamedParameter<std::string>("Head","X(A)0");

        FitResult fr(logFile);
        fr.plotSpline(head);    
    }
}
