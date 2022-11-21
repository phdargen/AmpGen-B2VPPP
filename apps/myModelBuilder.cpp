#include <iostream>
#include <sys/types.h> 
#include <dirent.h> 
#include <chrono>
#include <thread>
#include <algorithm>
#include <fstream> 
#include <sstream>
#include <ctime>

#include "AmpGen/Utilities.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Factory.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/LHCbStyle.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include <TH1.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

using namespace AmpGen;
using namespace std;

struct model {
    string name;
    int id;
    int gen;
    double nll;
    double chi2;    
} ;

std::vector<std::string> cmd( const std::string& command ){
    std::vector<std::string> output;
    FILE* proc = popen(command.c_str(), "r");
    char buf[4096];
    // Loop over the entries.
    while (!feof(proc) && fgets(buf, sizeof(buf), proc)) {
        output.push_back( buf );
    }
    pclose(proc);
    return output;
}

template <class iterator, class comparator>
std::vector<iterator> max_elements( iterator begin, 
                                   iterator end, 
                                   const unsigned int& N,
                                   comparator compare )
{
    std::vector<iterator> sorted;
    for( auto& it = begin; it != end ; ++it) sorted.push_back(it);
    std::sort( sorted.begin() , sorted.end() , compare );
    std::vector<iterator> partial( sorted.begin(), sorted.begin() + N );
    return partial;
}

std::map<std::string,bool> getRunningJobs(int clusterID){
  std::map<std::string,bool> running_jobs; 
  auto qstat_lines = cmd("condor_q " + to_string(clusterID) + " -format '%d\n' ProcId -run");
  for( auto& line : qstat_lines ) running_jobs[line] = true; 
  return running_jobs ; 
}

std::map<std::string,bool> getFailedJobs(int clusterID){
    std::map<std::string,bool> failed_jobs; 
    auto qstat_lines = cmd("condor_q " +  to_string(clusterID) + " -format '%d\n' ProcId -hold");
    for( auto& line : qstat_lines ) failed_jobs[line] = true; 
    return failed_jobs ; 
}

std::map<std::string,bool> getJobs(int clusterID){
    std::map<std::string,bool> jobs; 
    auto qstat_lines = cmd("condor_q " +  to_string(clusterID) + " -format '%d\n' ProcId");
    for( auto& line : qstat_lines ) jobs[line] = true; 
    return jobs ; 
}

bool readFile(const std::string& fname ){
    std::ifstream checkIsClosed( fname );
    if ( !checkIsClosed.is_open() || checkIsClosed.peek() == std::ifstream::traits_type::eof() ) return false;
    checkIsClosed.close();
    auto lines = vectorFromFile(fname);
    if( *lines.rbegin() != "End Log" ){
        ERROR("File not properly closed: " << *lines.rbegin() );
        return false; 
    }
    return true;
}

void killJobs(int clusterID){
    cmd("condor_hold " +  to_string(clusterID) + " -constraint 'JobStatus!=4' ");    
}

string makeReducedAmpList(MinuitParameterSet* mps, vector<string> ampsToRemove, const std::string& fname){
    
        std::ofstream outlog;
        //outlog << std::setprecision( 4 );
        outlog.open( fname );
        for (size_t i = 0; i < mps->size(); ++i ) {
            auto param = mps->at(i);
            //INFO(param->name());
            if(param->isHidden() || param->name().find( "::Spline" ) != std::string::npos )continue;

            bool found = false;
            for (auto& s : ampsToRemove) if( param->name().find(s) != std::string::npos ) found = true;
            if(found)continue;
            
            if(param->name().find( "_Re" ) != std::string::npos){
                outlog << replaceAll(param->name(),"_Re", "") << "  " << (int)param->flag() << " " 
                << param->mean() << " " << param->stepInit() << " " << param->minInit() << " " << param->maxInit() << "    " ;     
                
                auto param_im = mps->at(i+1);
                outlog << (int)param_im->flag() << " " << param_im->mean() << " " << param->stepInit() << " " << param_im->minInit() << " " << param_im->maxInit() << " " ;               
                i++;            
            }
//            else{
//                outlog << param->name() << "  " << (int)param->flag() << " " 
//                << param->mean() << " " <<  param->stepInit() << " "
//                << param->minInit() << " " << param->maxInit() << "    " ;               
//            }
            outlog << std::endl;
        }
        outlog.close();
        return fname;
}

vector<model> selectModel(string dir = "", MinuitParameterSet* mps = 0, int verbose = 0){
    
    string resultsFile = NamedParameter<string>("ResultsFile","result.root");
    TChain* chain = new TChain("Result");
    chain->Add((replaceAll(dir+"/"+resultsFile,".root","_*.root")).c_str());  
    
    const int n = chain->GetEntries();
    INFO("Analyzing " << n << " result files");
    
    double nll,chi2,sumFractions;
    int status, nPar,nAmps, seed;
    chain->SetBranchAddress("nll",&nll);
    chain->SetBranchAddress("chi2",&chi2);
    chain->SetBranchAddress("status",&status);
    chain->SetBranchAddress("nPar",&nPar);
    chain->SetBranchAddress("nAmps",&nAmps);
    chain->SetBranchAddress("seed",&seed);
    
    chain->GetEntry(0);
    double min_nll = nll; 
    double max_nll = nll; 
    double min_chi2 = chi2;   
    int min_nPar = nPar;  

    int min_i_LL = 0;  
    int min_i_chi2 = 0;  
    int min_i_BIC = 0;  
    
    for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
        
        if(nll>max_nll)max_nll=nll;
        
        if(nll<min_nll){ 
            min_nll=nll;
            min_i_LL=i;
            min_nPar = nPar;
        }
        
        if(chi2<min_chi2){
            min_chi2=chi2;      
            min_i_chi2=i;
        }
    }
    
    TH1D * h_nll = new TH1D("h_nll","h_nll",50,-0.5,min(max_nll-min_nll,3000.5));      
    h_nll->SetLineColor(kBlue);  
    TH1D * h_chi2 = new TH1D("h_chi2","h_chi2",50,0,10);  
    h_chi2->SetLineColor(kBlue);  

    vector<model> models;
    
    for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
        
        h_nll->Fill(nll-min_nll);
        h_chi2->Fill(chi2);
        
        double AIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) ;
        double BIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) * log(30000) ;

        model m = {.name=replaceAll(mps->at(seed*2)->name(),"_Re",""), .id = seed ,.nll=nll, .chi2=chi2};
        models.push_back(m);
        
        if(verbose){
            cout << endl << "model " << seed << endl;
            if(mps != 0)INFO("Added amplitude " << replaceAll(mps->at(seed*2)->name(),"_Re",""));  
            cout << "n2ll = " << nll << " (" << sqrt(nll-min_nll) << " sigma)" << endl;
            cout << "chi2 = " << chi2 << endl;
            cout << "status = " << status << endl;
            cout << "nAmps = " << nAmps << endl;
            cout << "nPar = " << nPar << endl;
            cout << "AIC = " << AIC << " ( " << sqrt( AIC ) << " sigma) " << endl;
            cout << "BIC = " << BIC << " ( " << sqrt( BIC ) << " sigma) " << endl;
        }
    }

    if(verbose){
        TCanvas* c = new TCanvas("c");    
        h_nll->Draw("h");
        c->Print((dir+"/"+"n2ll.pdf").c_str());  
        
        h_chi2->Draw("h");
        c->Print((dir+"/"+"chi2.pdf").c_str());

        chain->GetEntry(min_i_LL);
        cout << "==============================================" << endl;
        cout << endl << "Best model (LL): " << seed << endl;
        if(mps != 0)INFO("Added amplitude " << replaceAll(mps->at(seed*2)->name(),"_Re",""));  
        cout << "n2ll = " << nll << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "status = " << status << endl;
        cout << "nAmps = " << nAmps << endl;
        cout << "nPar = " << nPar << endl;

        chain->GetEntry(min_i_chi2);
        cout << "==============================================" << endl;
        cout << endl << "Best model (chi2): " << seed << endl;
        if(mps != 0)INFO("Added amplitude " << replaceAll(mps->at(seed*2)->name(),"_Re",""));  
        cout << "n2ll = " << nll << endl;
        cout << "chi2 = " << chi2 << endl;
        cout << "status = " << status << endl;
        cout << "nAmps = " << nAmps << endl;
        cout << "nPar = " << nPar << endl;
        delete c;
    }
    
    delete chain;
    delete h_nll;
    delete h_chi2;
    
    return models;
}

void analyzeModelBuilder(int generation, int nFits, string project_dir){
    LHCbStyle();
    
    string resultsFile = NamedParameter<string>("ResultsFile","result.root");
    resultsFile = replaceAll(resultsFile,".root","");
    TChain* chain = new TChain("Result");

    for (int i=0; i< generation; i++ ) {
        for (int j=0; j< nFits; j++ ) { 
            if(std::ifstream((project_dir + "/v" + to_string(i) + "/" +resultsFile+"_"+to_string(j)+".root").c_str()).good())
                chain->Add( (project_dir + "/v" + to_string(i) + "/" +resultsFile+"_"+to_string(j)+".root").c_str());          
        }
    }
    
    INFO("Analyzing " << chain->GetEntries() << " result files");
    
    double nll,chi2,sumFractions;
    int status, nPar,nAmps, seed;
    chain->SetBranchAddress("nll",&nll);
    chain->SetBranchAddress("chi2",&chi2);
    chain->SetBranchAddress("status",&status);
    chain->SetBranchAddress("nPar",&nPar);
    chain->SetBranchAddress("nAmps",&nAmps);
    chain->SetBranchAddress("seed",&seed);
    chain->SetBranchAddress("Sum_Bp",&sumFractions);

    TCanvas* c = new TCanvas();    
    vector<TGraph*> graphs_chi2;
    vector<int> bestInGeneration;
    
    int n = 0;
    for (int i=0; i< generation; i++ ) {
        
        double x_gen[nFits];
        double y_chi2[nFits];
        
        int best = 0;
        double min_chi2 = 999999.;
        
        for (int j=0; j< nFits; j++ ) {        
            if(!std::ifstream((project_dir + "/v" + to_string(i) + "/" +resultsFile+"_"+to_string(j)+".root").c_str()).good())continue;
            chain->GetEntry(n);
            x_gen[j] = (double)i;
            y_chi2[j] = chi2;     
            
            if(chi2 < min_chi2){
                min_chi2 = chi2;
                best = n;
            }
            n++;
        }

        TGraph* g_chi2 = new TGraph(nFits,x_gen,y_chi2); 
        g_chi2->SetMinimum(0.1);
        g_chi2->SetMaximum(3);
        g_chi2->SetMarkerStyle(20);
        g_chi2->SetMarkerSize(1.2);
        g_chi2->SetMarkerColor(kBlue);
        g_chi2->SetLineColor(kBlue);     
        g_chi2->GetXaxis()->SetLimits(-0.5,generation+0.5);
        g_chi2->SetTitle(";Iteration ; #chi^{2}/ndf");
        graphs_chi2.push_back(g_chi2);
        
        bestInGeneration.push_back(best);
    }

    graphs_chi2[0]->Draw("AP");
    for(int i=1; i < graphs_chi2.size(); i++)graphs_chi2[i]->Draw("P");
    c->Print((project_dir+"/"+"chi2_scanAll.pdf").c_str());
    
    double x_best[generation];
    double y_chi2_best[generation];
    double y_sumFractions_best[generation];

    for (int i=0; i< generation; i++ ){
        chain->GetEntry(bestInGeneration[i]);
        
        x_best[i] = (double)i;
        y_chi2_best[i] = chi2;  
        y_sumFractions_best[i] = sumFractions * 100.;  
    }
    
    TGraph* g_chi2 = new TGraph(generation,x_best,y_chi2_best); 
    g_chi2->SetMinimum(0.5);
    g_chi2->SetMaximum(3);
    g_chi2->SetMarkerStyle(20);
    g_chi2->SetMarkerSize(1.2);
    g_chi2->SetMarkerColor(kBlue);
    g_chi2->SetLineColor(kBlue);     
    g_chi2->GetXaxis()->SetLimits(-0.5,generation+0.5);
    g_chi2->SetTitle(";Iteration ; #chi^{2}/ndf");
    g_chi2->Draw("AP");
    c->Print((project_dir+"/"+"chi2_scan.pdf").c_str());    
    
    TGraph* g_sumFractions = new TGraph(generation,x_best,y_sumFractions_best); 
    g_sumFractions->SetMinimum(0);
    g_sumFractions->SetMaximum(200);
    g_sumFractions->SetMarkerStyle(20);
    g_sumFractions->SetMarkerSize(1.2);
    g_sumFractions->SetMarkerColor(kBlue);
    g_sumFractions->SetLineColor(kBlue);     
    g_sumFractions->GetXaxis()->SetLimits(-0.5,generation+0.5);
    g_sumFractions->SetTitle(";Iteration ; Sum of fit fractions [%]");
    g_sumFractions->Draw("AP");
    c->Print((project_dir+"/"+"sumFrac_scan.pdf").c_str());   
    
}

int main(int argc, char* argv[] )
{
  gStyle->SetOptStat(0);
    
  time_t startTime = time(0);
  OptionsParser::setArgs( argc, argv );
  bool submit = NamedParameter<bool>("submit",true);

  string project_dir  = NamedParameter<std::string>("Project","modelSel");
  string in_dir  = NamedParameter<std::string>("in_dir",".");    
  string starterModel  = NamedParameter<std::string>("StarterModel","core.txt");
  string addAmpList = NamedParameter<string>("addAmpList","addAmpList.txt"); // Same as in main_options.txt !
  string modelSelectionFoM  = NamedParameter<std::string>("modelSelectionFoM","chi2");
  double delta_chi2  = NamedParameter<double>("delta_chi2",0.0);
  double delta_nll  = NamedParameter<double>("delta_nll",9);

  int fixedFitsToDo = NamedParameter<int>("fixedFitsToDo", -1);  
  int nWorstAmpsToRemoveFromList = NamedParameter<int>("nWorstAmpsToRemoveFromList", 0);  

  int generation = NamedParameter<int>("generation", 0);  
  int maxGeneration = NamedParameter<int>("maxGeneration", 15);  
  double maxTimePerGen  = NamedParameter<double>("maxTimePerGen",4);

  int status=0;
  status &= system( ("mkdir -p " + project_dir).c_str() );
  if( status != 0 ) ERROR("Building project directories failed");
  
  bool improved = true;
  model bestModel = {.name="start", .id = -1,.gen=generation, .nll=99999, .chi2=99999};
    
  do{  
      time_t startTimeGen = time(0);
      cout << endl <<  "==============================================" << endl;    
      INFO("Model builder:: Starting generation " << generation);
      
      status &= system( ("mkdir -p " + project_dir+"/v"+ to_string(generation)).c_str() );
      if( status != 0 ) ERROR("Building project directories failed");
      cmd("kinit -R");

      MinuitParameterSet* addAmp = new MinuitParameterSet();
      addAmp->loadFromFile(in_dir + "/" + addAmpList);  
      int fitsToDo = fixedFitsToDo < 0 ? addAmp->size()/2 : fixedFitsToDo;
      
      const std::string command = "condor_submit o=" + starterModel + " i=" + in_dir + " d=" + project_dir+"/v"+to_string(generation) 
      + " l=" + addAmpList + " n=" + to_string(fitsToDo) + " submitModelBuilder.sub";  
      
      int clusterID; 
      if(submit){  
          INFO(command);
          auto status2 = cmd( command.c_str() );  
          for( auto& s : status2) {
              vector<string> cluster_str;
              if(s.find("cluster") != std::string::npos){
                  cluster_str = split(s,' ');
                  clusterID =  atoi((cluster_str[cluster_str.size()-1]).c_str());
              }
          }
          auto nJobs = getJobs(clusterID); 
          INFO(nJobs.size() << "/" << fitsToDo << " jobs submitted with clusterID " << clusterID);
      }
 
      unsigned int nCompleted = 0 ;     
      unsigned int nFailed = 0 ;       
      if(submit)do {
          nCompleted = 0 ;     
          for (unsigned i=0; i< fitsToDo; i++) {
              if(readFile(project_dir+"/v"+to_string(generation)+"/log_"+to_string(i)+".txt"))nCompleted++;
          }
          if( nCompleted < fitsToDo ){
             auto runningJobs = getRunningJobs(clusterID); 
             nFailed = getFailedJobs(clusterID).size() ;     

             cout << endl << " Time since start " << (time(0) - startTime)/60.0 << " min." << endl;
             INFO(fitsToDo - runningJobs.size() - nCompleted - nFailed << " jobs idle");
             INFO(runningJobs.size() << " jobs running");
             INFO(nCompleted << " jobs completed");
             INFO(nFailed << " jobs failed");   
              
             if(nFailed>fitsToDo/4){
                 ERROR("Too many failed jobs");
                 //killAllJobs();
                 return 0;
             }
             
             sleep(300);
              
             if((time(0) - startTimeGen)/60./60. > maxTimePerGen){
                 ERROR(fitsToDo - nCompleted - nFailed  << "jobs exceed time limit, will terminate them");
                 killJobs(clusterID); 
              }
             continue;
          }
          else INFO("Jobs finished after " << (time(0) - startTimeGen)/60 << " mins");
      } while( nCompleted + nFailed < fitsToDo );
      

      vector<model> newModels = selectModel(project_dir+"/v"+to_string(generation), addAmp);    
      if(modelSelectionFoM=="chi2")sort( newModels.begin(), newModels.end(), []( auto& m1, auto& m2 ){ return m1.chi2 < m2.chi2 ; } );
      else sort( newModels.begin(), newModels.end(), []( auto& m1, auto& m2 ){ return m1.nll < m2.nll ; } );

      INFO("SelectModel::Selected amplitude " << newModels[0].name);
      INFO("id = " << newModels[0].id);
      INFO("nll =  " << newModels[0].nll);
      INFO("chi2 = " << newModels[0].chi2);

      if( (newModels[0].chi2 + delta_chi2 > bestModel.chi2 && modelSelectionFoM=="chi2" ) || ( newModels[0].nll + delta_nll > bestModel.nll && modelSelectionFoM=="nll")   ){
          improved=false;
          INFO("Model builder found no improvement after generation " << generation);
          INFO("old chi2 = " << bestModel.chi2);
          INFO("old nll =  " << bestModel.nll);
      }
      else {
          bestModel = newModels[0]; 
          bestModel.gen = generation; 
          
          in_dir = project_dir + "/v" + to_string(bestModel.gen);
          starterModel = "model_" + to_string(bestModel.id)+".txt";
          
          vector<string> ampsToRemove = {bestModel.name};
          for(unsigned i=newModels.size()- min((int)newModels.size(),nWorstAmpsToRemoveFromList); i<newModels.size();++i){ampsToRemove.push_back(newModels[i].name);}
          makeReducedAmpList(addAmp, ampsToRemove, project_dir+"/v"+ to_string(generation) + "/" + addAmpList );
      }
      INFO("Completed generation " << generation);
      generation++;
      
  }while( generation < maxGeneration && improved );

  cout << endl << "==============================================" << endl;    
  INFO("Model building completed after generation " << generation-1);
  INFO("Best model: " << project_dir + "/v" + to_string(bestModel.gen) + "/model_" + to_string(bestModel.id));  

  analyzeModelBuilder(generation,100,project_dir);

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;   
}

