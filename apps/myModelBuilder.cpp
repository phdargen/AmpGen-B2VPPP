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

std::map<std::string,bool> getRunningJobs(){
  std::map<std::string,bool> running_jobs; 
  auto qstat_lines = cmd("condor_q -format '%d\n' ProcId -run");
  for( auto& line : qstat_lines ) running_jobs[line] = true; 
  return running_jobs ; 
}

std::map<std::string,bool> getFailedJobs(){
    std::map<std::string,bool> failed_jobs; 
    auto qstat_lines = cmd("condor_q -format '%d\n' ProcId -hold");
    for( auto& line : qstat_lines ) failed_jobs[line] = true; 
    return failed_jobs ; 
}

std::map<std::string,bool> getJobs(){
    std::map<std::string,bool> jobs; 
    auto qstat_lines = cmd("condor_q -format '%d\n' ProcId");
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

model selectModel(string mode = "LL", string dir = "", MinuitParameterSet* mps = 0, int verbose = 0){
    
    string resultsFile = NamedParameter<string>("ResultsFile","result.root");
    TChain* chain = new TChain("Result");
    chain->Add((replaceAll(dir+"/"+resultsFile,".root","_*.root")).c_str());  
    
    double nll,chi2,sumFractions;
    unsigned status, nPar,nAmps, seed;
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
    
    const int n = chain->GetEntries();
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
    
    for (int i = 0; i<n; i++) {
        chain->GetEntry(i);
        
        h_nll->Fill(nll-min_nll);
        h_chi2->Fill(chi2);
        
        double AIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) ;
        double BIC = (nll-min_nll) + 2. * ((double) nPar - (double) min_nPar) * log(30000) ;
        
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
    
    TCanvas* c = new TCanvas("c");    
    h_nll->Draw("h");
    c->Print((dir+"/"+"n2ll.eps").c_str());  
    
    h_chi2->Draw("h");
    c->Print((dir+"/"+"chi2.eps").c_str());

    if(verbose){
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
    }
    
    if(mode=="chi2") chain->GetEntry(min_i_chi2);
    else chain->GetEntry(min_i_LL);
    model bestModel = {.name=replaceAll(mps->at(seed*2)->name(),"_Re",""), .id = seed ,.nll=nll, .chi2=chi2};
    
    return bestModel;
}



int main(int argc, char* argv[] )
{
  gStyle->SetOptStat(0);
  LHCbStyle();
  
  time_t startTime = time(0);
  OptionsParser::setArgs( argc, argv );
  bool submit = NamedParameter<bool>("submit",true);

  string project_dir  = NamedParameter<std::string>("Project","modelSel");
  string in_dir  = NamedParameter<std::string>("in_dir",".");    
  string starterModel  = NamedParameter<std::string>("StarterModel","core.txt");
  auto addAmpList = NamedParameter<string>("addAmpList","addAmpList.txt");
  string modelSelectionFoM  = NamedParameter<std::string>("modelSelectionFoM","chi2");
    
  int status=0;
  status &= system( ("mkdir -p " + project_dir).c_str() );
  if( status != 0 ) ERROR("Building project directories failed");
  
  int fitsToDo = NamedParameter<int>("fitsToDo", -1);  
  MinuitParameterSet* addAmp = new MinuitParameterSet();
  addAmp->loadFromFile(addAmpList);  
  fitsToDo = fitsToDo < 0 ? addAmp->size()/2 : fitsToDo;
    
  int generation = NamedParameter<int>("generation", 0);  
  auto maxGeneration = NamedParameter<unsigned int>("maxGeneration", 3);  
  bool improved = true;
  model bestModel = {.name="start", .id = -1,.gen=generation, .nll=99999, .chi2=99999};
    
  do{  
      cout << endl <<  "==============================================" << endl;    
      INFO("Model builder:: Starting generation " << generation);
      
      status &= system( ("mkdir -p " + project_dir+"/v"+ to_string(generation)).c_str() );
      if( status != 0 ) ERROR("Building project directories failed");
      cmd("kinit -R");
      
      const std::string command = "condor_submit o=" + starterModel + " i=" + in_dir + " d=" + project_dir+"/v"+to_string(generation) + " n=" + to_string(fitsToDo) + " submitModelBuilder.sub";  
      if(submit){  
          auto status2 = cmd( command.c_str() );  
          for( auto& s : status2) {
              INFO( s ); 
              //vector<string>cluster_id; 
              //if(s.find("cluster"))cluster_id = split(s,' ');
              //INFO(cluster_id[cluster_id.size()-1]);
          }
          auto nJobs = getJobs(); 
          INFO(nJobs.size() << "/" << fitsToDo << " jobs submitted");
      }
 
      unsigned int nCompleted = 0 ;     
      unsigned int nFailed = 0 ;       
      do {
          nCompleted = 0 ;     
          nFailed = 0 ;     
          for (unsigned i=0; i< fitsToDo; i++) {
              if(readFile(project_dir+"/v"+to_string(generation)+"/log_"+to_string(i)+".txt"))nCompleted++;
          }
          if( nCompleted < fitsToDo ){
             auto runningJobs = getRunningJobs(); 
             auto failedJobs = getFailedJobs(); 
             nFailed = failedJobs.size();
              
             cout << endl << " Time since start " << (time(0) - startTime)/60.0 << " min." << endl;
             INFO(fitsToDo - runningJobs.size() - nCompleted - failedJobs.size() << " jobs idle");
             INFO(runningJobs.size() << " jobs running");
             INFO(nCompleted << " jobs completed");
             INFO(failedJobs.size() << " jobs failed");   
              
             if(failedJobs.size()>fitsToDo/4){
                 ERROR("Too many failed jobs");
                 //killAllJobs();
                 return 0;
             }
             
             sleep(300);
             continue;
          }
      } while( nCompleted + nFailed < fitsToDo );
      
      auto newModel = selectModel(modelSelectionFoM,project_dir+"/v"+to_string(generation), addAmp); 
      INFO("SelectModel::Selected amplitude " << newModel.name);
      INFO("nll =  " << newModel.nll);
      INFO("chi2 = " << newModel.chi2);

      if(newModel.chi2 + 0.05 > bestModel.chi2){
          improved=false;
          INFO("Model builder found no improvement after generation " << generation);
          INFO("old chi2 = " << bestModel.chi2);
          INFO("old nll =  " << bestModel.nll);
      }
      else {
          bestModel = newModel; 
          bestModel.gen = generation; 
          in_dir = project_dir + "/v" + to_string(bestModel.gen);
          starterModel = "model_" + to_string(bestModel.id)+".txt";
          //addAmpList = makeReducedAmpList();
      }
      INFO("Completed generation " << generation);
      generation++;
      
  }while( generation < maxGeneration && improved );

  cout << endl << "==============================================" << endl;    
  INFO("Model building completed after generation " << generation-1);
  INFO("Best model: " << project_dir + "/v" + to_string(bestModel.gen) + "/model_" + to_string(bestModel.id));  
    
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;   
}

