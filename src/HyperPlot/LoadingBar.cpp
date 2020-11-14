#include "AmpGen/HyperPlot/LoadingBar.h"

LoadingBar::LoadingBar(int nIterations, int nSteps) :
  _nIterations(nIterations),
  _nSteps     (nSteps),
  _prevStep   (-1)
{

    _tstart = time(0);


}


void LoadingBar::printTimeDiff(double minutes){

  if (minutes > 60.0){
    int hours = floor( minutes / 60.0 );
    int mins  = floor( minutes - hours*60.0 );
    std::cout << hours << "h " << mins; 
  }
  else if (minutes < 60.0 && minutes > 1.0){
    int mins    = floor( minutes );
    int seconds = floor( double( minutes - mins )*60.0 );
    std::cout << mins << "m " << seconds; 
  }
  else{
    int seconds = floor( (minutes)*600.0 );
    std::cout << seconds*0.1 << "s";
  }
}

void LoadingBar::update(int i){

  double frac = double(i+1)/double(_nIterations);
  int steps = floor(_nSteps*frac);
  
  int per = floor(100*frac);

 
  if (_prevStep != per){
  
    double delT = difftime(time(0), _tstart);
  
    double detTMins = delT / 60.0;
    double evtsPerMin = (double(i)/detTMins);
    
    double totalTimeMins = double(_nIterations)/evtsPerMin;
    double remainingTime = totalTimeMins*(1.0-frac);

    std::cout << "\r  [";
   
    for (int j = 0; j < _nSteps; j++){
      if (j < steps) std::cout << "=";
      else           std::cout << " ";
    }
    
    std::cout << "]   " << per << "% Complete. ";
    printTimeDiff(remainingTime);
    std::cout << " left.              " << std::flush;

    _prevStep = per;
    
  }
  
  if ( frac == 1.0 )  {

    double delT = difftime(time(0), _tstart);
    double detTMins = delT / 60.0;
    double evtsPerMin = (double(_nIterations)/detTMins);

    std::cout << "\r  Completed in ";
    printTimeDiff(detTMins);
    std::cout << " which is " << floor(evtsPerMin) << " iterations/min                " << std::flush << std::endl;

  }
}

LoadingBar::~LoadingBar(){

}
