#ifndef LOADINGBAR_HH
#define LOADINGBAR_HH

#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>      /* printf */
#include <math.h>       /* floor */
#include <unistd.h>

class LoadingBar{
  
  int _nIterations;
  int _nSteps;
  int _prevStep;
  
  time_t _tstart;  

  public:
  
  LoadingBar(int nIterations, int nSteps = 20);
  
  void printTimeDiff(double minutes);


  void update(int i);
  
  ~LoadingBar();

};


#endif


