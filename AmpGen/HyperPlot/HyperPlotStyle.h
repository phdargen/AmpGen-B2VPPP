/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Holds a n-dim vector of names that can be associtated
 *  to a n-dim Histogram (for axis titles etc). Also possible
 *  to store the units.
 **/

 
#ifndef HYPERPLOTSTYLE_HH
#define HYPERPLOTSTYLE_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"

// Root includes
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TColor.h"


// std includes
#include <iomanip>
#include <iostream>

class HyperPlotStyle {
  
  static int* palette_pulls;


  static void LHCbStyle(Bool_t colzPlot=kFALSE, Int_t NCont=25);
  static void makePalettes(); 

  public:
  
  HyperPlotStyle(); 
  
  static void init();

  static void setPalette(TString name);

  ~HyperPlotStyle(); 

}; 


#endif















