/*
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Global functions that are useful for various things
 */

#ifndef GLOBAL_FUNCTIONS_HH
#define GLOBAL_FUNCTIONS_HH

#include <string>
#include <sstream>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <cmath>
#include "TMath.h"
#include <limits>

///Check if a double is inf or -inf
///
inline bool isInf(double val)
{
  //bool decision = val >= std::numeric_limits<double>::min() && val <= std::numeric_limits<double>::max();
  if ( (val != 0.0) && (val + val == val) && (val == val) ) return 1;
  return 0;
}

///Check if a double is NaN
///
inline bool isNaN(double val)
{
  if (val != val) return 1;
  return 0;
}
 
///Check if a double is NaN or ±inf
///
inline bool isError(double val){
  return isInf(val) || isNaN(val);
}

///Convert things into strings
///
template <class typeToString> std::string makeString(const typeToString& thingToString)
{
  std::ostringstream ss;
  ss << thingToString;
  return ss.str();
}

///Convert an array into a vector
///
template <class arrayType> std::vector<arrayType> gVectorFromArray(arrayType* a, int size){

  std::vector<arrayType> returnVector;
  for (int i = 0; i < size; i++) returnVector.push_back(a[i]);
  return returnVector;

}

///Solve a quadratic equation
///
inline void gQuadraticSolver(double a, double b, double c, double& sol1, double& sol2){

  double root = b*b - 4.0*a*c;
  if (root >= 0.0){
    sol1 = (-b + sqrt(root))/(2.0*a);
    sol2 = (-b - sqrt(root))/(2.0*a);
  }
  else{
    sol1 = -99999999.9;
    sol2 = -99999999.9;
  }

}

///Divide two TH1D's (don't trust ROOT)
///
inline TH1D gDivideTH1D(TH1D const* hist1, TH1D const* hist2, TString name){

  TH1D outputHist(*hist1 );
  outputHist.SetName(name);

  for(int bin = 1; bin <= hist1->GetNbinsX(); bin++){
    double cont1  = hist1->GetBinContent(bin);
    double cont2  = hist2->GetBinContent(bin);
    double err1   = hist1->GetBinError(bin);
    double err2   = hist2->GetBinError(bin);
    double frac1  = err1/cont1;
    double frac2  = err2/cont2;
    double cont   = cont1/cont2;
    double err    = cont*sqrt(frac1*frac1 + frac2*frac2);
    if (cont != cont) {cont = 0.0;err = 0.0;}
    if (cont2 == 0.0) {cont = 0.0;err = 0.0;}
    if (err != err) {err = 0.0;}
    outputHist.SetBinContent(bin, cont);
    outputHist.SetBinError(bin, err);
  }
  
  return outputHist;
  
}

///Divide two TH2D's (don't trust ROOT)
///
inline TH2D gDivideTH2D(TH2D const* hist1, TH2D const* hist2, TString name){

  TH2D outputHist(*hist1 );
  outputHist.SetName(name);

  for(int binX = 1; binX <= hist1->GetNbinsX(); binX++){
    for(int binY = 1; binY <= hist1->GetNbinsY(); binY++){
    double cont1  = hist1->GetBinContent(binX, binY);
    double cont2  = hist2->GetBinContent(binX, binY);
    double err1   = hist1->GetBinError(binX, binY);
    double err2   = hist2->GetBinError(binX, binY);
    double frac1  = err1/cont1;
    double frac2  = err2/cont2;
    double cont   = cont1/cont2;
    double err    = cont*sqrt(frac1*frac1 + frac2*frac2);
    if (cont != cont) {cont = 0.0;err = 0.0;}
    if (cont2 == 0.0) {cont = 0.0;err = 0.0;}
    if (err != err) {err = 0.0;}
    outputHist.SetBinContent(binX, binY, cont);
    outputHist.SetBinError(binX, binY, err);
    }
  }
  
  return outputHist;
  
}

///Make a hard copy of a TH1D - probably not needed
///
inline TH1D hardCopyTH1D(TH1D const* hist, TString name){
  TH1D temp(name,name,hist->GetNbinsX(),hist->GetXaxis()->GetXbins()->GetArray());
  for(int bin = 1; bin <= hist->GetNbinsX(); bin++){
    temp.SetBinContent(bin, hist->GetBinContent(bin));
    temp.SetBinError  (bin, hist->GetBinError  (bin));
  }
  return temp;
}

///Make a pull histogram between two TH1D's
///
inline TH1D gPullTH1D(TH1D const* hist1, TH1D const* hist2, TString name){

  TH1D outputHist( *hist1 );
  outputHist.SetName(name);

  for(int bin = 1; bin <= hist1->GetNbinsX(); bin++){
    double cont1  = hist1->GetBinContent(bin);
    double cont2  = hist2->GetBinContent(bin);
    double err1   = hist1->GetBinError(bin);
    double err2   = hist2->GetBinError(bin);
    double err1minus2 = sqrt(err1*err1 + err2*err2);
    double cont   = (cont1 - cont2)/err1minus2;
    double err    = 0.0;
    outputHist.SetBinContent(bin, cont);
    outputHist.SetBinError(bin, err);
  }
  
  return outputHist;
  
}

///Multiply two TH1D's
///
inline TH1D gMultiplyTH1D(TH1D const* hist1, TH1D const* hist2, TString name){


  TH1D outputHist( *hist1 );

  for(int bin = 1; bin <= hist1->GetNbinsX(); bin++){
    double cont1  = hist1->GetBinContent(bin);
    double cont2  = hist2->GetBinContent(bin);
    double err1   = hist1->GetBinError(bin);
    double err2   = hist2->GetBinError(bin);
    double frac1  = err1/cont1;
    double frac2  = err2/cont2;
    double cont   = cont1*cont2;
    double err    = cont*sqrt(frac1*frac1 + frac2*frac2);
    outputHist.SetBinContent(bin, cont);
    outputHist.SetBinError(bin, err);
  }
  
  outputHist.SetName(name);

  return outputHist;

}

///Multiply two TH2D's
///
inline TH2D gMultiplyTH2D(TH2D const* hist1, TH2D const* hist2, TString name){

  TH2D outputHist( *hist1 );

  for(int binX = 1; binX <= hist1->GetNbinsX(); binX++){
    for(int binY = 1; binY <= hist1->GetNbinsY(); binY++){
      double cont1  = hist1->GetBinContent(binX, binY);
      double cont2  = hist2->GetBinContent(binX, binY);
      double err1   = hist1->GetBinError(binX, binY);
      double err2   = hist2->GetBinError(binX, binY);
      double frac1  = err1/cont1;
      double frac2  = err2/cont2;
      double cont   = cont1*cont2;
      double err    = cont*sqrt(frac1*frac1 + frac2*frac2);
      outputHist.SetBinContent(binX, binY, cont);
      outputHist.SetBinError(binX, binY, err);
    }
  }

  outputHist.SetName(name);
  
  return outputHist;

}

///Square root of a TH1D
///
inline TH1D gSqrtTH1D(TH1D const* hist1, TString name){


  TH1D outputHist( *hist1 );

  for(int bin = 1; bin <= hist1->GetNbinsX(); bin++){
    double cont  = hist1->GetBinContent(bin);
    double err   = hist1->GetBinError(bin);
    outputHist.SetBinContent(bin, sqrt(cont));
    outputHist.SetBinError(bin, err/(2.0*sqrt(cont)));
  }
  
  outputHist.SetName(name);

  return outputHist;

}

///Add two TH1D's
///
inline TH1D gAddTH1D(TH1D const* hist1, TH1D const* hist2, TString name){

  TH1D outputHist(name,name, hist1->GetNbinsX(), hist1->GetXaxis()->GetXbins()->GetArray() );

  for(int bin = 1; bin <= hist1->GetNbinsX(); bin++){
    double cont1  = hist1->GetBinContent(bin);
    double cont2  = hist2->GetBinContent(bin);
    double err1   = hist1->GetBinError(bin);
    double err2   = hist2->GetBinError(bin);
    double cont   = cont1 + cont2;
    double err    = sqrt(err1*err1 + err2*err2);
    outputHist.SetBinContent(bin, cont);
    outputHist.SetBinError(bin, err);
  }
  
  return outputHist;
}

///Function for printing the status of a loop
///
inline bool printInterationStatus(int interation, int total){
  
  int frequency = pow( 10, floor(log10(total)) - 1);

  double nPrints = double(total)/double(frequency); 
  
  if      (nPrints > 50) frequency *= 5;
  else if (nPrints > 20) frequency *= 2;
  
  if (interation < 4){
    //std::cout << "Interation " << interation << " of " << total << std::endl;
    return true;
  }
  else if (interation % frequency == 0) {
    //std::cout << "Interation " << interation << " of " << total << std::endl;
    return true;
  }
  return false;
}

///Add two TH2D's
///
inline TH2D gAddTH2D(TH2D const* hist1, TH2D const* hist2, TString name){

  TH2D outputHist(*hist1);
  outputHist.SetName(name);

  for(int binX = 1; binX <= hist1->GetNbinsX(); binX++){
    for(int binY = 1; binY <= hist1->GetNbinsY(); binY++){
      double cont1  = hist1->GetBinContent(binX, binY);
      double cont2  = hist2->GetBinContent(binX, binY);
      double err1   = hist1->GetBinError(binX, binY);
      double err2   = hist2->GetBinError(binX, binY);
      double cont   = cont1 + cont2;
      double err    = sqrt(err1*err1 + err2*err2);
      outputHist.SetBinContent(binX, binY, cont);
      outputHist.SetBinError(binX, binY, err);
    }
  }
  
  return outputHist;
}

///Minus two TH1D's (hist1 - hist2)
///
inline TH1D gMinusTH1D(TH1D const* hist1, TH1D const* hist2, TString name){


  TH1D outputHist(name,name, hist1->GetNbinsX(), hist1->GetXaxis()->GetXbins()->GetArray() );
  outputHist.SetName(name);

  for(int bin = 1; bin <= hist1->GetNbinsX(); bin++){
    double cont1  = hist1->GetBinContent(bin);
    double cont2  = hist2->GetBinContent(bin);
    double err1   = hist1->GetBinError(bin);
    double err2   = hist2->GetBinError(bin);
    double cont   = cont1 - cont2;
    double err    = sqrt(err1*err1 + err2*err2);
    outputHist.SetBinContent(bin, cont);
    outputHist.SetBinError(bin, err);
  }
  
  return outputHist;
}

///Minus two TH2D's (hist1 - hist2)
///
inline TH2D gMinusTH2D(TH2D const* hist1, TH2D const* hist2, TString name){

  TH2D outputHist(*hist1);
    outputHist.SetName(name);

  for(int binX = 1; binX <= hist1->GetNbinsX(); binX++){
    for(int binY = 1; binY <= hist1->GetNbinsY(); binY++){
      double cont1  = hist1->GetBinContent(binX, binY);
      double cont2  = hist2->GetBinContent(binX, binY);
      double err1   = hist1->GetBinError(binX, binY);
      double err2   = hist2->GetBinError(binX, binY);
      double cont   = cont1 - cont2;
      double err    = sqrt(err1*err1 + err2*err2);
      outputHist.SetBinContent(binX, binY, cont);
      outputHist.SetBinError(binX, binY, err);
    }
  }
  
  return outputHist;
}

///Convert degrees to radians
///
inline double degreesToRadians(double degrees) { return (degrees*TMath::Pi())/180.0; }

///Convert radians to degrees
///
inline double radiansToDegrees(double radians) { return (radians/TMath::Pi())*180.0; }

///Shift angle to the domain [0, 360] / [0,2pi] (use degrees = 0 for radians)
///
inline double shiftAngleToDomain(double angle, int degrees = 1){
  
  if (degrees == 1){
    if ((angle >= 0.0) && (angle < 360.0)) {
      return angle;
    }
    else {
      double val = angle - 360.0*floor(angle/360.0);
      return val;
    }
  }
  if (degrees == 0){
    if ((angle >= 0.0) && (angle < 2.0*TMath::Pi())) {
      return angle;
    }
    else {
      double val = angle - 2.0*TMath::Pi()*floor(angle/(2.0*TMath::Pi()));
      return val;
    }
  }
  
  return -9999.9;

}

#endif
