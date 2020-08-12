#ifndef _LHCBSTYLE_H
#define _LHCBSTYLE_H


#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TColor.h"

#include <iomanip>
#include <iostream>


/*
//

TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb official plots style");
//TPaveText *lhcbName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
TText *lhcbLabel = new TText();
TLatex *lhcbLatex = new TLatex();
TPaveText *lhcb7TeVPrelimL = new TPaveText();
TPaveText *lhcb0_9TeVPrelimL = new TPaveText();

void LHCbStyle(Bool_t colzPlot=kFALSE, Int_t NCont=25)
{
  // ////////////////////////////////////////////////////////////////////
  // // PURPOSE:
  // //
  // // This header file defines a reasonable style for (black-and-white)
  // // "publication quality" ROOT plots. The default settings contain
  // // many features that are either not desirable for printing on white
  // // paper or impair the general readibility of plots.
  // //
  // // USAGE:
  //
  // Simply include the line
  //   #include "LHCbStyle.h"
  // at the beginning of your root program (and make sure LHCbStyle.h is
  // in a location accessible to the compiler). Then add the line
  //   lhcbStyle();
  // somewhere at the beginning of your main().
  //
  // SOME COMMENTS:
  //
  // Statistics and fit boxes:
  //
  // "Decorative" items around the histogram are kept to a minimum.
  // In particular there is no box with statistics or fit information.
  // You can easily change this either by editing your private copy
  // of this style file or by calls to "gStyle" in your macro.
  // For example,
  //   gStyle->SetOptFit(1011);
  // will add some fit information.
  //
  // Font:
  //
  // The font is chosen to be 62, i.e.helvetica-bold-r-normal with
  // precision 2. Font is of course a matter of taste, but most people
  // will probably agree that Helvetica bold gives close to optimal
  // readibility in presentations. It appears to be the ROOT default,
  // and since there are still some features in ROOT that simply won't
  // respond to any font requests, it is the wise choice to avoid
  // ugly font mixtures on the same plot... The precision of the font (2)
  // is chosen in order to have a rotatable and scalable font. Be sure
  // to use true-type fonts! I.e.
  // Unix.*.Root.UseTTFonts: true  in your .rootrc file.
  //
  // "Landscape histograms":
  //
  // The style here is designed for more or less quadratic plots.
  // For very long histograms, adjustements are needed. For instance,
  // for a canvas with 1x5 histograms:
  //  TCanvas* c1 = new TCanvas("c1", "L0 muons", 600, 800);
  //  c1->Divide(1,5);
  // adaptions like the following will be needed:
  //  gStyle->SetTickLength(0.05,"x");
  //  gStyle->SetTickLength(0.01,"y");
  //  gStyle->SetLabelSize(0.15,"x");
  //  gStyle->SetLabelSize(0.1,"y");
  //  gStyle->SetStatW(0.15);
  //  gStyle->SetStatH(0.5);
  //
  ////////////////////////////////////////////////////////////////////

  //gROOT->Reset();

  //TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb official plots style");

  // use helvetica-bold-r-normal, precision 2 (rotatable)
Int_t lhcbFont = 62;
// line thickness
Int_t lhcbWidth = 3;

// use plain black on white colors
lhcbStyle->SetFrameBorderMode(0);
lhcbStyle->SetCanvasBorderMode(0);
lhcbStyle->SetPadBorderMode(0);
lhcbStyle->SetPadColor(0);
lhcbStyle->SetCanvasColor(0);
lhcbStyle->SetStatColor(0);
lhcbStyle->SetPalette(1);
lhcbStyle->SetTitleColor(0);
lhcbStyle->SetFillColor(0);

// set the paper & margin sizes
lhcbStyle->SetPaperSize(20,26);
lhcbStyle->SetPadTopMargin(0.05);
lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots!!
lhcbStyle->SetPadBottomMargin(0.16);
lhcbStyle->SetPadLeftMargin(0.14);

// use large fonts
lhcbStyle->SetTextFont(lhcbFont);
lhcbStyle->SetTextSize(0.08);
lhcbStyle->SetLabelFont(lhcbFont,"x");
lhcbStyle->SetLabelFont(lhcbFont,"y");
lhcbStyle->SetLabelFont(lhcbFont,"z");
lhcbStyle->SetLabelSize(0.05,"x");
lhcbStyle->SetLabelSize(0.05,"y");
lhcbStyle->SetLabelSize(0.05,"z");
lhcbStyle->SetTitleFont(lhcbFont);
lhcbStyle->SetTitleSize(0.06,"x");
lhcbStyle->SetTitleSize(0.06,"y");
lhcbStyle->SetTitleSize(0.06,"z");
lhcbStyle->SetTitleColor(1,"x");

// use bold lines and markers
lhcbStyle->SetLineWidth(lhcbWidth);
lhcbStyle->SetFrameLineWidth(lhcbWidth);
lhcbStyle->SetHistLineWidth(lhcbWidth);
lhcbStyle->SetFuncWidth(lhcbWidth);
lhcbStyle->SetGridWidth(lhcbWidth);
lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
lhcbStyle->SetMarkerStyle(8);
lhcbStyle->SetMarkerSize(0.5);

// label offsets
lhcbStyle->SetLabelOffset(0.015);

// by default, do not display histogram decorations:
lhcbStyle->SetOptStat(0);
//lhcbStyle->SetOptStat(1110);  // show only nent, mean, rms
lhcbStyle->SetOptTitle(0);
lhcbStyle->SetOptFit(0);
//lhcbStyle->SetOptFit(1011); // show probability, parameters and errors

// look of the statistics box:
lhcbStyle->SetStatBorderSize(1);
lhcbStyle->SetStatFont(lhcbFont);
lhcbStyle->SetStatFontSize(0.05);
lhcbStyle->SetStatX(0.9);
lhcbStyle->SetStatY(0.9);
lhcbStyle->SetStatW(0.25);
lhcbStyle->SetStatH(0.15);

// put tick marks on top and RHS of plots
lhcbStyle->SetPadTickX(1);
lhcbStyle->SetPadTickY(1);
//lhcbStyle->SetPadTickX(0);
//lhcbStyle->SetPadTickY(0);

// histogram divisions: only 5 in x to avoid label overlaps
lhcbStyle->SetNdivisions(505,"x");
lhcbStyle->SetNdivisions(510,"y");

if (colzPlot) {
lhcbStyle->SetNumberContours(NCont<100?NCont:100);
lhcbStyle->SetPadRightMargin(0.15);
lhcbStyle->SetPaintTextFormat(".3f");
lhcbStyle->SetHistMinimumZero(); // Extra line to allow zero bin to be plotted in 2D Text histo
}


gROOT->SetStyle("lhcbStyle");
gROOT->ForceStyle();


//TPaveText *lhcbName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
//lhcbName->SetFillColor(0);
//lhcbName->SetTextAlign(12);
//lhcbName->SetBorderSize(0);
//lhcbName->AddText("LHCb");

//TText *lhcbLabel = new TText();
lhcbLabel->SetTextFont(lhcbFont);
lhcbLabel->SetTextColor(1);
lhcbLabel->SetTextSize(0.04);
lhcbLabel->SetTextAlign(12);

//TLatex *lhcbLatex = new TLatex();
lhcbLatex->SetTextFont(lhcbFont);
lhcbLatex->SetTextColor(1);
lhcbLatex->SetTextSize(0.04);
lhcbLatex->SetTextAlign(12);

//  TPaveText *lhcb7TeVPrelimL = new TPaveText(lhcbStyle->GetPadLeftMargin() + 0.05,
//  0.78 - lhcbStyle->SetPadTopMargin(0.05),
//  lhcbStyle->GetPadLeftMargin() + 0.30,
//  0.88 - lhcbStyle->SetPadTopMargin(0.05),
//  "BRNDC");


lhcb7TeVPrelimL->SetX1NDC(lhcbStyle->GetPadLeftMargin() + 0.05);
lhcb7TeVPrelimL->SetX2NDC(lhcbStyle->GetPadLeftMargin() + 0.30);
lhcb7TeVPrelimL->SetY1NDC(0.78 - lhcbStyle->GetPadTopMargin());
lhcb7TeVPrelimL->SetY2NDC(0.98 - lhcbStyle->GetPadTopMargin());
lhcb7TeVPrelimL->SetFillColor(0);
lhcb7TeVPrelimL->SetTextAlign(12);
lhcb7TeVPrelimL->SetBorderSize(0);
lhcb7TeVPrelimL->SetTextSize(0.07);
lhcb7TeVPrelimL->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 7 TeV Data}}");

lhcb0_9TeVPrelimL->SetX1NDC(lhcbStyle->GetPadLeftMargin() + 0.05);
lhcb0_9TeVPrelimL->SetX2NDC(lhcbStyle->GetPadLeftMargin() + 0.30);
lhcb0_9TeVPrelimL->SetY1NDC(0.78 - lhcbStyle->GetPadTopMargin());
lhcb0_9TeVPrelimL->SetY2NDC(0.98 - lhcbStyle->GetPadTopMargin());
lhcb0_9TeVPrelimL->SetFillColor(0);
lhcb0_9TeVPrelimL->SetTextAlign(12);
lhcb0_9TeVPrelimL->SetBorderSize(0);
lhcb0_9TeVPrelimL->SetTextSize(0.07);
lhcb0_9TeVPrelimL->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 900 GeV Data}}");

}

*/

// all users - please change the name of this file to lhcbStyle.C
// Commits to lhcbdocs svn of .C files are not allowed

void LHCbStyle(Bool_t colzPlot=kFALSE, Int_t NCont=25)
{

  // define names for colours
  Int_t black  = 1;
  Int_t red    = 2;
  Int_t green  = 3;
  Int_t blue   = 4;
  Int_t yellow = 5;
  Int_t magenta= 6;
  Int_t cyan   = 7;
  Int_t purple = 9;


  ////////////////////////////////////////////////////////////////////
  // PURPOSE:
  //
  // This macro defines a standard style for (black-and-white)
  // "publication quality" LHCb ROOT plots.
  //
  // USAGE:
  //
  // Include the lines
  //   gROOT->ProcessLine(".L lhcbstyle.C");
  //   lhcbStyle();
  // at the beginning of your root macro.
  //
  // Example usage is given in myPlot.C
  //
  // COMMENTS:
  //
  // Font:
  //
  // The font is chosen to be 132, this is Times New Roman (like the text of
  //  your document) with precision 2.
  //
  // "Landscape histograms":
  //
  // The style here is designed for more or less square plots.
  // For longer histograms, or canvas with many pads, adjustements are needed.
  // For instance, for a canvas with 1x5 histograms:
  //  TCanvas* c1 = new TCanvas("c1", "L0 muons", 600, 800);
  //  c1->Divide(1,5);
  //  Adaptions like the following will be needed:
  //  gStyle->SetTickLength(0.05,"x");
  //  gStyle->SetTickLength(0.01,"y");
  //  gStyle->SetLabelSize(0.15,"x");
  //  gStyle->SetLabelSize(0.1,"y");
  //  gStyle->SetStatW(0.15);
  //  gStyle->SetStatH(0.5);
  //
  // Authors: Thomas Schietinger, Andrew Powell, Chris Parkes, Niels Tuning
  // Maintained by Editorial board member (currently Niels)
  ///////////////////////////////////////////////////////////////////

  // Use times new roman, precision 2
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06;

  // use plain black on white colors
  gROOT->SetStyle("Plain");
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");

  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);
  lhcbStyle->SetLegendFont(132);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  //int colors[8] = {0,5,7,3,6,2,4,1};
  //lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.14);

  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  //lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerStyle(8);
  //lhcbStyle->SetMarkerSize(1.0);
  lhcbStyle->SetMarkerSize(0.5); //Sam found these too big

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0);
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);

  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");

  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  // add LHCb label
  TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                      0.87 - gStyle->GetPadTopMargin(),
                                      gStyle->GetPadLeftMargin() + 0.20,
                                      0.95 - gStyle->GetPadTopMargin(),
                                      "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);

  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  std::cout << "-------------------------" << std::endl;
  std::cout << "Set LHCb Style - Feb 2012" << std::endl;
  std::cout << "-------------------------" << std::endl;

}



#endif
