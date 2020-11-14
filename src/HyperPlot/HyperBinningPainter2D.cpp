#include "AmpGen/HyperPlot/HyperBinningPainter2D.h"

/** Construct a 2D HyperBinningPainter for a given HyperVolumeBinning and HyperPointSet.
This option will draw the HyperPointSet as dots over the binning scheme
*/
HyperBinningPainter2D::HyperBinningPainter2D(BinningBase* binning, HyperPointSet* hyperPoints) :
  HyperBinningPainter(0),
  _binning(binning),
  _hyperPoints(hyperPoints)
{
  if (getBinning().getDimension() != 2){
    ERROR_LOG << "You have given a 2D painter a different dimensionality of binning and/or data. I'll probably crash soon";  
  }
  if (hyperPoints != 0){
    if(hyperPoints->getDimension() !=2){
      ERROR_LOG << "You have given a 2D painter a different dimensionality of binning and/or data. I'll probably crash soon"; 
    }
  }
  WELCOME_LOG << "Good day from the HyperBinningPainter2D() Constructor";  
}

/** Construct a 2D HyperBinningPainter for a given HyperBinningHistogram */
HyperBinningPainter2D::HyperBinningPainter2D(HyperHistogram* histogram) :
  HyperBinningPainter(histogram),
  _binning(0),
  _hyperPoints(0)
{
  if (getBinning().getDimension() != 2){
    ERROR_LOG << "You have given a 2D painter a different dimensionality of binning and/or data. I'll probably crash soon";  
  }
  WELCOME_LOG << "Good day from the HyperBinningPainter2D() Constructor";  
}

/** get the binning  (works for either constructor) */
const BinningBase& HyperBinningPainter2D::getBinning(){
  
  if (_binning != 0) return *_binning;
  return _histogram->getBinning();

}

/** add the bin edges to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawBinEdges(RootPlotter2D* plotter){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawBinEdge(plotter, i);
}

/** add the bin edges to the Plotter (for one HyperVolume) */
void HyperBinningPainter2D::drawBinEdge(RootPlotter2D* plotter, int bin){
  
  for (int i = 0; i < getBinning().getBinHyperVolume(bin).size(); i++){
    HyperCuboid temp (getBinning().getBinHyperVolume(bin).getHyperCuboid(i));
    drawBinEdge(plotter, &temp);
  }
}

/** add the bin edges to the Plotter (for one HyperCuboid) */
void HyperBinningPainter2D::drawBinEdge(RootPlotter2D* plotter, HyperCuboid* bin){

  TLine* leftLine   = new TLine(bin->getLowCorner().at(0) , bin->getLowCorner().at(1), bin->getLowCorner().at(0) , bin->getHighCorner().at(1));
  TLine* rightLine  = new TLine(bin->getHighCorner().at(0), bin->getLowCorner().at(1), bin->getHighCorner().at(0), bin->getHighCorner().at(1));
  TLine* topLine    = new TLine(bin->getLowCorner().at(0) , bin->getHighCorner().at(1), bin->getHighCorner().at(0) , bin->getHighCorner().at(1));
  TLine* bottomLine = new TLine(bin->getLowCorner().at(0) , bin->getLowCorner().at(1), bin->getHighCorner().at(0) , bin->getLowCorner().at(1));
  
  leftLine  ->SetLineWidth(1);
  rightLine ->SetLineWidth(1);
  topLine   ->SetLineWidth(1);
  bottomLine->SetLineWidth(1);

  leftLine  ->SetLineColor(kBlack);
  rightLine ->SetLineColor(kBlack);
  topLine   ->SetLineColor(kBlack);
  bottomLine->SetLineColor(kBlack);

  //
  leftLine  ->SetLineStyle(3);
  //rightLine ->SetLineStyle(3);
  //topLine   ->SetLineStyle(3);
  bottomLine->SetLineStyle(3);

  plotter->addObject(leftLine);
  //plotter->addObject(rightLine);
  //plotter->addObject(topLine);
  plotter->addObject(bottomLine);

}


/** draw edges between bins with different contents */
void HyperBinningPainter2D::drawBinEdges2(RootPlotter2D* plotter){

  double minBinWidthX = 10e50;
  double minBinWidthY = 10e50;

  for(int i = 0; i < getBinning().getNumBins(); i++){
    HyperCuboid vol = getBinning().getBinHyperVolume(i).getLimits();
    double binWidthX = vol.getWidth(0);
    double binWidthY = vol.getWidth(1);
    if (binWidthX < minBinWidthX) minBinWidthX = binWidthX;
    if (binWidthY < minBinWidthY) minBinWidthY = binWidthY;
  }

  for(int i = 0; i < getBinning().getNumBins(); i++) drawBinEdge2(plotter, i, minBinWidthX, minBinWidthY);
}

/** draw edges between bins with different contents */
void HyperBinningPainter2D::drawBinEdge2(RootPlotter2D* plotter, int bin, double minWidX, double minWidY){
  
  for (int i = 0; i < getBinning().getBinHyperVolume(bin).size(); i++){
    HyperCuboid temp (getBinning().getBinHyperVolume(bin).getHyperCuboid(i));
    drawBinEdge2(plotter, &temp, minWidX, minWidY);
  }
}

/** draw edges between bins with different contents */
void HyperBinningPainter2D::drawBinEdge2(RootPlotter2D* plotter, HyperCuboid* bin, double minWidX, double minWidY){

  drawBinEdge2(plotter, bin, 0 , minWidX, minWidY);
  drawBinEdge2(plotter, bin, 1 , minWidX, minWidY);
  drawBinEdge2(plotter, bin, 2 , minWidX, minWidY);
  drawBinEdge2(plotter, bin, 3 , minWidX, minWidY);
  
}

/** draw edges between bins with different contents */
void HyperBinningPainter2D::drawBinEdge2(RootPlotter2D* plotter, HyperCuboid* bin, int edge, double minWidX, double minWidY){
  
  double binContent = _histogram->getVal( bin->getCenter() );
  TLine* line = 0;

  if (edge == 0) {  //left edge

    HyperPoint cornerLow (bin->getLowCorner().at(0)-minWidX*0.5, bin->getLowCorner ().at(1)+minWidY*0.05);
    HyperPoint cornerHigh(bin->getLowCorner().at(0)-minWidX*0.5, bin->getHighCorner().at(1)-minWidY*0.05);   
    int    binNumLow  = getBinning().getBinNum( cornerLow  ); 
    int    binNumHigh = getBinning().getBinNum( cornerHigh ); 
    
    if (binNumLow  == -1 ) return;
    if (binNumHigh == -1 ) return;
    
    double binConLow  = _histogram->getBinContent(binNumLow );
    double binConHigh = _histogram->getBinContent(binNumHigh);
    
    HyperCuboid binLow  = getBinning().getBinHyperVolume(binNumLow ).at(0);
    HyperCuboid binHigh = getBinning().getBinHyperVolume(binNumHigh).at(0);

    if (binNumLow == binNumHigh && binContent == binConLow) return;
    else if (binNumLow == binNumHigh && binContent != binConLow) {
      line = new TLine( bin->getLowCorner().at(0), bin->getLowCorner ().at(1), bin->getLowCorner().at(0), bin->getHighCorner().at(1) );
    }
    else if (binConLow != binContent){
      line = new TLine( bin->getLowCorner().at(0), bin->getLowCorner ().at(1), bin->getLowCorner().at(0), binLow.getHighCorner().at(1) );
    }
    else if (binConHigh != binContent){
      line = new TLine( bin->getLowCorner().at(0), binHigh.getLowCorner ().at(1), bin->getLowCorner().at(0), bin->getHighCorner().at(1) );
    }

  }
  
  if (edge == 1) {  //right edge
  
    HyperPoint cornerLow (bin->getHighCorner().at(0)+minWidX*0.5, bin->getLowCorner ().at(1)+minWidY*0.05);
    HyperPoint cornerHigh(bin->getHighCorner().at(0)+minWidX*0.5, bin->getHighCorner().at(1)-minWidY*0.05);   
    int    binNumLow  = getBinning().getBinNum( cornerLow  ); 
    int    binNumHigh = getBinning().getBinNum( cornerHigh ); 
  
    if (binNumLow  == -1 ) return;
    if (binNumHigh == -1 ) return;
  
    double binConLow  = _histogram->getBinContent(binNumLow );
    double binConHigh = _histogram->getBinContent(binNumHigh);
    
    HyperCuboid binLow  = getBinning().getBinHyperVolume(binNumLow ).at(0);
    HyperCuboid binHigh = getBinning().getBinHyperVolume(binNumHigh).at(0);
  
    if (binNumLow == binNumHigh && binContent == binConLow) return;
    else if (binNumLow == binNumHigh && binContent != binConLow) {
      line = new TLine( bin->getHighCorner().at(0), bin->getLowCorner ().at(1), bin->getHighCorner().at(0), bin->getHighCorner().at(1) );
    }
    else if (binConLow != binContent){
      line = new TLine( bin->getHighCorner().at(0), bin->getLowCorner ().at(1), bin->getHighCorner().at(0), binLow.getHighCorner().at(1) );
    }
    else if (binConHigh != binContent){
      line = new TLine( bin->getHighCorner().at(0), binHigh.getLowCorner ().at(1), bin->getHighCorner().at(0), bin->getHighCorner().at(1) );
    }
  
  }
  
  if (edge == 2) {  //top edge
  
    HyperPoint cornerLow (bin->getLowCorner ().at(0)+minWidX*0.05, bin->getHighCorner().at(1)+minWidY*0.5);
    HyperPoint cornerHigh(bin->getHighCorner().at(0)-minWidX*0.05, bin->getHighCorner().at(1)+minWidY*0.5);   
    int    binNumLow  = getBinning().getBinNum( cornerLow  ); 
    int    binNumHigh = getBinning().getBinNum( cornerHigh ); 
  
    if (binNumLow  == -1 ) return;
    if (binNumHigh == -1 ) return;
  
    double binConLow  = _histogram->getBinContent(binNumLow );
    double binConHigh = _histogram->getBinContent(binNumHigh);
    
    HyperCuboid binLow  = getBinning().getBinHyperVolume(binNumLow ).at(0);
    HyperCuboid binHigh = getBinning().getBinHyperVolume(binNumHigh).at(0);
  
    if (binNumLow == binNumHigh && binContent == binConLow) return;
    else if (binNumLow == binNumHigh && binContent != binConLow) {
      line = new TLine( bin->getLowCorner ().at(0), bin->getHighCorner().at(1), bin->getHighCorner().at(0), bin->getHighCorner().at(1) );
    }
    else if (binConLow != binContent){
      line = new TLine( bin->getLowCorner ().at(0), bin->getHighCorner().at(1), binLow.getHighCorner().at(0), bin->getHighCorner().at(1)  );
    }
    else if (binConHigh != binContent){
      line = new TLine( binHigh.getLowCorner ().at(0), bin->getHighCorner().at(1), bin->getHighCorner().at(0), bin->getHighCorner().at(1)  );
    }
  
  }
  
  if (edge == 3) {  //bottom edge
  
    HyperPoint cornerLow (bin->getLowCorner ().at(0)+minWidX*0.05, bin->getLowCorner().at(1)-minWidY*0.5);
    HyperPoint cornerHigh(bin->getHighCorner().at(0)-minWidX*0.05, bin->getLowCorner().at(1)-minWidY*0.5);   
    int    binNumLow  = getBinning().getBinNum( cornerLow  ); 
    int    binNumHigh = getBinning().getBinNum( cornerHigh ); 
  
    if (binNumLow  == -1 ) return;
    if (binNumHigh == -1 ) return;     
    double binConLow  = _histogram->getBinContent(binNumLow );
    double binConHigh = _histogram->getBinContent(binNumHigh);
    
    HyperCuboid binLow  = getBinning().getBinHyperVolume(binNumLow ).at(0);
    HyperCuboid binHigh = getBinning().getBinHyperVolume(binNumHigh).at(0);
  
    if (binNumLow == binNumHigh && binContent == binConLow) return;
    else if (binNumLow == binNumHigh && binContent != binConLow) {
      line = new TLine( bin->getLowCorner ().at(0), bin->getLowCorner().at(1), bin->getHighCorner().at(0), bin->getLowCorner().at(1) );
    }
    else if (binConLow != binContent){
      line = new TLine( bin->getLowCorner ().at(0), bin->getLowCorner().at(1), binLow.getHighCorner().at(0), bin->getLowCorner().at(1)  );
    }
    else if (binConHigh != binContent){
      line = new TLine( binHigh.getLowCorner ().at(0), bin->getLowCorner().at(1), bin->getHighCorner().at(0), bin->getLowCorner().at(1)  );
    }
  
  }
  
  if (line == 0) return;

  line  ->SetLineWidth(1);
  line  ->SetLineColor(kBlack);
  line  ->SetLineStyle(1);

  plotter->addObject(line);


}


/** Code stolen from the THistPainter - choose the colour to
make the bin based on the colour scale of the TH2D */
int HyperBinningPainter2D::getFillColour(double binContents){
  double min = _histogram->getMin();
  double max = _histogram->getMax();
  if (_density == true){
    min = _histogram->getMinDensity();
    max = _histogram->getMaxDensity();
  }
  
  if (binContents < min || binContents > max) return kWhite;

  int ndivz = gStyle->GetNumberContours();
  double dz = max - min;
  double scale = ndivz/dz;
  int ncolors  = gStyle->GetNumberOfColors();

  int color = int(0.01+(binContents-min)*scale);
  int theColor = int((color+0.99)*double(ncolors)/double(ndivz));
  if (theColor > ncolors-1) theColor = ncolors-1;

  return gStyle->GetColorPalette(theColor);
}

/** Draw bin numbers on the all the bins */
void HyperBinningPainter2D::drawBinNumbers(RootPlotter2D* plotter){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawBinNumbers(plotter, i);
}


/** Draw bin number on a single bin */
void HyperBinningPainter2D::drawBinNumbers(RootPlotter2D* plotter, int bin){
  HyperPoint center = getBinning().getBinHyperVolume(bin).getAverageCenter();
  TString label = ""; label += bin;
  plotter->addText( label ,center.at(0), center.at(1), 2, 2, 0.02, 0);

}

/** Draw bin contents on the all the bins */
void HyperBinningPainter2D::drawBinCont(RootPlotter2D* plotter){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawBinCont(plotter, i);
}


/** Draw bin content on a single bin */
void HyperBinningPainter2D::drawBinCont(RootPlotter2D* plotter, int bin){
  HyperPoint center = getBinning().getBinHyperVolume(bin).getAverageCenter();
  double binCont = _histogram->getBinContent(bin );
  TString label = ""; label += binCont;
  plotter->addText( label ,center.at(0), center.at(1), 2, 2, 0.02, 0);

}


/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBins(RootPlotter2D* plotter, bool hashNeg ){
  for(int i = 0; i < getBinning().getNumBins(); i++) drawFilledBin(plotter, i, hashNeg);
}

/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBin(RootPlotter2D* plotter, int bin, bool hashNeg){
  
  double binContent = _histogram->getBinContent(bin);
  if (_density == true ) binContent = _histogram->getFrequencyDensity(bin);
  
  bool neg = false;
  if (hashNeg == true && binContent < 0.0) {
    binContent = -binContent;
    neg = true;
  }

  for (int i = 0; i < getBinning().getBinHyperVolume(bin).size(); i++){
    HyperCuboid temp (getBinning().getBinHyperVolume(bin).getHyperCuboid(i));
    drawFilledBin(plotter, &temp, binContent);
    if (neg){
      drawFilledBin(plotter, &temp, 1, 3344);
    }
  }

}


/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBin(RootPlotter2D* plotter, HyperCuboid* bin, double binContents){
  
  int fillColor = getFillColour(binContents);
  int fillStyle = 1001;

  drawFilledBin(plotter, bin, fillColor, fillStyle);

}

/** add filled bins to the Plotter (for all HyperVolumes) */
void HyperBinningPainter2D::drawFilledBin(RootPlotter2D* plotter, HyperCuboid* bin, int fillColor, int fillStyle){

  TBox* box = new TBox(bin->getLowCorner().at(0) , bin->getLowCorner().at(1), bin->getHighCorner().at(0) , bin->getHighCorner().at(1));
  box->SetFillColor( fillColor );
  box->SetLineWidth(0.0);
  box->SetFillStyle( fillStyle );
  box->SetLineColor( fillColor );
  plotter->addObject(box);

}


/** If HyperPoints provided, add to TH2D */
void HyperBinningPainter2D::addHyperPoints(TH2D* histogram){
  for(unsigned int i = 0; i < _hyperPoints->size(); i++){
    histogram->Fill(_hyperPoints->at(i).at(0), _hyperPoints->at(i).at(1));    
  }  
}

/** Draw the HyperBinningHistogram  */
void HyperBinningPainter2D::draw(TString path, TString option){
  
  bool drawBinEd1     = option.Contains("Edges1" );
  bool drawBinEd2     = option.Contains("Edges2" );
  bool drawBinNums    = option.Contains("BinNums");
  bool drawHashedNeg  = option.Contains("HashNeg");
  bool drawBinContOp  = option.Contains("Text");

  double x_min = getBinning().getMin(0);
  double x_max = getBinning().getMax(0);
  double y_min = getBinning().getMin(1);
  double y_max = getBinning().getMax(1);
  
  TString xtitle =  _histogram->getNames().getAxisString(0);
  TString ytitle =  _histogram->getNames().getAxisString(1);

  TH2D* histogram = new TH2D("2DHyperBinningPlot", "2DHyperBinningPlot", 100, x_min, x_max, 100, y_min, y_max);
  histogram->GetXaxis()->SetTitle( xtitle );
  histogram->GetYaxis()->SetTitle( ytitle );

  double min = _histogram->getMin();
  double max = _histogram->getMax();
  if (_density == true){
    min = _histogram->getMinDensity();
    max = _histogram->getMaxDensity();
  }
  
  if (drawHashedNeg){
    min = min<0.0?0.0:min;
    max = fabs(max)>fabs(min)?fabs(max):fabs(min);
  }

  if (_histogram   != 0){
    //if you want the scale to be correct on the zaxis, we need to fill the
    //histogram with events between min-max (annoying)
    for(int i = 1; i <= 100; i++) for (int j = 0; j <= 100; ++j) histogram->SetBinContent(i, j, max); 
    histogram->SetBinContent(2, 2, min); 
    histogram->SetMaximum(max);
    histogram->SetMinimum(min); 
  } 

  histogram->SetContour(gStyle->GetNumberContours());
  //std::cout << "Am setting number of contours to " << gStyle->GetNumberContours() << std::endl;
  RootPlotter2D plotter(histogram);

  if (_hyperPoints != 0) addHyperPoints(histogram);
  if (_binning     != 0) drawBinEdges  (&plotter);
  if (_histogram   != 0) {
    TBox* box = new TBox(getBinning().getMin(0) , getBinning().getMin(1), getBinning().getMax(0) , getBinning().getMax(1));
    box->SetFillColor( 0 );
    box->SetLineWidth(0.0);
    plotter.addObject(box)  ;
    drawFilledBins(&plotter, drawHashedNeg);
    if (drawBinEd1    ) drawBinEdges  (&plotter);
    if (drawBinEd2    ) drawBinEdges2 (&plotter);
    if (drawBinNums   ) drawBinNumbers(&plotter);
    if (drawBinContOp ) drawBinCont   (&plotter);

  }

  plotter.plot(path, "COLZ");

  delete histogram;
}

/** Constructor  */
HyperBinningPainter2D::~HyperBinningPainter2D(){
  
}
