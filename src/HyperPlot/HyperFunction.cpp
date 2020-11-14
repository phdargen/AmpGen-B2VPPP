#include "AmpGen/HyperPlot/HyperFunction.h"


HyperFunction::HyperFunction() :
  _limits(HyperPoint(0), HyperPoint(0))
{




}

HyperFunction::HyperFunction(const HyperCuboid& limits) :
  _limits(limits)
{

}


///Reweight a HyperPointSet by the HyperFunction. If weights already
///exist, the existing weights are mulitplied by the HyperFunction
///evaluation. If not, a the HyperFunction evaluation is added as the
///zeroth weight.
void HyperFunction::reweightDataset(HyperPointSet& points){

  int npoints = points.size(); 
 
  for (int i = 0; i < npoints; i++){

    HyperPoint& point = points.at(i);

    double val = this->getVal(point);

    int nW = point.numWeights();

    for (int w = 0; w < nW; w++){
      double oldW = point.getWeight(w);
      double newW = oldW*val;
      point.setWeight(w, newW);
    }

    if (nW == 0) point.addWeight(val);

  }


}

void HyperFunction::setFuncLimits(const HyperCuboid& limits){
  _limits = limits;
}


TH2D HyperFunction::make2DFuncSlice   (TString name, int sliceDimX, int sliceDimY, const HyperPoint& slicePoint, int nbins) const{
  
  if (_limits.getDimension() == 0){
    ERROR_LOG << "You need to set the domain of the fuction with setFuncLimits(const HyperCuboid& limits) to draw slices" << std::endl;
  }
  if (slicePoint.getDimension() != _limits.getDimension()){
    ERROR_LOG << "Your slice point has different dimensions to the function domain" << std::endl;    
  }
  if (_limits.inVolume(slicePoint) == false){
    ERROR_LOG << "Your slice point in not within the function domain" << std::endl;    
  }  
  
  double minX = _limits.getLowCorner ().at(sliceDimX);
  double maxX = _limits.getHighCorner().at(sliceDimX);

  double minY = _limits.getLowCorner ().at(sliceDimY);
  double maxY = _limits.getHighCorner().at(sliceDimY);

  TH2D hist(name, name, nbins, minX, maxX, nbins, minY, maxY);
  
  for (int i = 1; i <= nbins; i++){
    for (int j = 1; j <= nbins; j++){
      double x = hist.GetXaxis()->GetBinCenter(i);
      double y = hist.GetYaxis()->GetBinCenter(j);
      HyperPoint point(slicePoint);
      point.at(sliceDimX) = x;
      point.at(sliceDimY) = y;
      double val = getVal(point);
      hist.SetBinContent(i,j,val);
    }
  }
  
  return hist;

}

void HyperFunction::draw2DFuncSlice   (TString path, int sliceDimX, int sliceDimY, const HyperPoint& slicePoint, int nbins) const{
  
  TH2D hist = make2DFuncSlice("temp",  sliceDimX,  sliceDimY, slicePoint, nbins);
  RootPlotter2D plotter(&hist);
  plotter.plot(path, "COLZ");

}

void HyperFunction::draw2DFuncSliceSet(TString path, int sliceDimX, int sliceDimY, int sliceSetDim, int nSlices, const HyperPoint& slicePoint, int nbins) const{

  if (_limits.getDimension() == 0){
    ERROR_LOG << "You need to set the domain of the fuction with setFuncLimits(const HyperCuboid& limits) to draw slices" << std::endl;
  }
  if (slicePoint.getDimension() != _limits.getDimension()){
    ERROR_LOG << "Your slice point has different dimensions to the function domain" << std::endl;    
  }
  if (_limits.inVolume(slicePoint) == false){
    ERROR_LOG << "Your slice point in not within the function domain" << std::endl;    
  }  

  HyperPoint slicePointCp(slicePoint);

  double min = _limits.getLowCorner ().at(sliceSetDim);
  double max = _limits.getHighCorner().at(sliceSetDim);
  double width = (max - min)/double(nSlices);
  
  for (int i = 0; i < nSlices; i++){
    double val = min + width*(i + 0.5);
    slicePointCp.at(sliceSetDim) = val;
    
    TString uniquePath = path;
    uniquePath += "_sliceNum";
    uniquePath +=  i;
    draw2DFuncSlice(uniquePath, sliceDimX, sliceDimY, slicePointCp, nbins);

  }
  

}

void HyperFunction::draw2DFuncSliceSet(TString path, int sliceDimX, int sliceDimY, int nSlices, const HyperPoint& slicePoint, int nbins) const{
  

  
  for (int i = 0; i < slicePoint.getDimension(); i++){

    if (i == sliceDimX) continue;
    if (i == sliceDimY) continue;
    
    TString thsPath = path;
    thsPath += "_scanDim";
    thsPath += i;

    draw2DFuncSliceSet(thsPath, sliceDimX, sliceDimY, i, nSlices, slicePoint, nbins);

  }
  

}

void HyperFunction::draw2DFuncSliceSet(TString path, int nSlices, const HyperPoint& slicePoint, int nbins) const{
  

  
  for (int i = 0; i < slicePoint.getDimension(); i++){
    for (int j = 0; j < slicePoint.getDimension(); j++){
      
      if (i >= j) continue;

      TString thsPath = path;
      thsPath += "_";
      thsPath += i;
      thsPath += "vs";
      thsPath += j;

      draw2DFuncSliceSet(thsPath, i, j, nSlices, slicePoint, nbins);
    }
  }
  

}



double HyperFunction::getDifference(const HyperFunction& other, const HyperPoint& point){
  double val      =       getVal(point);
  double valother = other.getVal(point);
  return val - valother;
}


void HyperFunction::fillCorrelations(TH2D& hist, const HyperFunction& other, const HyperPointSet& points){

  for (unsigned i = 0; i < points.size(); i++){
    double val      =       getVal(points.at(i));
    double valother = other.getVal(points.at(i));
    hist.Fill( val, valother, points.at(i).getWeight() );    
  }
  
  int nBinsX = hist.GetXaxis()->GetNbins();
  int nBinsY = hist.GetYaxis()->GetNbins();

  for (int i = 1; i <= nBinsX; i++){
    double sumY = 0.0;
    for (int j = 1; j <= nBinsY; j++){
      sumY += hist.GetBinContent(i,j); 
    } 
    for (int j = 1; j <= nBinsY; j++){
      double val = hist.GetBinContent(i,j);
      val = (val / sumY)*100.0;
      if (sumY == 0.0) val = 0.0;
      hist.SetBinContent(i, j, val);
    }        
  }  

}








