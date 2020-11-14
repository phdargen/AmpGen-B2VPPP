#include "AmpGen/HyperPlot/HyperPointSet.h"

///Standard constuctor where only the dimensionality is given
///
HyperPointSet::HyperPointSet(int dimension) : 
  _dimension(dimension)
{
  WELCOME_LOG << "Hello from the HyperPointSet() Constructor";
}

///Constuctor for loading exisiting HyperPointSet from a file
///
HyperPointSet::HyperPointSet(TString path) :
  _dimension(-1)
{

  load(path);
 
}

///Constuctor for a HyperPointSet that repeats the same point
///a total of npoints times.
HyperPointSet::HyperPointSet(int npoints, const HyperPoint& point) :
  _dimension(point.getDimension())
{

  for (int i = 0; i < npoints; i++) push_back(point);

}


///Constuctor for a HyperPointSet that takes dimensionality from
///a given HyperPoint and adds that point to the HyperPointSet
HyperPointSet::HyperPointSet(const HyperPoint& point) :
  _dimension(point.getDimension())
{
  push_back(point);
}

///Constuctor for a HyperPointSet that takes dimensionality from
///a given HyperPoint and adds that point and another to the HyperPointSet
HyperPointSet::HyperPointSet(const HyperPoint& point1, const HyperPoint& point2) :
  _dimension(point1.getDimension())
{
  push_back(point1);
  push_back(point2);
}

///Constuctor for a HyperPointSet that takes dimensionality from
///a given HyperPoint and adds that point and another 2 to the HyperPointSet
HyperPointSet::HyperPointSet(const HyperPoint& point1, const HyperPoint& point2, const HyperPoint& point3) :
  _dimension(point1.getDimension())
{
  push_back(point1);
  push_back(point2);
  push_back(point3);
}

///Is this HyperPoint compatible with this HyperPointSet? (same dimension)
///
bool HyperPointSet::compatible(const HyperPoint& other, bool printError) const{

  if (_dimension == other.getDimension()) return true;

  if (printError) ERROR_LOG << "This HyperPoint is NOT compatible with this HyperPointSet";
  return false;
}

///Get a point from the HyperPointSet
///
const HyperPoint& HyperPointSet::at(int i) const      { 
  return _points.at(i);
}

///Get a point from the HyperPointSet
///
HyperPoint& HyperPointSet::at(int i)                  { 
  return _points.at(i);
}

///Get the number of points in the HyperPointSet
///
unsigned int HyperPointSet::size() const              { 
  return _points.size();
}


///Add new HyperPoint to the HyperPointSet
///
void HyperPointSet::push_back(const HyperPoint& point){ 
  if (compatible(point) == false) return;
  _points.push_back(point); 
} 


///Find the Sum of weights for all HyperPoints in the HyperPointSet
///
double HyperPointSet::getSumW() const{
  
  double sumW = 0.0;

  for(unsigned i = 0; i < size(); i++){
    sumW += at(i).getWeight();
  }
  return sumW;
}

///Find the Sum of weights squared for all HyperPoints in the HyperPointSet
///
double HyperPointSet::getSumW2() const{

  double sumW2 = 0.0;
  double w = 0.0;
  for(unsigned i = 0; i < size(); i++){
    w = at(i).getWeight();
    sumW2 += w*w;
  }
  return sumW2;

}

///Use the gramSchmidtProcess to make the first nDim
///vectors orthanormal. If there are less than this,
///only size() vectors will be returned (but they will
///be orthanormal).
HyperPointSet HyperPointSet::gramSchmidtProcess() const{

  int nPoints = size();
  if ( nPoints > getDimension() ) nPoints = getDimension();
  
  if (linearlyIndependant() == 0) {
    ERROR_LOG << "Trying to perform the Gram Schmidt Process on a set of L.D. vectors";
  }

  HyperPointSet orthogonalPoints(at(0));

  for (int i = 1; i < nPoints; i++){
    HyperPoint nextPoint = at(i);
    for (int j = 0; j < i; j++){
      nextPoint = nextPoint - orthogonalPoints.at(j).project(at(i));
    }
    orthogonalPoints.push_back(nextPoint);
  }
  
  orthogonalPoints.normalise();
  return orthogonalPoints;

}

///Check if the HyperPoints in the HyperPointSet are linearly 
///independent (this isn't implemented yet!!!)
bool HyperPointSet::linearlyIndependant() const{
  //need to implement this...
  return 1;
}

///Normalise all points in the HyperPointSet.
///
void HyperPointSet::normalise(){

  for(unsigned i = 0; i < size(); i++){
    at(i) = at(i)/at(i).norm();
  }

}

///Sort all the HyperPoints in acending order.
///First sort by the first element, then by the second
///etc.
void HyperPointSet::sort(){

  std::sort(_points.begin(), _points.end());

}

///Remove all duplicate events. IMPORTANT
///must be followed by sort()
void HyperPointSet::removeDuplicates(){

  //INFO_LOG << "Removing all duplicates from the HyperPointSet - did you remember to call sort() first?!?";

  HyperPointSet hyperPointSet(getDimension());
  
  HyperPoint currentHyperPoint = _points.at(0); 
  hyperPointSet.push_back(currentHyperPoint);

  for(unsigned int i = 1; i < size(); i++){
    if ( _points.at(i) != currentHyperPoint ) {
      currentHyperPoint = _points.at(i); 
      hyperPointSet.push_back(currentHyperPoint);
    }
  }
  
  _points.clear();
  *this = hyperPointSet;

}

///Add all the points from another HyperPointSet to this one.
///
void HyperPointSet::addHyperPointSet(const HyperPointSet& other){
  
  if (_dimension != other._dimension){
    ERROR_LOG << "Removing all duplicates from the HyperPointSet - did you remember to call sort() first?!?";

  }

  for (unsigned int i = 0; i < other.size(); i++){
    this->push_back(other.at(i));
  }

}

///Get the correlation between two variables.
///
double HyperPointSet::getCorrelation(int i, int j) const{

  WidthFinder statsi;
  WidthFinder statsj;

  for (unsigned evt = 0; evt < size(); evt++){
    statsi.add(at(evt).at(i));
    statsj.add(at(evt).at(j));
  }
  
  double meani  = statsi.mean();
  double meanj  = statsj.mean();
  double widthi = statsi.width();
  double widthj = statsj.width();

  double correlation = 0.0;

  for (unsigned evt = 0; evt < size(); evt++){
    correlation += (at(evt).at(i) - meani)*(at(evt).at(j) - meanj);
  } 
  
  correlation = correlation / (widthi * widthj * size());
  
  return correlation;
}

///Get the covarience between two variables.
///
double HyperPointSet::getCovarience(int i, int j) const{

  MeanFinder statsi;
  MeanFinder statsj;

  for (unsigned evt = 0; evt < size(); evt++){
    statsi.add(at(evt).at(i));
    statsj.add(at(evt).at(j));
  }
  
  double meani  = statsi.mean();
  double meanj  = statsj.mean();

  double correlation = 0.0;

  for (unsigned evt = 0; evt < size(); evt++){
    correlation += (at(evt).at(i) - meani)*(at(evt).at(j) - meanj);
  } 
  
  correlation = correlation / size();
  
  return correlation;
}

///Get the covarience matrix associated with the HyperPointSet.
///
TMatrixD HyperPointSet::getCovarienceMatrix() const{

  TMatrixD matrix(getDimension(), getDimension());

  for (int i = 0; i < getDimension(); i++){
    for (int j = 0; j < getDimension(); j++){
      matrix(i,j) = getCovarience(i,j);
    }    
  }
  
  return matrix;
}

///Get the correlation matrix associated with the HyperPointSet.
///
TMatrixD HyperPointSet::getCorrelationMatrix() const{

  TMatrixD matrix(getDimension(), getDimension());

  for (int i = 0; i < getDimension(); i++){
    for (int j = 0; j < getDimension(); j++){
      matrix(i,j) = getCorrelation(i,j);
    }    
  }
  
  return matrix;
}

///Print out the entire HyperPointSet.
///
void HyperPointSet::print(std::ostream& os) const{
  
  os << "---------- HyperPointSet ------------" << std::endl;
  for (unsigned i = 0; i < size(); i++){
    os << "HyperPoint " << i << ":   ";
    at(i).print(os, 0);
    os << std::endl;
  }

}

///Save the HyperPointSet to a file (opened with RECREATE).
///
void HyperPointSet::save(TString path){

  TFile* file = new TFile(path, "RECREATE");
  
  if (file == 0) {
    ERROR_LOG << "Cannot open file at " << path;
    return;
  }

  save();

  file->Write();
  file->Close();  

}

///Save the HyperPointSet to the file currently in scope.
///
void HyperPointSet::save(){
 
  TTree* tree = new TTree("HyperPointSet", "HyperPointSet");
  
  if (tree == 0) {
    ERROR_LOG << "Cannot create tree";
    return;
  }

  std::vector<double>  values;
  std::vector<double>  weights;

  tree->Branch("values" , &values );
  tree->Branch("weights", &weights);
  
  for(unsigned i = 0; i < size(); i++ ){
    weights = at(i).getWeightsVector();
    values  = at(i).getVector();
    tree->Fill();
  }
  
  tree->ResetBranchAddresses();

  tree->Write();

}

///Load exisiting HyperPointSet from a file.
///HyperPoints are appended to the end of the
//HyperPointSet, so multiple files can be loaded
void HyperPointSet::load(TString path){
  
  TFile* file = new TFile(path, "READ");
  if (file == 0) {
    ERROR_LOG << "Cannot open file at " << path;
    return;
  }

  TTree* tree = (TTree*)file->Get("HyperPointSet");
  if (tree == 0) {
    ERROR_LOG << "Cannot open tree";
    return;
  }
  
  std::vector<double>*  values  = 0;
  std::vector<double>*  weights = 0;

  tree->SetBranchAddress("values" , &values );
  tree->SetBranchAddress("weights", &weights);
  
  int nEntries = tree->GetEntries();

  for(int i = 0; i < nEntries; i++ ){
    tree->GetEntry(i);
    
    if (_dimension == -1) _dimension = values->size();
    if ( int(values->size()) != _dimension) ERROR_LOG << "You are loading HyperPoints that do not have the correct dimensionality for this HyperPointSet";

    HyperPoint point(*values);
    for(unsigned w = 0; w < weights->size(); w++){
      point.addWeight(weights->at(w));
    }
    this->push_back(point);
  }
  
  tree->ResetBranchAddresses();
  file->Close();

}

///Get a HyperCuboid that surrounds the points
//HyperCuboid HyperPointSet::getLimits() const{
//
//  HyperCuboid limits(getDimension());
//
//  for (int i = 0; i < getDimension(); i++){
//    MinMaxFinder minMax;
//    for (unsigned j = 0; j < size(); j++){
//      minMax.add( this->at(j).at(i) );
//    }
//    limits.getLowCorner() .at(i) = minMax.getMin();
//    limits.getHighCorner().at(i) = minMax.getMax();
//  }
//  
//  return limits;
//}

///Get a HyperPoint that gives the minimum in each dim
HyperPoint HyperPointSet::getMin() const{
  HyperPoint val(getDimension());

  for (int i = 0; i < getDimension(); i++){
    MinMaxFinder minMax;
    for (unsigned j = 0; j < size(); j++){
      minMax.add( this->at(j).at(i) );
    }
    val.at(i) = minMax.getMin();
  }
  return val;
}

///Get a HyperPoint that gives the maximum in each dim
HyperPoint HyperPointSet::getMax() const{
  HyperPoint val(getDimension());

  for (int i = 0; i < getDimension(); i++){
    MinMaxFinder minMax;
    for (unsigned j = 0; j < size(); j++){
      minMax.add( this->at(j).at(i) );
    }
    val.at(i) = minMax.getMax();
  }
  return val;
}


///Get the arithmatic mean of all points.
///This does not use the weights.
HyperPoint HyperPointSet::mean() const{
  HyperPoint temp(getDimension());
  for (unsigned i = 0; i < size(); i++){
    temp = temp + at(i);
  }
  return temp/double(size());
}

///Get the geometric mean of all points.
///This does not use the weights.
HyperPoint HyperPointSet::geometricMean() const{

  HyperPoint geoMean(getDimension());
  for (unsigned i = 0; i < size(); i++){
    HyperPoint temp(getDimension());
    for (int j= 0; j < getDimension(); j++){
      temp.at(j) = log( this->at(i).at(j) );
    }
    geoMean = geoMean + temp;
  }
  
  geoMean = geoMean/double(size());

  for (int j= 0; j < getDimension(); j++){
    geoMean.at(j) = exp( geoMean.at(j) );
  }
  
  return geoMean;

}

///Get the harmonic mean of all points.
///This does not use the weights.
HyperPoint HyperPointSet::harmonicMean() const{

  HyperPoint geoMean(getDimension());
  for (unsigned i = 0; i < size(); i++){
    HyperPoint temp(getDimension());
    for (int j= 0; j < getDimension(); j++){
      temp.at(j) = 1.0 / this->at(i).at(j);
    }
    geoMean = geoMean + temp;
  }
  
  geoMean = geoMean/double(size());

  for (int j= 0; j < getDimension(); j++){
    geoMean.at(j) = 1.0 / geoMean.at(j);
  }
  
  return geoMean;

}

///Destructor
///
HyperPointSet::~HyperPointSet(){
  GOODBYE_LOG << "Goodbye from the HyperPointSet() Destructor";
}
