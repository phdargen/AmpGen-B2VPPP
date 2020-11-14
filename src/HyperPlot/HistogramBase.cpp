#include "AmpGen/HyperPlot/HistogramBase.h"

#include "AmpGen/HyperPlot/StatisticsFinder.h"


///Construct a histogram base with a specified number of bins
///
HistogramBase::HistogramBase(int nBins) :
  _nBins(nBins),
  _binContents(nBins+1,0.0),
  _sumW2      (nBins+1,0.0),
  _min(-999.9),
  _max(-999.9),
  _minDensity(-999.9),
  _maxDensity(-999.9)
{
  WELCOME_LOG << "Good day from the HistogramBase() Constructor"; 
}

///Update the number of bins and set all contents to zero
///
void HistogramBase::resetBinContents(int nBins){
  _nBins = nBins;
  _binContents = std::vector<double>(nBins+1,0.0);
  _sumW2       = std::vector<double>(nBins+1,0.0);
}

void HistogramBase::clear(){
  for (unsigned i = 0; i < _binContents.size(); i++){
    _binContents.at(i) = 0.0;
    _sumW2      .at(i) = 0.0;
  }
}



///
void HistogramBase::reserveCapacity(int nElements){
  _binContents.reserve(nElements+1);
  _sumW2      .reserve(nElements+1);
}


///Merge one HistogramBase with another
///
void HistogramBase::merge( const HistogramBase& other ){
  
  double overflowCont = _binContents.at(_nBins);
  overflowCont += other._binContents.at(other._nBins);
  
  double overflowSumW2 = _sumW2.at(_nBins);
  overflowSumW2 += other._sumW2.at(other._nBins);  
  
  int capacityNeeded = _nBins + other._nBins + 1;

  reserveCapacity(capacityNeeded);

  _binContents.at(_nBins) = other._binContents.at(0);
  _sumW2      .at(_nBins) = other._sumW2      .at(0);
  
  for (int i = 1; i < other._nBins; i++){
    _binContents.push_back( other._binContents.at(i) );
    _sumW2      .push_back( other._sumW2      .at(i) );    
  }

  _binContents.push_back( overflowCont  );
  _sumW2      .push_back( overflowSumW2 );

  _nBins += other._nBins;

}



/// Save the contents, sumw2, and bin numbers in a TTree
/// to an open TFile
void HistogramBase::saveBase(){

  TTree tree("HistogramBase", "HistogramBase");

  int    binNumber  = -1;
  double binContent = 0.0;
  double sumW2      = 0.0;

  tree.Branch("binNumber" , &binNumber );
  tree.Branch("binContent", &binContent);
  tree.Branch("sumW2"     , &sumW2     );

  for(unsigned int bin = 0; bin < _binContents.size(); bin++ ){
    binNumber  = bin;
    binContent = _binContents.at(bin);
    sumW2      = _sumW2      .at(bin);
    tree.Fill();
  }
  
  tree.Write();

}

/// Save the contents, sumw2, and bin numbers in a TTree
/// into the ROOT file specified (opened using RECREATE)
void HistogramBase::saveBase(TString filename){

  TFile file(filename, "RECREATE");

  saveBase();

  file.Write();
  file.Close();

}

/// Load the contents, sumw2, and bin numbers from a TTree
/// in the ROOT file specified (opened using READ)
void HistogramBase::loadBase(TString filename){

  TFile file(filename, "READ");
  TTree* tree = (TTree*)file.Get("HistogramBase");

  int    binNumber  = -1;
  double binContent = 0.0;
  double sumW2      = 0.0;
  
  int nEntries = tree->GetEntries();

  tree->SetBranchAddress("binNumber" , &binNumber );
  tree->SetBranchAddress("binContent", &binContent);
  tree->SetBranchAddress("sumW2"     , &sumW2     );
  
  this->resetBinContents(nEntries - 1);


  for(int ent = 0; ent < nEntries; ent++){
    tree->GetEntry(ent);
    
    VERBOSE_LOG << "Bin = " << binNumber << "       Content = " << binContent << "        SumW2 = " << sumW2;

    _binContents.at((int)binNumber) = binContent;
    _sumW2      .at((int)binNumber) = sumW2;
  }
  
  file.Close();

}

///Destructor
///
HistogramBase::~HistogramBase(){
  GOODBYE_LOG << "Goodbye from the HistogramBase() Constructor"; 
}

///Fill the specifed bin with a specifed weight.
///Note that the bin numbers run from 0 to nBins-1 inclusive
/// so that nBins is overflow/underflow
void HistogramBase::fillBase(int binNum, double weight){
  binNum = checkBinNumber(binNum);
  _binContents[binNum] += weight;
  _sumW2[binNum]       += weight*weight;
}

///Check if it's a valid bin number, if not
///return the overflow/underflow bin
int HistogramBase::checkBinNumber(int bin) const{
  if(bin == -1){
    //This is also a valid way to add to the underflow/overflow bin with no error message
    return _nBins;
  }
  if(bin > _nBins || bin < 0) {
    ERROR_LOG << "Bin " << bin << " does not exist! Adding to underflow/overflow bin " << _nBins;
    return _nBins;
  }
  return bin;
}

///If the user hasn't specified a minimum using setMin, then
///loop over the bins and find the minimum
double HistogramBase::getMin() const{
  if (_min != -999.9) return _min;
  MinMaxFinder stats;
  for (int i = 0; i < _nBins; i++){
    stats.add( _binContents[i] );
  }
  return stats.getMin();
}

///If the user hasn't specified a minimum using setMax, then
///loop over the bins and find the maximum
double HistogramBase::getMax() const{
  if (_max != -999.9) return _max;
  MinMaxFinder stats;
  for (int i = 0; i < _nBins; i++){
    stats.add( _binContents[i] );
  }
  return stats.getMax();
}

///If the user hasn't specified a minimum using setMinDensity, then
///loop over the bins and find the minimum density
double HistogramBase::getMinDensity() const{
  if (_minDensity != -999.9) return _minDensity;
  MinMaxFinder stats;
  for (int i = 0; i < _nBins; i++){
    stats.add( getFrequencyDensity(i) );
  }
  return stats.getMin();
}

///If the user hasn't specified a maximum using setMaxDensity, then
///loop over the bins and find the maximum density
double HistogramBase::getMaxDensity() const{
  if (_maxDensity != -999.9) return _maxDensity;
  MinMaxFinder stats;
  for (int i = 0; i < _nBins; i++){
    stats.add( getFrequencyDensity(i) );
  }
  return stats.getMax();
}

///Divide this hitogram by another.
///Histograms must have the same number of bins
void HistogramBase::divide(const HistogramBase& other){
 
  if (other._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for(unsigned int i = 0; i < _binContents.size(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);
    double frac1sq  = var1/(binCont1*binCont1);
    double frac2sq  = var2/(binCont2*binCont2);
    
    double content = binCont1 / binCont2;
    double var     = content*content*(frac1sq + frac2sq);
    
    if (binCont2 == 0.0) content = 0.0;

    _binContents.at(i) = content;
    _sumW2.at(i)       = var;

  }
  
  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;
}

///Multiply this hitogram by another.
///Histograms must have the same number of bins
void HistogramBase::multiply(const HistogramBase& other){

  if (other._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for(unsigned int i = 0; i < _binContents.size(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);
    double frac1sq  = var1/(binCont1*binCont1);
    double frac2sq  = var2/(binCont2*binCont2);
    
    double content = binCont1 * binCont2;
    double var     = content*content*(frac1sq + frac2sq);

    _binContents.at(i) = content;
    _sumW2.at(i)       = var;

  }

  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;

}

///Add this hitogram to another.
///Histograms must have the same number of bins
void HistogramBase::add(const HistogramBase& other){

  if (other._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for(unsigned int i = 0; i < _binContents.size(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);
    
    double content = binCont1 + binCont2;
    double var     = var1 + var2;

    _binContents.at(i) = content;
    _sumW2.at(i)       = var;

  }

  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;

}

///Subtract other histogram from this one.
///Histograms must have the same number of bins
void HistogramBase::minus(const HistogramBase& other){

  if (other._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for(unsigned int i = 0; i < _binContents.size(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);
    
    double content = binCont1 - binCont2;
    double var     = var1 + var2;

    _binContents.at(i) = content;
    _sumW2.at(i)       = var;

  }

  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;

}

///Find pulls between this histogram and another.
///Histograms must have the same number of bins
void HistogramBase::pulls(const HistogramBase& other){

  if (other._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for( int i = 0; i < getNBins(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);

    double content = binCont1 - binCont2;
    double var     = var1 + var2;

    _binContents.at(i) = content/sqrt(var);
    _sumW2.at(i)       = 0.0;

  }
  
  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;  
}


///Find pulls between two histograms.
///Histograms must have the same number of bins
void HistogramBase::pulls(const HistogramBase& other1, const HistogramBase& other2){

  if (other1._binContents.size() != _binContents.size() || other2._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for( int i = 0; i < getNBins(); i++){

    double binCont1 = other1._binContents.at(i);
    double binCont2 = other2._binContents.at(i);
    double var1     = other1._sumW2.at(i);
    double var2     = other2._sumW2.at(i);

    double content = binCont1 - binCont2;
    double var     = var1 + var2;

    _binContents.at(i) = content/sqrt(var);
    _sumW2.at(i)       = 0.0;

  }
  
  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;  
}


void HistogramBase::asymmetry(const HistogramBase& other){

  if (other._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for( int i = 0; i < getNBins(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);

    double content = (binCont1 - binCont2)/(binCont1 + binCont2);
    
    double dzdx  = (2.0*binCont2)/((binCont1 + binCont2)*(binCont1 + binCont2));
    double dzdy  = (2.0*binCont1)/((binCont1 + binCont2)*(binCont1 + binCont2));
    double varsq = dzdx*dzdx*var1*var1 + dzdy*dzdy*var2*var2;

    _binContents.at(i) = content;
    _sumW2.at(i)       = varsq;

  }
  
  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;  

}

void HistogramBase::asymmetry(const HistogramBase& other1, const HistogramBase& other2){

  if (other1._binContents.size() != _binContents.size() || other2._binContents.size() != _binContents.size()){
    ERROR_LOG << "Trying to divide histograms with different numbers of bins. Doing nothing.";
    return;
  }

  for( int i = 0; i < getNBins(); i++){

    double binCont1 = other1._binContents.at(i);
    double binCont2 = other2._binContents.at(i);
    double var1     = other1._sumW2.at(i);
    double var2     = other2._sumW2.at(i);

    double content = (binCont1 - binCont2)/(binCont1 + binCont2);
    
    double dzdx  = (2.0*binCont2)/((binCont1 + binCont2)*(binCont1 + binCont2));
    double dzdy  = (2.0*binCont1)/((binCont1 + binCont2)*(binCont1 + binCont2));
    double varsq = dzdx*dzdx*var1*var1 + dzdy*dzdy*var2*var2;

    _binContents.at(i) = content;
    _sumW2.at(i)       = varsq;

  }
  
  //reset min and max to default values
  _min = -999.9;
  _max = -999.9;  

}



///Randomise the histogram within it's Gaussian errors using
///the given seed.
///If value is < 0.0 then content is set to 0.0.
///errors remain the same before and after randomisation.
void HistogramBase::randomiseWithinErrors(int seed){

  TRandom3 random(seed);

  for( int i = 0; i < getNBins(); i++){

    double mean  = _binContents.at(i);
    double sigma = sqrt(_sumW2.at(i));

    double newContent = random.Gaus(mean, sigma);
    if (newContent < 0.0) newContent = 0.0;
    _binContents.at(i) = newContent;
  }

}

///Calculate the chi2 between this histogram and another, and 
///use this to find a p-value using the specified degrees of 
///freedom. If ndof isn't specified, nBins is used.
///Histograms must have the same number of bins
double HistogramBase::pvalue(const HistogramBase& other, int ndof) const{
  
  if (ndof < 0.0) ndof = getNBins();
  double chi2 = this->chi2(other);

  return TMath::Prob(chi2, ndof);

}

///Calculate the significance of the chi2 value obtained (in #sigma) given that
///it follows a chi2 distribution for ndof degrees of freedom. If
///ndof is not given it is assumed that ndof == nbins
double HistogramBase::chi2sig(const HistogramBase& other, int ndof) const{
  
  double pval = pvalue(other, ndof);
  double sigma = TMath::NormQuantile( pval/2.0 );
  return fabs(sigma);

}


///Calculate the chi2 between this histogram and another.
///Histograms must have the same number of bins
double HistogramBase::chi2(const HistogramBase& other) const{

  double chi2 = 0.0;

  for( int i = 0; i < getNBins(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);
    
    double content = binCont1 - binCont2;
    double var     = var1 + var2;

    double thisTerm = (content * content) / var;

    chi2 += thisTerm;

  }

  return chi2;

}

///Draw the distribution of pulls between two histograms.
///They are only drawn between +- pmLimits sigma, using
///nBins. The plot is saved with the name 'name'. 
void HistogramBase::drawPullHistogram(const HistogramBase& other, TString name, int nBins, double pmLimits) const{

  TH1D pulls ("pulls", "pulls", nBins, -pmLimits, pmLimits);
  pulls.GetXaxis()->SetTitle("pull");
  pulls.GetYaxis()->SetTitle("frequency");
  for( int i = 0; i < getNBins(); i++){

    double binCont1 =       _binContents.at(i);
    double binCont2 = other._binContents.at(i);
    double var1     =       _sumW2.at(i);
    double var2     = other._sumW2.at(i);
    
    double content = binCont1 - binCont2;
    double var     = var1 + var2;

    double pull = (content) / sqrt(var);
    pulls.Fill(pull);

  }

  RootPlotter1D plotter(&pulls);
  plotter.plot(name);

}

///Calcualte the integral (sum of bin contents, excluding the undeflow/overflow)
///
double HistogramBase::integral() const{
  double integral = 0.0;
  for(int i = 0; i < getNBins() ; i++){
    integral += _binContents.at(i);
  }
  return integral;
}

///Calcualte the error on the integral (excluding the undeflow/overflow)
///
double HistogramBase::integralError() const{
  double var = 0.0;
  for(int i = 0; i < getNBins() ; i++){
    var += _sumW2.at(i);
  }
  return sqrt(var);
}

///Normalise sum of bin contents to specified area.
///Errors are also scaled.
void HistogramBase::normalise(double area){

  double divisor  = integral()/area;
  double divisor2 = divisor*divisor;

  for( int i = 0; i < getNBins(); i++){
    _binContents.at(i) = _binContents.at(i)/divisor;
    _sumW2      .at(i) = _sumW2      .at(i)/divisor2;
  }

}

///Set the content of a bin (leaves error unchanged)
///
void HistogramBase::setBinContent(int bin, double val){
  bin = checkBinNumber(bin);
  _binContents[bin] = val;
}

///Set the error of a bin
///
void HistogramBase::setBinError  (int bin, double val){
  bin = checkBinNumber(bin);
  _sumW2[bin] = val*val;
}

///Get the content of a bin
///
double HistogramBase::getBinContent(int bin) const{
  bin = checkBinNumber(bin); 
  return _binContents[bin];
}

///Get the error of a bin
///
double HistogramBase::getBinError  (int bin) const{
  bin = checkBinNumber(bin);
  return sqrt(_sumW2[bin]);
}

///Get the bin volume. This is a virual function, as 
///the bin volume depends on the binning used (which isn't
///defined here).
double HistogramBase::getBinVolume(int bin) const{
  bin = checkBinNumber(bin);
  return 1.0;
}

///Get the frequency density in a bin - this is just
///the bin content divided by the bin volume
double HistogramBase::getFrequencyDensity(int bin) const{

  return getBinContent(bin)/getBinVolume(bin);

}

///Replace all the bin contents with frequency density.
///
void HistogramBase::makeFrequencyDensity(){
  
  for(int i = 0; i < _nBins; i++){
    double binVolume = this->getBinVolume(i);
    _binContents[i] = _binContents[i] / binVolume;
    _sumW2[i] = _sumW2[i] / (binVolume*binVolume);
  }

}

///Print the contents and errors of all the bins to the screen.
///
void HistogramBase::print(){

  for(int i = 0; i < _nBins; i++){
    INFO_LOG << "Bin Content " << i << ": " << _binContents[i] << "      SumW2: " << _sumW2[i];
  }
  INFO_LOG << "Overflow: " << _binContents[_nBins];

}

