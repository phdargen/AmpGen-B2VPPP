#include "AmpGen/HyperPlot/HyperBinningDiskRes.h"


///The empty constuctor. Must call the load function to associate this
///object to a file
HyperBinningDiskRes::HyperBinningDiskRes() :
  _file(0),
  _writeable(false),
  _tree(0),
  _cuboid(0),
  _linkedBins(new std::vector<int>()),
  _volumeNumber(-1),
  _currentEntry(-1),
  _treePrimVol(0),
  _primVolNum(-1)
{
  WELCOME_LOG << "Hello from the HyperBinningDiskRes() Constructor";
}

///The copy constuctor. Essentially just calls the empty constuctor, then
///loads from the same file as the other binning.
HyperBinningDiskRes::HyperBinningDiskRes(const HyperBinningDiskRes& other) :
  _file(0),
  _writeable(false),
  _tree(0),
  _cuboid(0),
  _linkedBins(new std::vector<int>()),
  _volumeNumber(-1),
  _currentEntry(-1),
  _treePrimVol(0),
  _primVolNum(-1)
{ 
  if (other._writeable == true){
    other._file->Write();
  }  
  load(other._file->GetName(), "READ");
}


///Set the dimension of the HyperBinningDiskRes. This can only be 
///called once, when it is known what the dimension is i.e. either
///when loading from a file, or the first time a HyperVolume is 
///added to an empty file
void HyperBinningDiskRes::setDimension(int dim){
  
  if (getDimension() == 0){
    HyperBinning::setDimension(dim);
    _cuboid = HyperCuboid(dim);
  }

}

///Get a HyperVolume from the tree and load it into memory
///
void HyperBinningDiskRes::getEntry(int volumeNumber) const{
  if (_tree == 0){
    ERROR_LOG << "HyperBinningDiskRes::getEntry - tree doesn't exist" << std::endl; 
    return;
  }
  if ( _currentEntry != volumeNumber ){
    _currentEntry = volumeNumber;
    _tree->GetEntry(volumeNumber);
  }
}

///Get a HyperVolume from its volume number
HyperVolume HyperBinningDiskRes::getHyperVolume(int volumeNumber) const{
  getEntry(volumeNumber);
  return HyperVolume(_cuboid);
} 

///Get all HyperVolumes linked to a specific volume number i.e.
///from the binning hiearcy
std::vector<int> HyperBinningDiskRes::getLinkedHyperVolumes( int volumeNumber ) const{
  getEntry(volumeNumber);
  return *_linkedBins;
}

///Create a clone of the object and return a pointer to it.
BinningBase* HyperBinningDiskRes::clone() const{
  
  BinningBase* ret = dynamic_cast<BinningBase*>(new HyperBinningDiskRes());

  if (_file != 0) {

    TString filenm = filename();
    _file->cd();

    if (_writeable == true){
      _file->Write();
    }
    
    ret->load(filenm, "READ");

  }

  return ret;

} 


///Add a primary volume number
///
void HyperBinningDiskRes::addPrimaryVolumeNumber(int volumeNumber){

  if (_writeable == false){
    ERROR_LOG << "HyperBinningDiskRes::addPrimaryVolumeNumber - Cannot write, you must have opened in READ mode" << std::endl;
    return;
  }
  if (_treePrimVol == 0){
    ERROR_LOG << "HyperBinningDiskRes::addPrimaryVolumeNumber - Cannot find tree" << std::endl;
    return;
  }  

  _primVolNum = volumeNumber;
  _treePrimVol->Fill();
}



///Get the number of HyperVolumes (this is not the number of bins)
///
int HyperBinningDiskRes::getNumHyperVolumes() const{
  if (_tree == 0){
    return 0;
  }
  return _tree->GetEntries();
}  


///Add a HyperVolume and its linked bins to the HyperBinning
///
bool HyperBinningDiskRes::addHyperVolume(const HyperVolume& hyperVolume, std::vector<int> linkedVolumes){
  
  if (_writeable == false){
    ERROR_LOG << "HyperBinningDiskRes::addHyperVolume - Cannot write, you must have opened in READ mode" << std::endl;
    return false;
  }

  //If dimension hasn't been set, and it's a writeable, we can now
  //set the dimesnion and create the tree

  if (getDimension() == 0 && _writeable == true){
    setDimension(hyperVolume.getDimension());
    createHyperBinningTree();
  }
  
  if (_tree == 0){
    ERROR_LOG << "HyperBinningDiskRes::addHyperVolume - Cannot find tree" << std::endl;
    return false;
  }

  int nVolumes = getNumHyperVolumes();
  for (int i = 0; i < hyperVolume.size(); i++){
    *_linkedBins  = linkedVolumes;
    _cuboid       = hyperVolume.at(i);
    _volumeNumber = nVolumes;
    _tree->Fill();
  }
  
  updateCash();
  
  return true;
}


///Get the number of primary volumes
///
int HyperBinningDiskRes::getNumPrimaryVolumes  () const{
  if (_treePrimVol == 0) return 0;
  return _treePrimVol->GetEntries();
}

///Get the primary volume numbers 
///
int HyperBinningDiskRes::getPrimaryVolumeNumber(int i) const{
  if (_treePrimVol == 0) {
    ERROR_LOG << "HyperBinningDiskRes::getPrimaryVolumeNumber - No tree found" << std::endl;
  }
  _treePrimVol->GetEntry(i);
  return _primVolNum;
}


void HyperBinningDiskRes::loadHyperBinningTree(){

  _tree = dynamic_cast<TTree*>(_file->Get("HyperBinning"));

  if (_tree == 0){
    ERROR_LOG << "Could not open TTree in HyperBinningDiskRes::load()";
    return;
  }
  
  //Figure out how many dimensions there are from the tree
  setDimension( getHyperBinningDimFromTree(_tree) );
  
  _tree->SetBranchAddress("binNumber" , &_volumeNumber);
  _tree->SetBranchAddress("linkedBins", &_linkedBins   );

  for (int i = 0; i < getDimension(); i++){
    TString lowCornerName  = "lowCorner_" ; lowCornerName  += i; 
    TString highCornerName = "highCorner_"; highCornerName += i;   
    _tree->SetBranchAddress(lowCornerName , &_cuboid.getLowCorner ().at(i) );
    _tree->SetBranchAddress(highCornerName, &_cuboid.getHighCorner().at(i) );
  }
  
  //_tree->SetMaxVirtualSize(1000);
  _currentEntry = -1;
  updateCash();

}



void HyperBinningDiskRes::loadPrimaryVolumeTree(){

  _treePrimVol = dynamic_cast<TTree*>( _file->Get("PrimaryVolumeNumbers") );
  
  if (_treePrimVol == 0){
    ERROR_LOG << "HyperBinningDiskRes::loadPrimaryVolumeNumbers - could not open TTree";
    return;
  }

  _treePrimVol->SetBranchAddress("volumeNumber", &_primVolNum);

}

void HyperBinningDiskRes::createHyperBinningTree(){
  
  _file->cd();
  _tree = new TTree("HyperBinning", "HyperBinning");

  if (_treePrimVol == 0){
    ERROR_LOG << "HyperBinningDiskRes::createHyperBinningTree - could not open TTree";
    return;
  }

  //Figure out how many dimensions there are from the tree
  int dim = getDimension();
  if (dim == 0){
    ERROR_LOG << "The dimesion has not yet been set, so I cannot createHyperBinningTree!!" << std::endl;
    return;
  }

  _tree->Branch("binNumber" , &_volumeNumber );
  _tree->Branch("linkedBins", &_linkedBins   );

  for (int i = 0; i < dim; i++){
    TString lowCornerName  = "lowCorner_" ; lowCornerName  += i; 
    TString highCornerName = "highCorner_"; highCornerName += i;   
    _tree->Branch(lowCornerName , &_cuboid.getLowCorner ().at(i) );
    _tree->Branch(highCornerName, &_cuboid.getHighCorner().at(i) );
  }
  
  //random access is really slow. Thought this might help
  //_tree->SetMaxVirtualSize(1000);

  _currentEntry = -1;
  updateCash();

}

void HyperBinningDiskRes::createPrimaryVolumeTree(){

  _file->cd();
  _treePrimVol = new TTree("PrimaryVolumeNumbers", "PrimaryVolumeNumbers");

  if (_treePrimVol == 0){
    ERROR_LOG << "HyperBinningDiskRes::createPrimaryVolumeTree - could not open TTree";
    return;
  }

  _treePrimVol->Branch("volumeNumber", &_primVolNum);

}

///Load HyperBinningDiskRes from a file
///
void HyperBinningDiskRes::load(TString filename, TString option){
  
  _writeable = false;
  if (option == "UPDATE" || option == "RECREATE") {
    _writeable = true;
  }
  
  if (_file != 0) _file->Close();

  _file = new TFile(filename, option);
  
  if (_file == 0){
    ERROR_LOG << "Could not open TFile in HyperBinningDiskRes::load(" << filename << ")";
    return;
  }
  
  if (option == "UPDATE" || option == "READ"){
    loadHyperBinningTree   ();
    loadPrimaryVolumeTree  ();
  }
  else if (option == "RECREATE"){
    //if we're opening a new file, we don't know the dimension yet.
    //wait until the dimension is set explicitly with setDimension().
    //createHyperBinningTree   (); 
    createPrimaryVolumeTree  ();    
  }
  else{
    ERROR_LOG << "There are only three load options (READ, RECREATE, UPDATE) - You have selected " << option << std::endl;
    return;
  }

  INFO_LOG << "Sucessfully attached Disk Resident HyperBinning to file " << filename << std::endl;

}

bool HyperBinningDiskRes::isDiskResident() const{
  return true;
}
TString HyperBinningDiskRes::filename() const{
  if (_file == 0){
    ERROR_LOG << "HyperBinningDiskRes::filename - there is no file associated to this object yet" << std::endl;
    return "";
  }
  return _file->GetName();
}

void HyperBinningDiskRes::reserveCapacity(int nElements){
  HyperBinning::reserveCapacity(nElements);  
  
}



///Destructor
///
HyperBinningDiskRes::~HyperBinningDiskRes(){
  GOODBYE_LOG << "Goodbye from the HyperBinningDiskRes() Constructor";

  if (_file != 0) {
    _file->cd();
    if (_writeable == true){
      _tree       ->Write();
      _treePrimVol->Write();
    }
    _file->Close();
    _file = 0;
  }

  delete _linkedBins;

}



