#include "AmpGen/HyperPlot/HyperName.h"

///Fill HyperName with the standard names {var_0, var_1, ...}
///
HyperName::HyperName(int dim){
  for(int i = 0; i < dim; i++){
    TString name = "var_"; name += i;
    _names.push_back( name );
    _units.push_back( " "  );
  }
}

///Fill HyperName with a std::vector<TString> of names
///
HyperName::HyperName(std::vector<TString> names){
  _names = names;
  _units = std::vector<TString>(names.size(), " ");

}

///Make a 1D HyperName with the given name
///
HyperName::HyperName(TString name0){
  _names.push_back(name0);
  _units.push_back(" ");
}

///Make a 2D HyperName with the given names
///
HyperName::HyperName(TString name0, TString name1){
  _names.push_back(name0);
  _names.push_back(name1);
  _units.push_back(" ");
  _units.push_back(" ");
}

///Make a 3D HyperName with the given names
///
HyperName::HyperName(TString name0, TString name1, TString name2){
  _names.push_back(name0);
  _names.push_back(name1);
  _names.push_back(name2);
  _units.push_back(" ");
  _units.push_back(" ");
  _units.push_back(" ");
}

///Make a 4D HyperName with the given names
///
HyperName::HyperName(TString name0, TString name1, TString name2, TString name3){
  _names.push_back(name0);
  _names.push_back(name1);
  _names.push_back(name2);
  _names.push_back(name3);
  _units.push_back(" ");
  _units.push_back(" ");
  _units.push_back(" ");
  _units.push_back(" ");
}

///Make a 5D HyperName with the given names
///
HyperName::HyperName(TString name0, TString name1, TString name2, TString name3, TString name4){
  _names.push_back(name0);
  _names.push_back(name1);
  _names.push_back(name2);
  _names.push_back(name3);
  _names.push_back(name4);
  _units.push_back(" ");
  _units.push_back(" ");
  _units.push_back(" ");
  _units.push_back(" ");
  _units.push_back(" ");
}

///The the units of every dimension to be the given TString
///
void HyperName::setUnits(TString units){
  for(unsigned i = 0; i < _units.size(); i++){
    _units.at(i) = units;
  }
}

///The the units of the given dimension to be the given TString
///
void HyperName::setUnits(TString units, int dim){
  if (dim < 0 || dim >= (int)_names.size()){
    std::cout << "Trying to set the units of the HyperName element that doesn't exist" << std::endl;
    return;
  }
  _units.at(dim) = units;
}

///Make a TString sutable for an axis title for a 
///given dimension.
TString HyperName::getAxisString(int dim){
  if (dim < 0 || dim >= (int)_names.size()){
    std::cout << "Trying call HyperName::getAxisString on an element that doesn't exist" << std::endl;
    return "ERROR";
  }  
  TString retVal = _names.at(dim);
  if (_units.at(dim) != " ") retVal += " [" + _units.at(dim) + "]";
  return retVal;
}


int HyperName::getDimension() const{
  return _names.size();
}

/**
  Slice (remove) the given dimesnions, and return new HyperName
*/
HyperName HyperName::slice(std::vector<int> sliceDims){
  
  int nStartingDims = getDimension();
  int nSliceDims    = sliceDims.size();

  HyperName names( nStartingDims - nSliceDims );
  int count = 0;

  for (int i = 0; i < nStartingDims; i++){

    bool doesExist = false;

    for (int j = 0; j < nSliceDims; j++){
      int dim = sliceDims.at(j);
      if (i == dim) {
        doesExist = true;
        break;
      }
    }

    if (doesExist == false) { 
      names.at(count) = at(i); count++; 
    }
  }

  return names;    

}



///Destructor
///
HyperName::~HyperName(){

}


