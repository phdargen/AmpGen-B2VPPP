#include "AmpGen/HyperPlot/HyperBinningAlgorithms.h"

///The empty constuctor which is private. This means
///it can only be called from a static member function
AlgOption::AlgOption() :
  _optionName(EMTPY),
  _bool(false),
  _int(0),
  _double(0.0),
  _hyperPoint(0),
  _hyperPointSet(0),
  _hyperName(0),
  _hyperFunc(0)
{
  
}

///Get the EMTPY AlgOption, which takes no arguments.
///  
AlgOption AlgOption::Empty              ( ){
  AlgOption algOption;
  algOption._optionName = EMTPY;
  return algOption;
} 

///Get the START_DIM AlgOption, which tells the binning
///algorithm what dimension to split first (not applicable
/// to all algorithms). 
AlgOption AlgOption::StartDimension     (int dim ){
  AlgOption algOption;
  algOption._optionName = START_DIM;
  algOption._int = dim;  
  return algOption;
} 

///Get the BINNING_DIMS AlgOption, which tells the binning
/// algorithm what dimensions its allowed to split
AlgOption AlgOption::BinningDimensions  (std::vector<int> dims){
  AlgOption algOption;
  algOption._optionName = BINNING_DIMS;
  algOption._intvector = dims;
  return algOption;
} 

///Get the RAND_SEED AlgOption, which tells the binning algorithm
///what random seed to use (not applicable to all algorithms). 
AlgOption AlgOption::RandomSeed         (int seed){
  AlgOption algOption;
  algOption._optionName = RAND_SEED;
  algOption._int = seed;
  return algOption;
} 

///Get the MIN_BIN_WIDTH AlgOption, which tells the binning algorithm
///the minimum bin width that is allowed for ALL dimensions
AlgOption AlgOption::MinBinWidth        (double width){
  AlgOption algOption;
  algOption._optionName = MIN_BIN_WIDTH;
  algOption._double = width;
  return algOption;
} 

///Get the MIN_BIN_WIDTH AlgOption, which tells the binning algorithm
///the minimum bin width that is allowed for EACH dimension
AlgOption AlgOption::MinBinWidth        (HyperPoint widths){
  AlgOption algOption;
  algOption._optionName = MIN_BIN_WIDTH;
  algOption._hyperPoint = widths;
  return algOption;
} 

///Get the MIN_BIN_CONTENT AlgOption, which tells the binning algorithm
///the minimum bin content allowed.
AlgOption AlgOption::MinBinContent      (double val ){
  AlgOption algOption;
  algOption._optionName = MIN_BIN_CONTENT;
  algOption._double = val;
  return algOption;
} 

///Get the MIN_SHADOW_BIN_CONTENT AlgOption, which tells the binning algorithm
///the minimum bin content allowed in the shadow dataset.
AlgOption AlgOption::MinShadowBinContent(double val ){
  AlgOption algOption;
  algOption._optionName = MIN_SHADOW_BIN_CONTENT;
  algOption._double = val;
  return algOption;
} 

///Get the USE_WEIGHTS AlgOption, which tells the binning algorithm
///if it should use the event weights
AlgOption AlgOption::UseWeights         (bool   val){
  AlgOption algOption;
  algOption._optionName = USE_WEIGHTS;
  algOption._bool = val;
  return algOption;
} 

///Get the USE_SHADOW_DATA AlgOption, which tells the binning algorithm
///if it should use the a shadow dataset. This allows each bin to have 
///a specified number of events from both the dataset and the shadow dataset.
AlgOption AlgOption::UseShadowData      (const HyperPointSet& data){
  AlgOption algOption;
  algOption._optionName = USE_SHADOW_DATA;
  algOption._hyperPointSet = &data;
  return algOption;
}


///Use this if you want to draw the HyperBinning after every iteration
///
AlgOption AlgOption::DrawAlgorithm      (TString path){
  AlgOption algOption;
  algOption._optionName = DRAW_ALGORITHM;
  algOption._string = path;
  return algOption;  
}

///Use this if you want to set the axis titles - can also do this later,
///unless you want axis titles for the DrawAlgorithm() option.
AlgOption AlgOption::AxisTitles         (HyperName name){
  AlgOption algOption;
  algOption._optionName = AXIS_NAMES;
  algOption._hyperName = name;
  return algOption;   
}

///Use this if you want to provide a HyperFunction
AlgOption AlgOption::UseFunction          (HyperFunction* func){
  AlgOption algOption;
  algOption._optionName = FUNC;
  algOption._hyperFunc = func;
  return algOption;   
}


AlgOption AlgOption::GridMultiplier     (int val){
  AlgOption algOption;
  algOption._optionName = GRID_MULTIPLIER;
  algOption._int = val;
  return algOption;   
}

AlgOption AlgOption::GridMultiplier     (HyperPoint val){
  AlgOption algOption;
  algOption._optionName = GRID_MULTIPLIER;
  algOption._hyperPoint = val;
  return algOption;   
}

AlgOption AlgOption::SnapToGrid         (bool val){
  AlgOption algOption;
  algOption._optionName = SNAP_TO_GRID;
  algOption._bool = val;
  return algOption;     
}


AlgOption AlgOption::NumPhaseBinPairs   (int val){
  AlgOption algOption;
  algOption._optionName = NUM_BIN_PAIRS;
  algOption._int = val;
  return algOption;      
}

AlgOption AlgOption::PhaseBinEdges   (std::vector<double> val){
  AlgOption algOption;
  algOption._optionName = PHASE_BIN_EDGES;
  algOption._doublevector = val;
  return algOption;      
}

AlgOption AlgOption::StartBinning       (const HyperBinning& val){
  AlgOption algOption;
  algOption._optionName = START_BINNING;
  algOption._hyperBinning = &val;
  return algOption;      
}



///Get the AlgOption::OptionName 
///
AlgOption::OptionName  AlgOption::getOptionName            (){
  return _optionName;
}

///Get the boolean option
///
bool                 AlgOption::getBoolOpt         (){
  return _bool;
}

///Get the integer option
///
int                  AlgOption::getIntOpt          (){
  return _int;
}

///Get the double member
///
double               AlgOption::getDoubleOpt       (){
  return _double;
}

///Get the std::vector<int>  member
///
std::vector<int>     AlgOption::getIntVectorOpt    (){
  return _intvector;
}

///Get the string  member
///
TString              AlgOption::getStringOpt       (){
  return _string;
}


///Get the std::vector<double> member
///
std::vector<double>  AlgOption::getDoubleVectorOpt (){
  return _doublevector;
}

///Get the HyperPoint member
///
HyperPoint           AlgOption::getHyperPointOpt   (){
  return _hyperPoint;
}

///Get the HyperPointSet member
///
const HyperPointSet&    AlgOption::getHyperPointSetOpt(){
  if (_hyperPointSet == 0){
    std::cout << "ERROR: No HyperPointSet hass been passed to AlgOption!" << std::endl;
    return *(new HyperPointSet(0));
  }
  return *_hyperPointSet;
}

HyperFunction*       AlgOption::getFuncOpt         (){
  return _hyperFunc;
}


///Get the HyperName member
///
HyperName           AlgOption::getHyperNameOpt   (){
  return _hyperName;
}

///Get the HyperBinning member
///
const HyperBinning*  AlgOption::getHyperBinningOpt (){
  return _hyperBinning; 
}


///Check if the OptionName is EMTPY
/// 
bool AlgOption::isEmpty(){
  return _optionName == EMTPY;
}

//Desctructor
AlgOption::~AlgOption(){

}


///look through the list of AlgOption's and see if
///one with a specific OptionName exsits. If so, return it,
///if not return the OptionName EMTPY
AlgOption HyperBinningAlgorithms::getOpt(AlgOption::OptionName name){

  for (unsigned i = 0; i < _algOptions.size(); i++) {
    if (_algOptions.at(i).getOptionName() == name) return _algOptions.at(i);
  }

  return AlgOption::Empty(); 
}

///look through the list of AlgOption's and see if
///one with a specific OptionName exsits. 
bool HyperBinningAlgorithms::optExist(AlgOption::OptionName name){

  return getOpt(name).isEmpty() == 0;

}

///Constuct a HyperBinningAlgorithms by selecting the algortihm
///you would like to use
HyperBinningAlgorithms::HyperBinningAlgorithms(Alg algorithm) :
  _alg(algorithm)
{
    
}

///Add an AlgOption which is passed to the binning algortihm
void HyperBinningAlgorithms::addAlgOption(AlgOption option){
  _algOptions.push_back(option);
}

///Get the HyperBinningMaker (the binning algorithm) with the
///chosen algorithm type, and with the chosen AlgOption's
HyperBinningMaker* HyperBinningAlgorithms::getHyperBinningMaker(HyperCuboid binningRange, HyperPointSet points){
     
  int dim = binningRange.getDimension();

  HyperBinningMaker* binnningMaker = 0;
  
  int startDimension = 0;

  if (optExist(AlgOption::START_DIM    )) 
    startDimension = getOpt(AlgOption::START_DIM  ).getIntOpt();  
  
  int numPhaseBinPairs = 3;
  if (optExist(AlgOption::NUM_BIN_PAIRS    )) 
    numPhaseBinPairs = getOpt(AlgOption::NUM_BIN_PAIRS  ).getIntOpt();



  //Chose the binning algorithm

  if (_alg == SMART ) {
    binnningMaker = new HyperBinningMakerSmart(binningRange, points, startDimension);
  }
  if (_alg == MINT ) {
    binnningMaker = new HyperBinningMakerMint(binningRange, points, startDimension);
  }
  if (_alg == MINT_SMART ) {
    binnningMaker = new HyperBinningMakerMintSmart(binningRange, points, startDimension);
  }
  if (_alg == MINT_RANDOM ) {
    binnningMaker = new HyperBinningMakerMintRandomise(binningRange, points);
  }
  if (_alg == SMART_RANDOM ) {
    binnningMaker = new HyperBinningMakerSmartRandomise(binningRange, points);
  }
  if (_alg == LIKELIHOOD ) {
    binnningMaker = new HyperBinningMakerLikelihood(binningRange, points);
  }
  if (_alg == SMART_LIKELIHOOD ) {
    binnningMaker = new HyperBinningMakerSmartLikelihood(binningRange, points);
  }  
  if (_alg == SMART_MULTI) {
    binnningMaker = new HyperBinningMakerMultiSmart(binningRange, points, startDimension);
  }
  if (_alg == FUNC_PHASE) {
    HyperBinningMakerPhaseBinning* phaseBinningMaker = new HyperBinningMakerPhaseBinning(binningRange, 0);
    phaseBinningMaker->setNumBinPairs(numPhaseBinPairs);

    if (optExist(AlgOption::PHASE_BIN_EDGES    )) 
      phaseBinningMaker->setBinEdges( getOpt(AlgOption::PHASE_BIN_EDGES  ).getDoubleVectorOpt() );

    binnningMaker = dynamic_cast<HyperBinningMaker*>(phaseBinningMaker);
  }
  //Now set the options


  if (optExist(AlgOption::RAND_SEED    )) 
    binnningMaker->setSeed  ( getOpt(AlgOption::RAND_SEED  ).getIntOpt() );

  if (optExist(AlgOption::USE_WEIGHTS  )) 
    binnningMaker->useEventWeights( getOpt(AlgOption::USE_WEIGHTS).getBoolOpt() );

  if (optExist(AlgOption::BINNING_DIMS )) 
    binnningMaker->setBinningDimensions( getOpt(AlgOption::BINNING_DIMS).getIntVectorOpt() );

  if (optExist(AlgOption::FUNC )) 
    binnningMaker->setHyperFunction( getOpt(AlgOption::FUNC).getFuncOpt() );

  if (optExist(AlgOption::MIN_BIN_WIDTH )){
    //can provide either a double or a HyperPoint here.

    HyperPoint widthsA = getOpt(AlgOption::MIN_BIN_WIDTH).getHyperPointOpt();
    double     widthsB = getOpt(AlgOption::MIN_BIN_WIDTH).getDoubleOpt    (); 
    if ( widthsA.getDimension() == 0 ) widthsA = HyperPoint(dim, widthsB);
    binnningMaker->setMinimumEdgeLength(widthsA);
  }

  if (optExist(AlgOption::USE_SHADOW_DATA )) 
    binnningMaker->addShadowHyperPointSet( getOpt(AlgOption::USE_SHADOW_DATA).getHyperPointSetOpt() );

  if (optExist(AlgOption::MIN_BIN_CONTENT )) 
    binnningMaker->setMinimumBinContent( getOpt(AlgOption::MIN_BIN_CONTENT).getDoubleOpt() );

  if (optExist(AlgOption::MIN_SHADOW_BIN_CONTENT )) 
    binnningMaker->setShadowMinimumBinContent( getOpt(AlgOption::MIN_BIN_CONTENT).getDoubleOpt() );

  if (optExist(AlgOption::DRAW_ALGORITHM )) 
    binnningMaker->drawAfterEachIteration( getOpt(AlgOption::DRAW_ALGORITHM).getStringOpt() );

  if (optExist(AlgOption::AXIS_NAMES )) 
    binnningMaker->setNames( getOpt(AlgOption::AXIS_NAMES).getHyperNameOpt() );
  
  if (optExist(AlgOption::SNAP_TO_GRID )){
    binnningMaker->useSnapToGrid( getOpt(AlgOption::SNAP_TO_GRID).getBoolOpt() );
  }

  if (optExist(AlgOption::GRID_MULTIPLIER )){

    HyperPoint valA = getOpt(AlgOption::GRID_MULTIPLIER).getHyperPointOpt();    
    int        valB = getOpt(AlgOption::GRID_MULTIPLIER).getIntOpt       ();

    if (valA.getDimension() != 0){
      binnningMaker->setGridMultiplier(valA);
    }
    else{
      binnningMaker->setGridMultiplier(valB);
    }
  }

  if (optExist(AlgOption::START_BINNING) ){
    binnningMaker->updateFromExistingHyperBinning( *getOpt(AlgOption::START_BINNING).getHyperBinningOpt() );
  }
  
  return binnningMaker;

}

///Destructor
///
HyperBinningAlgorithms::~HyperBinningAlgorithms(){

}






