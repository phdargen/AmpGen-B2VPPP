#include "AmpGen/HyperPlot/Weights.h"

/// Empty constructor
///
Weights::Weights()
{

}

/// Construct instance with a single weight (probably the most used
/// scenario)
Weights::Weights(double weight) :
  _weights(1, weight)
{

}

/// How many weights have been stored
///
int Weights::numWeights() const{
  return _weights.size();
}

/// print weights to specified ostream (default is cout)
///
void Weights::printWeight(std::ostream& os, int endline) const{
  if (numWeights() != 0){
    os << "weights = ";
    for (int i = 0; i < numWeights(); i++) {
      os << getWeight(i);
      if (i != (numWeights() -1)) os << ", ";
    }
  }
  if(endline == 1) os << std::endl;

}

/// Get weight i. If no weights exist just return 1.0.
///
double Weights::getWeight(int i) const{
  if ( numWeights() == 0 && i == 0) return 1.0; //If no weights, assume it's unweighted i.e. w = 1.0
  if (i >= numWeights()){
    ERROR_LOG << "There are not this many weights avaliable. Returning weight 0 - may crash if this doesn't exist";
    return _weights.at(0);
  }
  return _weights.at(i);
}


/// Set weight i to specified value - if this weight doesn't yet 
/// exist, keep adding weights of 1.0 until it does
void Weights::setWeight(int i, double w){ 

  if (i >= numWeights()) {
  	addWeight(1.0); 
  	setWeight(i, w);
  }
  else{
    _weights.at(i) = w; 
  }
}

/// Set weight i=0 to desired value
///
void Weights::setWeight(double w){ 
  setWeight(0, w);
}

/// Add weight to the end of the weight vector
///
void Weights::addWeight(const double& weight){
  _weights.push_back(weight);
}

Weights::~Weights(){

}

