/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 *  Holds a vector of weights that is associtated
 *  to a HyperPoint via inheritance (sometimes useful 
 *  to have more than one for systematic studies)
 **/


#ifndef WEIGHTS_HH
#define WEIGHTS_HH

// HyperPlot includes
#include "AmpGen/HyperPlot/MessageService.h"

// Root includes

// std includes


class Weights {

  protected:

  std::vector<double> _weights; /**< vector of weights */

  public:

  Weights();
  Weights(double weight);

  int numWeights() const;
  
  const std::vector<double>& getWeightsVector(){return _weights;} /**< return the vector of weights */

  void printWeight(std::ostream& os = std::cout, int endline=1) const;

  double getWeight(int i = 0) const;

  void setWeight(int i, double w);
  void setWeight(double w);

  void addWeight(const double& weight);

  virtual ~Weights();

};



#endif

