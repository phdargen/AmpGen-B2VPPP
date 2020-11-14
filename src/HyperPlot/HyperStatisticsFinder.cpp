#include "AmpGen/HyperPlot/HyperStatisticsFinder.h"


HyperStatisticsFinder::HyperStatisticsFinder(int dim, bool mean, bool width, bool widthError, bool keepOrderedEvents) :
  _dim(dim)
{
  
  for (int i = 0; i < _dim; i++){
    std::vector<StatisticsFinder> statsFinders;
    for (int j = 0; j < _dim; j++){
      statsFinders.push_back( StatisticsFinder(mean, width, widthError, keepOrderedEvents) );
    }
    _statisticsFinders.push_back(statsFinders);
  }

}

HyperStatisticsFinder::HyperStatisticsFinder(const HyperPoint& x, bool mean, bool width, bool widthError, bool keepOrderedEvents) :
  _dim(x.getDimension())
{
  
  for (int i = 0; i < _dim; i++){
    std::vector<StatisticsFinder> statsFinders;
    for (int j = 0; j < _dim; j++){
      statsFinders.push_back( StatisticsFinder(mean, width, widthError, keepOrderedEvents) );
    }
    _statisticsFinders.push_back(statsFinders);
  }
  
  add(x);

}

void HyperStatisticsFinder::add( const HyperPoint& x ){

  if (x.getDimension() != _dim) ERROR_LOG << "The HyperPoint you are adding is not of the correct dimension";
  
  for (int i = 0; i < _dim; i++){ 
    for (int j = 0; j < _dim; j++){
      double val = 0.0;
      if (i == j) val = x.at(i);
      if (i != j) val = x.at(i)*x.at(j);
      _statisticsFinders.at(i).at(j).add( val, x.getWeight() );
    }
  }  

}

double HyperStatisticsFinder::correlation(int i, int j) const{
  
  return covarience(i,j)/sqrt(covarience(i,i)*covarience(j,j));

}






double HyperStatisticsFinder::covarience(int i, int j) const{
  
  if (i == j) {
    return _statisticsFinders.at(i).at(i).varience();
  }

  const StatisticsFinder& statsI  = _statisticsFinders.at(i).at(i);
  const StatisticsFinder& statsJ  = _statisticsFinders.at(j).at(j);
  const StatisticsFinder& statsIJ = _statisticsFinders.at(i).at(j);
  
  double expI  = statsI .mean();
  double expJ  = statsJ .mean();
  double expIJ = statsIJ.mean();

  return expIJ - expI*expJ;

}


double HyperStatisticsFinder::mean(int i) const{
  return _statisticsFinders.at(i).at(i).mean();
}

double HyperStatisticsFinder::meanError(int i) const{
  return _statisticsFinders.at(i).at(i).meanError();
}

double HyperStatisticsFinder::width(int i) const{
  return _statisticsFinders.at(i).at(i).width();
}

double HyperStatisticsFinder::getMin(int i) const{
  return _statisticsFinders.at(i).at(i).getMin();
}

double HyperStatisticsFinder::getMax(int i) const{
  return _statisticsFinders.at(i).at(i).getMax();
}




HyperPoint HyperStatisticsFinder::mean() const{
  HyperPoint point(_dim);
  for (int i = 0; i < _dim; i++) point.at(i) = _statisticsFinders.at(i).at(i).mean();
  return point;
}

HyperPoint HyperStatisticsFinder::meanError() const{
  HyperPoint point(_dim);
  for (int i = 0; i < _dim; i++) point.at(i) = _statisticsFinders.at(i).at(i).meanError();
  return point;
}

HyperPoint HyperStatisticsFinder::width() const{
  HyperPoint point(_dim);
  for (int i = 0; i < _dim; i++) point.at(i) = _statisticsFinders.at(i).at(i).width();
  return point;
}

HyperPoint HyperStatisticsFinder::getMin() const{
  HyperPoint point(_dim);
  for (int i = 0; i < _dim; i++) point.at(i) = _statisticsFinders.at(i).at(i).getMin();
  return point;
}

HyperPoint HyperStatisticsFinder::getMax() const{
  HyperPoint point(_dim);
  for (int i = 0; i < _dim; i++) point.at(i) = _statisticsFinders.at(i).at(i).getMax();
  return point;
}


HyperStatisticsFinder::~HyperStatisticsFinder(){

}


