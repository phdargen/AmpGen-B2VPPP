/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Used to find statistics for a dataset. 
 * 
 **/
 
#ifndef STATISTICS_FINDER_HH
#define STATISTICS_FINDER_HH

#include "AmpGen/HyperPlot/MessageService.h"

class StatisticsFinder {

  double _min;     /**< The smallest member added to the StatisticsFinder */
  double _max;     /**< The largest  member added to the StatisticsFinder */
  
  double _nEvents; /**< The number of events added to the StatisticsFinder */
  double _wSum;    /**< The weighted sum of values   added to the StatisticsFinder */
  double _wSum2;   /**< The weighted sum of values^2 added to the StatisticsFinder */
  double _wSum3;   /**< The weighted sum of values^3 added to the StatisticsFinder */
  double _wSum4;   /**< The weighted sum of values^4 added to the StatisticsFinder */
  double _sumW;    /**< The sum of weights */

  //Can specify to keep track of the ordered events
  //allowing things like the median to be calculated

  int _keepOrderedEvents; /**<  Keep a list of the values added (so the median can be found) */
  mutable std::vector<double> _orderedEvents; /**<  list of the values added  */

  bool needOrderedEvents() const;
  void warnIfWeightedEvents() const;
  bool notEnoughInformation(const double& val) const;

  public:

  StatisticsFinder(bool mean = 1, bool width = 1, bool widthError = 1, bool keepOrderedEvents = 0);

  double median() const;
  
  double numEvents() const{return _nEvents;}  /**<  numer of events added to the StatisticsFinder  */

  void add(const double& x, const double& weight = 1.0);

  double mean() const;
  double meanError() const;
  double varience() const;
  double width() const;
  double widthError() const;
 
  double expX() const; 
  double expX2() const;
  double expX3() const;
  double expX4() const;
  
  double secondCentralMom() const;
  double fourthCentralMom() const;

  const double& getMin() const{return _min;}  /**<  min value added to the StatisticsFinder  */
  const double& getMax() const{return _max;}  /**<  max value added to the StatisticsFinder  */
  
  double range() const{return _max - _min;} /**<  max - min value added to the StatisticsFinder  */

  virtual ~StatisticsFinder();

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Make a StatisticsFinder, but tell it to only store 
 * enough information to calcuate the minimum and maximum
 * 
 **/
class MinMaxFinder : public StatisticsFinder{
  public:
  MinMaxFinder() : StatisticsFinder(0,0,0,0){} /**< Constructor */
  ~MinMaxFinder(){} /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Make a StatisticsFinder, but tell it to only store 
 * enough information to calcuate the mean
 * 
 **/
class MeanFinder : public StatisticsFinder{
  public:
  MeanFinder() : StatisticsFinder(1,0,0,0){} /**< Constructor */
  ~MeanFinder(){} /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Make a StatisticsFinder, but tell it to only store 
 * enough information to calcuate the mean and width
 * (error on the mean also comes for free)
 **/
class WidthFinder : public StatisticsFinder{
  public:
  WidthFinder() : StatisticsFinder(1,1,0,0){} /**< Constructor */
  ~WidthFinder(){} /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Make a StatisticsFinder, but tell it to only store 
 * enough information to calcuate the mean, width, 
 * and error on the width 
 **/
class WidthErrorFinder : public StatisticsFinder{
  public:
  WidthErrorFinder() : StatisticsFinder(1,1,1,0){} /**< Constructor */
  ~WidthErrorFinder(){} /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Make a StatisticsFinder, and let it store enough info
 * to calculate the median 
 **/
class MedianFinder : public StatisticsFinder{
  public:
  MedianFinder() : StatisticsFinder(1,1,1,1){} /**< Constructor */
  ~MedianFinder(){} /**< Destructor */

};

#endif
