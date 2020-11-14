/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Jan 2017
 *
 * This class is used to find statistics (e.g. mean, width...) of a multi dimensional (Hyper) dataset. 
 * To do so it uses a 2D array of the StatisticsFinder class. I suspect there is a computationally
 * less expensive way to do this, but this does the job for now.
 *
 **/

 
#ifndef HYPER_STATISTICS_FINDER_HH
#define HYPER_STATISTICS_FINDER_HH

// HyperPlot includes

#include "AmpGen/HyperPlot/MessageService.h"
#include "AmpGen/HyperPlot/StatisticsFinder.h"
#include "AmpGen/HyperPlot/HyperPoint.h"

class HyperStatisticsFinder {

  int _dim; 
  /**< dimension of the HyperStatisticsFinder */

  std::vector< std::vector < StatisticsFinder > > _statisticsFinders; 
  /**< matrix of StatisticsFinder's. Need a matrix rather than a vector in order to
  calculate the covarience between two dimensions. */
  
  public:

  HyperStatisticsFinder( int dimension, bool mean = 1, bool width = 1, bool widthError = 1, bool keepOrderedEvents = 0);
  /**<
  In the constuctor you decide what things you want to be stored once you start
  adding values. This will determine what statistics you are able to calcuate later.
  Also decide the dimensionality of each data point.
  */  

  HyperStatisticsFinder( const HyperPoint& point, bool mean = 1, bool width = 1, bool widthError = 1, bool keepOrderedEvents = 0);
  /**<
  In the constuctor you decide what things you want to be stored once you start
  adding values. This will determine what statistics you are able to calcuate later.
  Also decide the dimensionality of each data point by passing a HyperPoint. This
  point is also added to the HyperStatisticsFinder
  */  

  void  add( const HyperPoint& x );
  /**< add a HyperPoint to the HyperStatisticsFinder */
  double correlation(int i, int j) const;
  /**< get the correlation coefficient of dimesnions i and j */
  double covarience(int i, int j) const;
  /**< get the covarience of dimesnions i and j */

  double mean     (int i) const;
  /**< get the mean in a given dimension */
  double meanError(int i) const;
  /**< get the meanError in a given dimension */
  double width    (int i) const;
  /**< get the width in a given dimension */
  double getMin   (int i) const;
  /**< get the minimum in a given dimension */
  double getMax   (int i) const;
  /**< get the maximum in a given dimension */

  HyperPoint mean     () const;
  /**< get a HyperPoint filled with the mean of each dimension*/
  HyperPoint meanError() const;
  /**< get a HyperPoint filled with the error on the mean of each dimension*/
  HyperPoint width    () const;
  /**< get a HyperPoint filled with the width of each dimension*/
  HyperPoint getMin   () const;
  /**< get a HyperPoint filled with the minimum of each dimension*/
  HyperPoint getMax   () const;
  /**< get a HyperPoint filled with the maximum of each dimension*/

  const StatisticsFinder& getStatisticsFinder(int i) const{return _statisticsFinders.at(i).at(i);}
  /**< get the statistics finder for a specific dimesion */

  virtual ~HyperStatisticsFinder();
  /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Jan 2017
 *
 * Make a HyperStatisticsFinder, but tell it to only store 
 * enough information to calcuate the minimum and maximum
 * 
 **/
class HyperMinMaxFinder : public HyperStatisticsFinder{
  public:
  HyperMinMaxFinder(int dim) : HyperStatisticsFinder(dim,0,0,0,0){} /**< Construct with given dimension */
  HyperMinMaxFinder(const HyperPoint& point) : HyperStatisticsFinder(point,0,0,0,0){}  
  /**< Construct with dimension of given HyperPoint, and add this point to the statistics finder */
  ~HyperMinMaxFinder(){}

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Jan 2017
 *
 * Make a HyperStatisticsFinder, but tell it to only store 
 * enough information to calcuate the mean
 * 
 **/
class HyperMeanFinder : public HyperStatisticsFinder{
  public:
  HyperMeanFinder(int dim) : HyperStatisticsFinder(dim,1,0,0,0){} /**< Construct with given dimension */
  HyperMeanFinder(const HyperPoint& point) : HyperStatisticsFinder(point,1,0,0,0){} 
  /**< Construct with dimension of given HyperPoint, and add this point to the statistics finder */
  ~HyperMeanFinder(){} /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Jan 2017
 *
 * Make a HyperStatisticsFinder, but tell it to only store 
 * enough information to calcuate the mean and rms
 * 
 **/
class HyperWidthFinder : public HyperStatisticsFinder{
  public:
  HyperWidthFinder(int dim) : HyperStatisticsFinder(dim,1,1,0,0){}  /**< Construct with given dimension */
  HyperWidthFinder(const HyperPoint& point) : HyperStatisticsFinder(point,1,1,0,0){}  
  /**< Construct with dimension of given HyperPoint, and add this point to the statistics finder */
  ~HyperWidthFinder(){} /**< Destructor */

};

/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Jan 2017
 *
 * Make a HyperStatisticsFinder, but tell it to only store 
 * enough information to calcuate the mean and rms, and error on the rms
 * 
 **/
class HyperWidthErrorFinder : public HyperStatisticsFinder{
  public:
  HyperWidthErrorFinder(int dim) : HyperStatisticsFinder(dim,1,1,1,0){} /**< Construct with given dimension */
  HyperWidthErrorFinder(const HyperPoint& point) : HyperStatisticsFinder(point,1,1,1,0){} 
  /**< Construct with dimension of given HyperPoint, and add this point to the statistics finder */
  ~HyperWidthErrorFinder(){} /**< Destructor */

};

#endif
