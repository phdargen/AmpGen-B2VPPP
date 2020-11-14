/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Oct 2016
 *
 * A vector of HyperCuboids. This is identical to a HyperVolume, but will be used
 * for different situations. For instance a hyperVolume would be used to describe some
 * multi-dimensional volume e.g. a bin in a HyperVolumeHistogram. A HyperCuboidSet is
 * simply a container for a load of HyperCuboids. 
 *
 **/

#ifndef HYPERCUBOIDSET_HH
#define HYPERCUBOIDSET_HH

#include "AmpGen/HyperPlot/HyperVolume.h"

typedef HyperVolume HyperCuboidSet;


#endif


