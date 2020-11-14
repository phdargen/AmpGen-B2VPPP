/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Base class for all RootHistogram plotters - these just
 * make plotting root histograms a little easier.
 *
 **/

 
#ifndef PLOTTER_HH
#define PLOTTER_HH

#include "TH1.h"
#include "TObject.h"
#include "TLegend.h"
#include "TBox.h"
#include "TArc.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"
#include "AmpGen/HyperPlot/MessageService.h"


class Plotter {
  
  protected:


  //This class has ownership of the _objToPlot, but not the _histograms

  TPad*            _canvas;           /**< The canvas which add histograms and TObjects are plotted on */
  TLegend*         _legend;           /**< The legend (not really implemented yet) */
  double           _forcedMax;        /**< Force a maximum for plotting histograms */
  double           _forcedMin;        /**< Force a minimum for plotting histogram */
  std::vector<TObject*> _objToPlot;   /**< Vector of TObject%s to plot (e.g TLine or TMarker)  */
  std::vector<TObject*> _histograms;  /**< Vector of Histograms%s to plot (e.g TH1D or TH2D) */
  std::vector<int>      _colours;     /**< Vector of colours to draw the histograms */

  double _lMargin;  /**< left margin size (as a fracition of the pad width)  */
  double _rMargin;  /**< right margin size (as a fracition of the pad width)  */
  double _tMargin;  /**< top margin size (as a fracition of the pad width)  */
  double _bMargin;  /**< bottom margin size (as a fracition of the pad width)  */

  double _xAxisTitleOffset; /**< x-axis title offset  */
  double _yAxisTitleOffset; /**< y-axis title offset  */
  
  double _xAxisLabelOffset;  /**< x-axis label offset  */
  double _yAxisLabelOffset;  /**< y-axis label offset  */

  double _xAxisTickLength; /**< x-axis tick length  */
  double _yAxisTickLength; /**< y-axis tick length  */

  double _xAxisLabelSize;  /**< x-axis label size  */
  double _yAxisLabelSize;  /**< y-axis label size  */

  double _xAxisTitleSize;  /**< x-axis title size  */
  double _yAxisTitleSize;  /**< y-axis title size  */

  bool _histogramOwnership; /**< do I own the Histogram pointers given to me? */
  bool _objectOwnership;    /**< do I own the Objects pointers given to me? */
  
  bool _usePresetColours; /**< Do I use my colour scheme for histograms (true) or the original colours given (false) */
  bool _allImageFormats ; /**< If true, I will save the canvas in all formats when plot() is called */

  virtual void setCanvasDefaults(TPad* pad); 

  public:
    
  static TString s_imageformat;     /**< Primary   image format */
  static TString s_imageformat2;    /**< Secondary image format */
  static TString s_legend_position; /**< Legend positon */
  static int     s_plotterCount;    /**< Plotter count (good for making unique names) */
  static double  s_forcedMax;       /**< Force a maximum histogram value for all plotters */
  static double  s_forcedMin;       /**< Force a minimum histogram value for all plotters */

  int getColor(int i);  /**<  Get histogram colour i (these are preset)  */
  
  void setColor(int i, int color);  /**<  Set histogram colour i  */
  
  void usePresetColours(bool val= true){_usePresetColours = val;}  /**< Should I use the preset colours? */
  void allImageFormats (bool val= true){_allImageFormats  = val;}  
  /**< If switched on, the plotter will save the canvas to .pdf .png .eps and .C formats (as required by CDS) */
  
  void setHistogramOwnership(bool i=1){_histogramOwnership = i;}
  /**< Should I take ownership of the histograms (default is no) */
  void setObjectOwnership(bool i=1){_objectOwnership = i;}
  /**< Should I take ownership of the objects (default is yes) */

  Plotter(TString canvasName, double width, double height);

  Plotter(const Plotter& other);
 
  ///Draw the histograms and objects onto the canvas
  virtual void plot(TString plotDirectory, TString plotOptions = "", TPad* pad = 0, double scaleFactor = 1.0) = 0;
  void add(TObject* histogram);
  void addDot(double xpos, double ypos, double size, int colour = 1, TString shape = "circle", double sizeY = 0.0);
  void logX(bool log = 1);
  void logY(bool log = 1);
  void logZ(bool log = 1);

  //note that any object passed here will have ownership taken by the class
  void addObject(TObject* obj);
  
  int getNumObjects();

  void setImageFormat(TString format);
  

  TPad* getCanvas() { return _canvas; }
  /**< Get the canvas */

  void scaleTextSize(double scale){
    _xAxisLabelSize   *= scale;
    _yAxisLabelSize   *= scale;
    _xAxisTitleSize *= scale;
    _yAxisTitleSize *= scale;
  }
  /**< Scale the labels and titles on the axis by a constant factor */

  void scaleAxisTitleSize(double scale){
    _xAxisTitleSize *= scale;
    _yAxisTitleSize *= scale;
  }
  /**< Scale the titles on the axis by a constant factor */

  void scaleAxisTitleOffset(double scale){
    _xAxisTitleOffset *= scale;
    _yAxisTitleOffset *= scale;
  }
  /**< Scale the title offsets on the axis by a constant factor */

  void setXAxisLabelSize(double val ){_xAxisLabelSize = val;}          /**< set the x-axis label size */
  void setYAxisLabelSize(double val ){_yAxisLabelSize = val;}      /**< set the y-axis label size */
  void setXAxisTitleSize(double val ){_xAxisTitleOffset = val;}      /**< set the x-axis title size */
  void setYAxisTitleSize(double val ){_yAxisTitleOffset = val;}      /**< set the y-axis title size */
  void setXAxisLabelOffset(double val ){_xAxisLabelOffset = val;}      /**< set the x-axis label offset */
  void setYAxisLabelOffset(double val ){_yAxisLabelOffset = val;}      /**< set the y-axis label offset */
  void setXAxisTitleOffset(double val ){_xAxisTitleOffset = val;}      /**< set the x-axis title offset */
  void setYAxisTitleOffset(double val ){_yAxisTitleOffset = val;}      /**< set the y-axis title offset */
  void setXAxisTickLength(double val ){_xAxisTickLength = val;}      /**< set the x-axis tick length */
  void setYAxisTickLength(double val ){_yAxisTickLength = val;}      /**< set the y-axis tick length */
  
  void setPropertiesFromTH1(TH1* hist);
  /**< set _xAxisTitleOffset, _yAxisTitleOffset, _xAxisLabelOffset, 
  _yAxisLabelOffset, _xAxisTickLength , _yAxisTickLength , 
  _xAxisLabelSize  , _yAxisLabelSize  , _xAxisTitleSize  , 
  _yAxisTitleSize  from the histogram given.
  */


  void setMin(double min) { _forcedMin = min; } 
  /**< Set the forced minimum (histograms will be forced to draw with this min)*/
  void setMax(double max) { _forcedMax = max; } 
  /**< Set the forced maximum (histograms will be forced to draw with this max)*/

  void setBMargin(double val){_bMargin = val;} /**< set the left margin size (as a fracition of the pad width)  */
  void setLMargin(double val){_lMargin = val;} /**< set the right margin size (as a fracition of the pad width)  */
  void setRMargin(double val){_rMargin = val;} /**< set the top margin size (as a fracition of the pad width)  */
  void setTMargin(double val){_tMargin = val;} /**< set the bottom margin size (as a fracition of the pad width)  */

  virtual ~Plotter();
  
};

#endif
