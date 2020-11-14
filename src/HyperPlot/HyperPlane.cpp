#include "AmpGen/HyperPlot/HyperPlane.h"

///Constructor that requires nDim points to define the HyperPlane
///
///
HyperPlane::HyperPlane(const HyperPointSet& a) : 
  NPlane(a)
{ 
  WELCOME_LOG << "Hello from the HyperPlane() Constructor";

  if ( a.getDimension() < (int)a.size() ) ERROR_LOG << "You only need #dim points to define a hyperplane - I'll just use the first #dim";
  if ( a.getDimension() > (int)a.size() ) ERROR_LOG << "You have not provided enough points to define a hyperplane - I'll probably crash soon";

}

///Print the parameters that define the HyperPlane
///
///
void HyperPlane::print(std::ostream& os, int endline) const{

  os << "HyperPlane: v = ";
  _origin.print(os, 0);

  for (unsigned int i = 0; i < _v.size(); i++){
    os << " + ";
    _v.at(i).print(os, 0);
    os << "t_" << i;
  }
  if (endline) os << std::endl;

}

///Check to see if a point is within a HyperPlane. Distance between the 
///two has to be smaller than sqrt(1e-2) 
///
bool HyperPlane::pointInPlane(const HyperPoint& hyperPoint) const{
  
 
  HyperPoint interceptionParameters = findPlaneIntersectionParameters(hyperPoint);

  if (interceptionParameters.size() == 0) return false;
  HyperPoint pointInPlane           = getParametricPoint(interceptionParameters);
  
  HyperPoint vectorDifference = hyperPoint - pointInPlane;

  double difference = 0.0;
  for (int i = 0; i < getDimension(); i++) difference += fabs(vectorDifference.at(i));

  if( fabs(difference) < 1e-2 ) return true;
  
  return false;

}


HyperPoint HyperPlane::findPlaneIntersectionParameters(const HyperPoint& hyperPoint, int dimToOmmit) const{

  if (hyperPoint.getDimension() != getDimension()){
    ERROR_LOG << "Your HyperPlane and HyperPoint have different dimensionality. Returning 0";
    HyperPoint(0);
  }
  
  int oldErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 5000;

  int dim = getDimension();
  int reducedDim = dim -1;

  std::vector<int> linIndep;
  for (int i = 0; i < dim; i++){
    if (dimToOmmit != i) linIndep.push_back(i);
  }

  HyperPoint lhs = hyperPoint - getOrigin();
 
  TMatrixD vectorD(reducedDim, 1);
  for (int i = 0; i < reducedDim ; i++){
    vectorD(i,0) = lhs.at(linIndep.at(i));
  }

  TMatrixD matrix(reducedDim, reducedDim);
  for (int j = 0; j < reducedDim; j++){
    for (int i = 0; i < reducedDim; i++){
      matrix(i,j) = getDirection(j).at(linIndep.at(i));
    }
  }
  
  TDecompLU lu(matrix);
  
  //std::cout << "getting eigenVectors" << std::endl;
  //TVectorT<double> eigenValues;
  //matrix.EigenVectors(eigenValues);
  //std::cout << "done" << std::endl;



  if( !lu.Decompose() ) {
    
    if (dimToOmmit != (dim -1)) return findPlaneIntersectionParameters(hyperPoint, dimToOmmit + 1);
    VERBOSE_LOG << "Matrix is singular. Cannot find solution.";
    matrix.Print();
    
    gErrorIgnoreLevel = oldErrorLevel;

    return HyperPoint(0);
  }

  VERBOSE_LOG << "Matrix before inversion";

  double determinant = 0.0;
  matrix.Invert(&determinant);
  VERBOSE_LOG << "Determinant = " << determinant;

  VERBOSE_LOG << "Matrix after inversion";

  TMatrixD result = matrix*vectorD;
  
  HyperPoint interceptionParameter(reducedDim);
  for (int i = 0; i < reducedDim; i++){
    interceptionParameter.at(i) = result(i,0);
  }

  VERBOSE_LOG << "All done - returning parameters";
   
  gErrorIgnoreLevel = oldErrorLevel;

  return interceptionParameter;

}

HyperPoint HyperPlane::findIntersectionPoint(const HyperLine& hyperLine) const{
  
  HyperPoint lineIntersectionParameter = findLineIntersectionParameter(hyperLine);

  if (lineIntersectionParameter.getDimension()==0) return lineIntersectionParameter;

  return hyperLine.getParametricPoint(lineIntersectionParameter.at(0));

}

HyperPoint HyperPlane::findLineIntersectionParameter(const HyperLine& hyperLine) const{

  TMatrixD solution = findIntersectionSolution(hyperLine);
  if (solution.GetNcols() == 0) return HyperPoint(0);
  return HyperPoint(solution(0,0));

}

HyperPoint HyperPlane::findPlaneIntersectionParameters(const HyperLine& hyperLine) const{

  TMatrixD solution = findIntersectionSolution(hyperLine);
  if (solution.GetNcols() == 0) return HyperPoint(0);

  HyperPoint point(getDimension() - 1);
  
  for (int i = 1; i < getDimension(); i++){
    point.at(i-1) = solution(i,0);
  }  

  return point;
  
}

TMatrixD HyperPlane::findIntersectionSolution(const HyperLine& hyperLine) const{

  if (hyperLine.getDimension() != getDimension()){
    ERROR_LOG << "Your HyperPlane and HyperLine have different dimensionality. Returning 0";
    return TMatrixD(0,0);
  }

  int oldErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 5000;

  HyperPoint lhs = hyperLine.getOrigin() - getOrigin();

  VERBOSE_LOG << "Making Column Vector";  
  TMatrixD vectorD(getDimension(), 1);

  for (int i = 0; i < getDimension(); i++){
    vectorD(i,0) = lhs.at(i);
  }

  VERBOSE_LOG << "Making Matrix";  
  TMatrixD matrix(getDimension(), getDimension());
  
  VERBOSE_LOG << "Filling Column 0"; 
  for (int i = 0; i < getDimension(); i++){
    matrix(i,0) = -hyperLine.getDirection().at(i);
  }

  for (int j = 1; j < getDimension(); j++){
    VERBOSE_LOG << "Filling Column " << j; 
    for (int i = 0; i < getDimension(); i++){
      matrix(i,j) = getDirection(j - 1).at(i);
    }
  }
  
  TDecompLU lu(matrix);
  if( !lu.Decompose() ) {
    VERBOSE_LOG << "Matrix is singular. Cannot find solution.";
    gErrorIgnoreLevel = oldErrorLevel;
    return TMatrixD(0,0);
  }

  VERBOSE_LOG << "Matrix before inversion";
  //if (_verbose) matrix.Print();

  double determinant = 0.0;
  matrix.Invert(&determinant);
  VERBOSE_LOG << "Determinant = " << determinant;

  VERBOSE_LOG << "Matrix after inversion";
  //if (_verbose) matrix.Print();

  TMatrixD result = matrix*vectorD;

  gErrorIgnoreLevel = oldErrorLevel;

  return result;

}

///See if two planes are the same. Cannot just comapre the parameters.
///Instead see if the points that define plane A are within plane B.
///
bool HyperPlane::operator ==(const HyperPlane& other) const{
  if (other.getDimension() != this->getDimension()) {
    ERROR_LOG << "Trying to compare HyperPlanes of different dimensions, returning false";
    return false; 
  }
  if ( this->pointInPlane(other.getOrigin()) == 0) return false;

  for (int i = 0; i < getDimension() - 1; i++){
    HyperPoint parameters(getDimension() - 1);
    parameters.at(i) = 1;
    HyperPoint pointInOtherPlane = other.getParametricPoint(parameters);
    if ( this->pointInPlane(pointInOtherPlane) == 0) return false;
  }
  return true;
}


HyperPlane::~HyperPlane() { 
  GOODBYE_LOG << "Goodbye from the HyperPlane() Constructor";
}



