// $Id: FEInterpolation.h,v 1.6 2006-09-05 20:41:48 manav Exp $

#ifndef __fesystem_fe_interpolation_h__
#define __fesystem_fe_interpolation_h__


// C++ incliudes


// FESystem includes
#include "Interpolation/InterpolationBase.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"


class InterpolatinDriver;

// Forward declerations
class Elem;


/// A FE based interpolation class that derives from the InterpolationBase 
/// class. The basic purpose of this class is to create the interpolation
/// solution vector. These values are created based on the solution of the
/// system
class FEInterpolation: 
public InterpolationBase, public Driver::LinearAnalysisDriver 
{
 public:
  /// constructor
  FEInterpolation(InterpolationDriver& driver);

  /// destructor 
  ~FEInterpolation();

  /// returns the interpolated for the interpolated values for the given 
  /// global solution vector
  /// @param sol solution vector which needs to be interpolated
  /// @param interpolated_sol vector in which the interpolated solution will
  /// be returned
  virtual void getInterpolatedValues(NumericVector<double>& sol,
                                     NumericVector<double>& interpolated_sol);
};



#endif // __fesystem_fe_interpolation_h__
