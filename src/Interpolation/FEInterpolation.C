// $Id: FEInterpolation.C,v 1.6 2006-09-05 20:41:48 manav Exp $

// C++ includes

// FESystem includes
#include "Interpolation/FEInterpolation.h"
#include "Interpolation/InterpolationBase.h"
#include "Solvers/DirectLinearSolver.h"
//#include "Interpolation/FEInterpolationAnalysis.h"
#include "Interpolation/InterpolationDriver.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"

// libMesh inclueds



FEInterpolation::FEInterpolation(InterpolationDriver& driver):
InterpolationBase(driver),
LinearAnalysisDriver(1, driver.getFESystemController())
{
  
}



FEInterpolation::~FEInterpolation()
{
  
}



void FEInterpolation::getInterpolatedValues(NumericVector<double>& sol,
                                            NumericVector<double>& interpolated_sol)
{
  // here, the driver will be asked to solve, for which, the driver will ask 
  // the discipline for the matrices and vectors and finally solve 
  // the linear system of equations to get the interpolated value.
  
  // first set the local sol to the given vector
  
  this->sol_to_interpolate = &sol;
  
  // now having set this, ask the linear analysis driver to solve for the 
  // interpolated solution
  this->solveCurrentLoadCase();
  
  // finally, now that the system has been solved, fill the interpolated_sol 
  // with the calculated solution and return
  interpolated_sol = this->getVector("Solution");
}

