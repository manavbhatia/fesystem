// $Id: NonLinearAnalysisDriver.h,v 1.11.6.2 2007/05/08 05:18:27 manav Exp $

#ifndef __fesystem_non_linear_analysis_driver_h__
#define __fesystem_non_linear_analysis_driver_h__

// C++ includes


// FESystem includes
#include "AnalysisDriver/AnalysisDriver.h"


// forward declerations
namespace FESystem
{
  class FESystemController;
}


namespace Discipline
{
  class DisciplineInfo;
}


#ifndef NONLINEAR_ANALYSIS_DRIVER_ENUM_ID
#define NONLINEAR_ANALYSIS_DRIVER_ENUM_ID 2
#else
#error
#endif

#ifndef NONLINEAR_ANALYSIS_DRIVER_ENUM_NAME
#define NONLINEAR_ANALYSIS_DRIVER_ENUM_NAME "NONLINEAR_ANALYSIS_DRIVER"
#else
#error
#endif


namespace Driver
{
  
  DeclareEnumName(NONLINEAR_ANALYSIS_DRIVER, Driver::AnalysisDriverTypeEnum,
                  NONLINEAR_ANALYSIS_DRIVER_ENUM_ID, NONLINEAR_ANALYSIS_DRIVER_ENUM_NAME);
  
  
  class NonLinearAnalysisDriver : public Driver::AnalysisDriver
    {
public:
      /// constructor
      NonLinearAnalysisDriver(const unsigned int ID,
                              FESystem::FESystemController& controller);
      
      /// destructor
      virtual ~NonLinearAnalysisDriver();
      
      /// clears the data structures of this object
      virtual void clear();
      
      /// returns the current iteration number of the nonlinear solution
      inline 
        unsigned int currentNonlinearIterationNumber() const;
      
      
      /// calculate the residual
      /// @param sol solution iterate at which the residual is to be calculated
      /// @param res vector in which the residual will be returned
      void getResidualAtVector(NumericVector<double>& sol,
                               NumericVector<double>& res);
      
      
      /// calculate the jacobian
      /// @param sol solution iterate at which the residual is to be calculated
      /// @param jac matrix in which the jacobian will be returned
      void getJacobianAtVector(NumericVector<double>& sol,
                               SparseMatrix<double>& jac);
      
protected:
          
      /// this method initializes itself before the solutions are started for the driver
      virtual void initialize();
      
      /// this method solves for the current load case
      virtual void solveCurrentLoadCase();
      
      /// function will add the matrices and vectors for this analysis
      void addMatricesAndVectors();
     
      /// method sets the initial guess for the nonlinear solver
      void setInitialGuessForCurrentLoadCase();
      
/*       /// method applies boundary condition to the residual vector */
/*       void applyBoundaryConditionsToResidualForCurrentAnalysis(NumericVector<double>& vector); */
      
/*       /// method applies boundary condition to the residual vector */
/*       void applyBoundaryConditionsToJacobianForCurrentAnalysis(SparseMatrix<double>& matrix); */

      /// iteration data
      unsigned int iteration_number;
      
      /// total number of iterations
      unsigned int n_iterations;
    };
  
}
#endif // __fesystem_non_linear_analysis_driver_h__
