// $Id: EigenProblemAnalysisDriver.h,v 1.6 2006/11/13 00:09:25 manav Exp $

#ifndef __fesystem_eigen_problem_analysis_driver_h__
#define __fesystem_eigen_problem_analysis_driver_h__

// C++ includes


// FESystem includes
#include "AnalysisDriver/AnalysisDriver.h"

// libMesh includes


namespace FESystem
{
  class FESystemController;
}


#ifndef EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_ID
#define EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_ID 3
#else
#error
#endif

#ifndef EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_NAME
#define EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_NAME "EIGENPROBLEM_ANALYSIS_DRIVER"
#else
#error
#endif


namespace Driver
{

  DeclareEnumName(EIGENPROBLEM_ANALYSIS_DRIVER, Driver::AnalysisDriverTypeEnum,
                  EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_ID, EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_NAME);

  /// this class defines an analysis driver for a generic eigen problem. 
  /// Its basic purpose is to have 
  /// the basic solution routine for both analysis and sensitivity analysis for an eigenproblem.
  class EigenProblemAnalysisDriver: public Driver::AnalysisDriver
  {
public:
    /// constructor 
    EigenProblemAnalysisDriver(const unsigned int ID,
                               FESystem::FESystemController& controller);
    
    /// destructor 
    virtual ~EigenProblemAnalysisDriver();
    
protected:
    
    /// this method initializes itself before the solutions are started for the driver
    virtual void initialize();
    
    /// this method solves for the current load case
    virtual void solveCurrentLoadCase();
    
    /// function will add the matrices and vectors for this analysis
    virtual void addMatricesAndVectors();
		
    /// method to apply boundary conditions for current load case
    void applyBoundaryConditionsForCurrentAnalysis(SparseMatrix<double>& A_matrix,
                                                   SparseMatrix<double>* B_matrix);
    
    void getSubmatrixAfterBoundaryCondition(SparseMatrix<double>& A_matrix,
                                            SparseMatrix<double>* B_matrix,
                                            SparseMatrix<double>& A_matrix_sub,
                                            SparseMatrix<double>* B_matrix_sub);
    
    
  };
}


#endif // __fesystem_eigen_problem_analysis_driver_h__
