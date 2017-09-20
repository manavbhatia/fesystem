// $Id: LinearAnalysisDriver.h,v 1.9.6.1 2007/05/08 05:18:27 manav Exp $

#ifndef __fesystem_linear_analysis_driver_h__
#define __fesystem_linear_analysis_driver_h__

// C++ includes


// FESystem includes
#include "AnalysisDriver/AnalysisDriver.h"

// forward decleration
namespace FESystem
{
  class FESystemController;
}

#ifndef LINEAR_ANALYSIS_DRIVER_ENUM_ID
#define LINEAR_ANALYSIS_DRIVER_ENUM_ID 1
#else
#error
#endif

#ifndef LINEAR_ANALYSIS_DRIVER_ENUM_NAME
#define LINEAR_ANALYSIS_DRIVER_ENUM_NAME "LINEAR_ANALYSIS_DRIVER"
#else
#error
#endif


namespace Driver
{

  DeclareEnumName(LINEAR_ANALYSIS_DRIVER, Driver::AnalysisDriverTypeEnum,
                  LINEAR_ANALYSIS_DRIVER_ENUM_ID, LINEAR_ANALYSIS_DRIVER_ENUM_NAME);

  
  class LinearAnalysisDriver : public Driver::AnalysisDriver
  {
public:
    /// Constructor
    LinearAnalysisDriver(const unsigned int ID,
                         FESystem::FESystemController& controller);
    
    /// destructor
    virtual ~LinearAnalysisDriver();
    
protected:

    /// this method initializes itself before the solutions are started for the driver
    virtual void initialize();
      
    /// this method solves for the current load case
    virtual void solveCurrentLoadCase();
    
    /// function will add the matrices and vectors for this analysis
    virtual void addMatricesAndVectors();
		
/*     /// method to apply boundary conditions for the current load case */
/*     void applyBoundaryConditionsToRHSForCurrentAnalysis(NumericVector<double>& rhs); */
    
/*     /// method to apply boundary conditions for current load case */
/*     void applyBoundaryConditionsToSystemMatrixForCurrentAnalysis(SparseMatrix<double>& matrix); */
  };
}


#endif // __fesystem_linear_analysis_driver_h__
