// $Id: TransientAnalysisDriver.h,v 1.1.2.6 2008/02/25 04:18:27 manav Exp $

#ifndef __fesystem_transient_analysis_driver_h__
#define __fesystem_transient_analysis_driver_h__

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


#ifndef LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID
#define LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID 4
#else
#error
#endif

#ifndef LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_NAME
#define LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_NAME "LINEAR_TRANSIENT_ANALYSIS_DRIVER"
#else
#error
#endif


#ifndef NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID
#define NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID 5
#else
#error
#endif

#ifndef NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_NAME
#define NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_NAME "NONLINEAR_TRANSIENT_ANALYSIS_DRIVER"
#else
#error
#endif



#ifndef TRANSIENT_C1_MATRIX_ENUM_ID
#define TRANSIENT_C1_MATRIX_ENUM_ID 11
#else
#error
#endif


#ifndef TRANSIENT_C1_MATRIX_ENUM_NAME
#define TRANSIENT_C1_MATRIX_ENUM_NAME "TRANSIENT_C1_MATRIX"
#else
#error
#endif


#ifndef TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_ID
#define TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_ID 12
#else
#error
#endif


#ifndef TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_NAME
#define TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_NAME "TRANSIENT_C1_MATRIX_SENSITIVITY"
#else
#error
#endif



#ifndef TRANSIENT_C2_MATRIX_ENUM_ID
#define TRANSIENT_C2_MATRIX_ENUM_ID 13
#else
#error
#endif


#ifndef TRANSIENT_C2_MATRIX_ENUM_NAME
#define TRANSIENT_C2_MATRIX_ENUM_NAME "TRANSIENT_C2_MATRIX"
#else
#error
#endif


#ifndef TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_ID
#define TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_ID 14
#else
#error
#endif


#ifndef TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_NAME
#define TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_NAME "TRANSIENT_C2_MATRIX_SENSITIVITY"
#else
#error
#endif



namespace Driver
{
  
  DeclareEnumName(LINEAR_TRANSIENT_ANALYSIS_DRIVER, Driver::AnalysisDriverTypeEnum,
                  LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID, 
                  LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_NAME);
  
  DeclareEnumName(NONLINEAR_TRANSIENT_ANALYSIS_DRIVER, Driver::AnalysisDriverTypeEnum,
                  NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID, 
                  NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_NAME);
  
  DeclareEnumName(TRANSIENT_C1_MATRIX, Driver::AnalysisDriverQtyEnum,
                  TRANSIENT_C1_MATRIX_ENUM_ID, 
                  TRANSIENT_C1_MATRIX_ENUM_NAME);
  
  DeclareEnumName(TRANSIENT_C1_MATRIX_SENSITIVITY, Driver::AnalysisDriverQtyEnum,
                  TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_ID, 
                  TRANSIENT_C1_MATRIX_SENSITIVITY_ENUM_NAME);
  
  DeclareEnumName(TRANSIENT_C2_MATRIX, Driver::AnalysisDriverQtyEnum,
                  TRANSIENT_C2_MATRIX_ENUM_ID, 
                  TRANSIENT_C2_MATRIX_ENUM_NAME);
  
  DeclareEnumName(TRANSIENT_C2_MATRIX_SENSITIVITY, Driver::AnalysisDriverQtyEnum,
                  TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_ID, 
                  TRANSIENT_C2_MATRIX_SENSITIVITY_ENUM_NAME);
  
  
  /// this forms the basis of all transient analysis drivers.
  class TransientAnalysisDriver : public Driver::AnalysisDriver
  {
  public:
    /// constructor
    TransientAnalysisDriver(const unsigned int ID,
                            FESystem::FESystemController& controller,
                            const unsigned int driver_enum_ID);
    
    /// destructor
    virtual ~TransientAnalysisDriver();
    
    /// @returns the current time for calculation of system quantitites. This 
    /// will not necessarily be the same as the simulation time, since the solver 
    /// might be going through sub-iterations for each time step. For actual simulation time
    /// that has been traversed, use the method currentSimulationTime();
    double getCurrentAnalysisTime() const;
    
    /// current simulation time is the time that has been solved for the system
    double getSimulatedTime() const;
    
    /// @returns the iteration number for which solutions are available
    unsigned int getSimulatedIterationNumber() const;
    
    /// @returns the iteration number for which current quantities are being calculated
    unsigned int getCurrentIterationNumber() const;
    
    /// clears the data structures of this object
    virtual void clear();
    
    void setNewIteration(const double time, 
                         const unsigned int iter_num,
                         const std::vector<NumericVector<double>*>& iter_sols);
    
    /// calculate the coefficient matrices. The location of the matrix is also the order
    /// of the time derivative vector
    virtual void getCoefficientMatrices(const double time,
                                        std::vector<SparseMatrix<double>* >& matrices) =0;
    
    
    /// calculate the coefficient matrices. The location of the matrix is also the order
    /// of the time derivative vector
    virtual void getForceVector(const double time, NumericVector<double>& vector) =0;
    
    
  protected:
    
    /// method sets the initial guess for the nonlinear solver
    void setInitialConditionForCurrentLoadCase(std::vector<NumericVector<double>*>& vector);
    
    /// current time of analysis 
    double current_analysis_time; 
    
    /// vector of time instants for the transient iterations
    std::vector<double> time_values;
  };
  
  
  
  /// this is the base for nonlinear transient analysis drivers of various orders
  class LinearTransientAnalysisDriver : public Driver::TransientAnalysisDriver
  {
  public:
    /// constructor
    LinearTransientAnalysisDriver(const unsigned int ID,
                                  FESystem::FESystemController& controller);
    
    /// destructor
    virtual ~LinearTransientAnalysisDriver();
    
    
    /// calculate the coefficient matrices. The location of the matrix is also the order
    /// of the time derivative vector
    virtual void getCoefficientMatrices(const double time,
                                        std::vector<SparseMatrix<double>* >& matrices);
    
    
    /// calculate the coefficient matrices. The location of the matrix is also the order
    /// of the time derivative vector
    virtual void getForceVector(const double time, NumericVector<double>& vector);
    
    
  protected:
    
    /// this method initializes itself before the solutions are started for the driver
    virtual void initialize();

    /// this method solves for the current load case
    virtual void solveCurrentLoadCase();
    
    /// function will add the matrices and vectors for this analysis
    virtual void addMatricesAndVectors();
  };
  
  
  
  
  /// this is the base for nonlinear transient analysis drivers of various orders
  class NonlinearTransientAnalysisDriver : public Driver::TransientAnalysisDriver
  {
  public:
    /// constructor
    NonlinearTransientAnalysisDriver(const unsigned int ID,
                                     FESystem::FESystemController& controller);
    
    /// destructor
    virtual ~NonlinearTransientAnalysisDriver();
    
    
    /// calculate the coefficient matrices. The location of the matrix is also the order
    /// of the time derivative vector
    void getJacobianMatrix(const double time,
                           std::vector<NumericVector<double>*>& sol_vec,
                           SparseMatrix<double>& matrix);
    
    
    /// calculate the coefficient matrices. The location of the matrix is also the order
    /// of the time derivative vector
    void getCoefficientMatrices(const double time,
                                std::vector<NumericVector<double>*>& sol_vec,
                                std::vector<SparseMatrix<double>* >& matrices);
    
    /// calculates the force vector at the specified time and solution vectors specified 
    /// in vecs
    void getForceVector(const double time,
                        std::vector<NumericVector<double>* >& vecs,                                                                 
                        NumericVector<double>& res);
    
    /// calculates the nonlinear residual at the specified time and solution vectors specified 
    /// in vecs
    void getResidualVector(const double time,
                           std::vector<NumericVector<double>* >& vecs,                                                                 
                           NumericVector<double>& res);
    
    /// calculates the force vector for sensitivity analysis problems
    virtual void getForceVector(const double time, NumericVector<double>& vec);
    
    /// calculates the coefficient matrices for sensitivity analysis for the 
    /// specified time
    virtual void getCoefficientMatrices(const double time, 
                                        std::vector<SparseMatrix<double>* >& matrices);
    
  protected:
    
    /// this method initializes itself before the solutions are started for the driver
    virtual void initialize();

    /// this method solves for the current load case
    virtual void solveCurrentLoadCase();
    
    /// function will add the matrices and vectors for this analysis
    virtual void addMatricesAndVectors();
  };
  
}



#endif // __fesystem_transient_analysis_driver_h__


