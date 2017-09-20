// $Id: TransientSolver.h,v 1.1.2.3 2008-02-25 04:32:39 manav Exp $

#ifndef __fesystem_transient_solver_h__
#define __fesystem_transient_solver_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/FESystemSolverBase.h"

// Forward Declerations
template <typename DataType> class NumericVector;
template <typename DataType> class SparseMatrix;

namespace Solver
{
  class LinearSolverInfo;
  class NonlinearSolverInfo;
  class TransientSolverInfoBase;
  class LinearTransientSolverInfo;
  class NonlinearTransientSolverInfo;
}

#ifndef LINEAR_TRANSIENT_SOLVER_ENUM_ID
#define LINEAR_TRANSIENT_SOLVER_ENUM_ID 4
#else
#error
#endif

#ifndef LINEAR_TRANSIENT_SOLVER_ENUM_NAME
#define LINEAR_TRANSIENT_SOLVER_ENUM_NAME "LINEAR_TRANSIENT_SOLVER"
#else
#error
#endif



#ifndef NONLINEAR_TRANSIENT_SOLVER_ENUM_ID
#define NONLINEAR_TRANSIENT_SOLVER_ENUM_ID 5
#else
#error
#endif

#ifndef NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME
#define NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME "NONLINEAR_TRANSIENT_SOLVER"
#else
#error
#endif


namespace Solver
{
  DeclareEnumName(LINEAR_TRANSIENT_SOLVER, Solver::SolverClassEnum,
                  LINEAR_TRANSIENT_SOLVER_ENUM_ID,
                  LINEAR_TRANSIENT_SOLVER_ENUM_NAME);
  
  DeclareEnumName(NONLINEAR_TRANSIENT_SOLVER, Solver::SolverClassEnum,
                  NONLINEAR_TRANSIENT_SOLVER_ENUM_ID,
                  NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME);
  
  DeclareEnumClass(TransientSolverKindEnum);
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class TransientSolverBase : public FESystemSolverBase
  {
public:

    /// constructor
    TransientSolverBase(const Solver::TransientSolverInfoBase& info);
    
    // destructor
    virtual ~TransientSolverBase();
    
    /// this method clears the data structures of this object. This should be called 
    /// each time the used finishes using this object.
    virtual void clear();
    
    /// method to solve the eigen system
    virtual void solve() = 0;
    
    /// @returns the current time of the solver integration
    virtual double getSimulatedTime() = 0;

    /// @returns the current time of the solver integration
    virtual double getCurrentTime() = 0;

    /// @returns the current time step of the solver integration
    virtual double getCurrentStepSize() = 0;

    /// @returns the current time iteration of the solver integration
    virtual unsigned int getSimulatedIterationNumber() = 0;

    /// @returns the current time iteration of the solver integration
    virtual unsigned int getCurrentIterationNumber() = 0;
    
    /// this attaches the matrix to be used for jacobian evaluation
    void attachCoefficientMatrices(std::vector<SparseMatrix<double>* > & matrices);
    
    /// @returns the time step size of the solver
    double getTimeStepSize() const;

/*     /// this attaches the matrix to be used for jacobian evaluation */
/*     void attachCoefficientMatrixSensitivity(std::vector<SparseMatrix<double>* > & matrices); */

    /// set initial conditions for the analysis
    virtual void setInitialCondition(std::vector<NumericVector<double>*>& state_vecs) = 0;
    
/*     /// set initial conditions for the sensitivity analysis case */
/*     virtual void setSensitivityInitialCondition(std::vector<NumericVector<double>*>& state_vecs) =0; */

/*     /// sets simultaneous sensitivity on or off (off by default) */
/*     void setSimultaneousSensitivity(const bool sensitivity); */
    
protected:
      
/*     /// if simultaneous sensitivity is to be performed */
/*     bool simultaneous_sensitivity; */

    /// TransientSolverInfo that stores information about this solver
    const Solver::TransientSolverInfoBase& transient_solver_info;

    /// coefficient matrices stored according to the order of the time derivative of the
    /// solution vector for which they act as coefficient
    std::vector<SparseMatrix<double>*> coefficient_matrices;

/*     /// coefficient matrix sensitivities stored according to the order of the time derivative of the */
/*     /// solution vector for which they act as coefficient. This is used only if  */
/*     /// simultaneous sensitivity is turned on */
/*     std::vector<SparseMatrix<double>*> coefficient_matrix_sensitivity; */
  };
  
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class LinearTransientSolverBase : public TransientSolverBase
    {
public:
      
      /// constructor
      LinearTransientSolverBase(const Solver::LinearTransientSolverInfo& transient_info,
                                const Solver::LinearSolverInfo& linear_info);
      
      // destructor
      virtual ~LinearTransientSolverBase();
      
      /// sets if the matrices are time dependent or time independent
      void setCoefficientMatrixTimeDependence(const bool dependence);
      
      virtual void clear();
      
protected:
        
//      /// TransientSolverInfo that stores information about this solver
//      const Solver::LinearTransientSolverInfo& transient_solver_info;
      
      /// LinearSolverInfo that stores information about this solver
      const Solver::LinearSolverInfo& linear_solver_info;
      
      /// if the matrices are dependent on time or not
      bool time_independent_coefficient_matrices;
    };

  
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class NonlinearTransientSolverBase : public TransientSolverBase
    {
public:
      
      /// constructor
      NonlinearTransientSolverBase(const Solver::NonlinearTransientSolverInfo& transient_info,
                                   const Solver::NonlinearSolverInfo& nonlinear_info,
                                   const Solver::LinearSolverInfo& linear_info);
      
      // destructor
      virtual ~NonlinearTransientSolverBase();
      
      virtual void clear();
      
      /// this attaches the matrix to be used for jacobian evaluation
      void attachJacobianMatrix(SparseMatrix<double>& joc_mat);
      
      /// evaluates the residual vector at the specified solution
      /// @param sol vector at which the residual should be evaluated
      /// @param res vector in which the residual is stored
      void evaluateRHSFunction(const double current_time,
                               NumericVector<double>& sol,
                               NumericVector<double>& res);
      
      /// evaluates the residual vector at the specified solution
      /// @param sol vector at which the residual should be evaluated
      /// @param jac matrix in which the jacobian is stored
      void evaluateRHSJacobian(const double current_time,
                               NumericVector<double>& sol,
                               SparseMatrix<double>& jac);
      
protected:
      
//      /// TransientSolverInfo that stores information about this solver
//      const Solver::NonlinearTransientSolverInfo& transient_solver_info;
      
      /// LinearSolverInfo that stores information about this solver
      const Solver::NonlinearSolverInfo& nonlinear_solver_info;
      
      /// LinearSolverInfo that stores information about this solver
      const Solver::LinearSolverInfo& linear_solver_info;
      
      /// jacobian matrix
      SparseMatrix<double>* jacobian_matrix;
    };
  
  
}



#endif // __fesystem_transient_solver_h__
