// $Id: EulerTransientSolver.h,v 1.1.2.1 2008-02-25 04:31:21 manav Exp $

#ifndef __fesystem_euler_transient_solver_h__
#define __fesystem_euler_transient_solver_h__

// C++ includes

// FESystem includes
#include "Solvers/TransientSolver.h"
#include "Utilities/AutoptrVector.h"

// forward declerations 
namespace Solver
{
  class EulerTransientSolverInfo;
  
  namespace EulerSolver
  {
    void getResidual(NumericVector<double>& sol,
                     NumericVector<double>& res,
                     void* ctx);
    void getJacobian(NumericVector<double>& sol,
                     SparseMatrix<double>& jac,
                     void* ctx);
  }
}

template <typename T> class NumericVector;
template <typename T> class SparseVector;


#ifndef EULER_LINEAR_TRANSIENT_SOLVER_ENUM_ID
#define EULER_LINEAR_TRANSIENT_SOLVER_ENUM_ID 7
#else
#error
#endif

#ifndef EULER_LINEAR_TRANSIENT_SOLVER_ENUM_NAME
#define EULER_LINEAR_TRANSIENT_SOLVER_ENUM_NAME "EULER_LINEAR_TRANSIENT_SOLVER"
#else
#error
#endif



#ifndef EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID
#define EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID 8
#else
#error
#endif

#ifndef EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME
#define EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME "EULER_NONLINEAR_TRANSIENT_SOLVER"
#else
#error
#endif






namespace Solver
{
  DeclareEnumName(EULER_LINEAR_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  EULER_LINEAR_TRANSIENT_SOLVER_ENUM_ID,
                  EULER_LINEAR_TRANSIENT_SOLVER_ENUM_NAME);
  
  
  DeclareEnumName(EULER_NONLINEAR_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID,
                  EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME);
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class EulerLinearTransientSolver : public Solver::LinearTransientSolverBase
  {
  public:
    
    /// constructor
    EulerLinearTransientSolver(const Solver::LinearTransientSolverInfo& transient_info,
                                 const Solver::LinearSolverInfo& linear_info);
    
    // destructor
    virtual ~EulerLinearTransientSolver();
    
    /// this method clears the data structures of this object. This should be called 
    /// each time the used finishes using this object.
    virtual void clear();
    
    /// method to solve the eigen system
    virtual void solve();
    
    /// @returns the current time of the solver integration   
    virtual double getSimulatedTime(); 
    
    /// @returns the current time of the solver integration
    virtual double getCurrentTime();
    
    /// @returns the current time step of the solver integration
    virtual double getCurrentStepSize();
    
    /// @returns the current time iteration of the solver integration 
    virtual unsigned int getSimulatedIterationNumber(); 
    
    /// @returns the current time iteration of the solver integration 
    virtual unsigned int getCurrentIterationNumber(); 
    
    /// set initial conditions for the analysis 
    virtual void setInitialCondition(std::vector<NumericVector<double>*>& state_vecs); 
    
  protected:
    
    /// simulated time 
    double simulated_time;
    
    /// current time
    double current_time;
    
    /// current step size
    double time_step_size;
    
    /// simulated iteration number
    unsigned int simulated_iteration_number;
    
    /// linear solver
    std::auto_ptr<Solver::LinearSolver> linear_solver;
    
    /// internal state vectors, arranged in increasing order of 
    /// time derivative. The (n-1)th iterates are stored here
    FESystemUtility::AutoPtrVector<NumericVector<double> > current_state_vectors;
  };
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class EulerNonlinearTransientSolver : public Solver::NonlinearTransientSolverBase
  {
  public:
    
    /// constructor
    EulerNonlinearTransientSolver(const Solver::NonlinearTransientSolverInfo& transient_info,
                                    const Solver::NonlinearSolverInfo& nonlinear_info,
                                    const Solver::LinearSolverInfo& linear_info);
    
    // destructor
    virtual ~EulerNonlinearTransientSolver();
    
    /// this method clears the data structures of this object. This should be called 
    /// each time the used finishes using this object.
    virtual void clear();
    
    /// method to solve the eigen system
    virtual void solve();
    
    /// @returns the current time of the solver integration   
    virtual double getSimulatedTime(); 
    
    /// @returns the current time of the solver integration
    virtual double getCurrentTime();
    
    /// @returns the current time step of the solver integration
    virtual double getCurrentStepSize();
    
    /// @returns the current time iteration of the solver integration 
    virtual unsigned int getSimulatedIterationNumber(); 
    
    /// @returns the current time iteration of the solver integration 
    virtual unsigned int getCurrentIterationNumber(); 
    
    /// set initial conditions for the analysis 
    virtual void setInitialCondition(std::vector<NumericVector<double>*>& state_vecs); 
    
    friend void Solver::EulerSolver::getResidual(NumericVector<double>& sol,
                                                   NumericVector<double>& res,
                                                   void* ctx);
    
    friend void Solver::EulerSolver::getJacobian(NumericVector<double>& sol,
                                                   SparseMatrix<double>& jac,
                                                   void* ctx);
    
    
  protected:
    
    /// simulated time 
    double simulated_time;
    
    /// current time
    double current_time;
    
    /// current step size
    double time_step_size;
    
    /// simulated iteration number
    unsigned int simulated_iteration_number;
    
    /// nonlinear solver
    std::auto_ptr<Solver::NonlinearSolver> nonlinear_solver;
        
    /// internal state vectors, arranged in increasing order of 
    /// time derivative. The (n-1)th iterates are stored here
    FESystemUtility::AutoPtrVector<NumericVector<double> > last_iter_state_vectors;
    
    /// internal state vectors, arranged in increasing order of 
    /// time derivative. The (n-1)th iterates are stored here
    FESystemUtility::AutoPtrVector<NumericVector<double> > current_state_vectors;
    
    /// internal state vectors, arranged in increasing order of 
    /// time derivative. The nth iterates are stored here
    FESystemUtility::AutoPtrVector<NumericVector<double> > new_state_vectors;
    
  };
  
  
  
  namespace EulerSolver
  {
        
    /// given all other matrices and vectors, 
    /// this calculated the RHS of the linear dynamic system 
    template <typename EulerInfoType>
    void initEulerForceVector(const EulerInfoType& solver_info,
                                std::vector<double>& constant_vector,
                                std::vector<SparseMatrix<double>*>& coeff_matrices,
                                std::vector<NumericVector<double>*>& state_vectors,
                                NumericVector<double>& external_force_vec,
                                NumericVector<double>& scratch_vec1,
                                NumericVector<double>& scratch_vec2,
                                NumericVector<double>& newmark_force_vec);
    
    
    /// function to be called by the nonlinear solver for the residual function
    void getResidual(NumericVector<double>& sol,
                     NumericVector<double>& res,
                     void* ctx);
    
    /// function to be called by the nonlinear solver for a jacobian matrix
    void getJacobian(NumericVector<double>& sol,
                     SparseMatrix<double>& jac,
                     void* ctx);
    
  }
}


#endif // __fesystem_euler_transient_solver_h__

