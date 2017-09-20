// $Id: NewmarkTransientSolver.h,v 1.1.2.1 2007-05-15 20:38:53 manav Exp $

#ifndef __fesystem_newmark_transient_solver_h__
#define __fesystem_newmark_transient_solver_h__

// C++ includes

// FESystem includes
#include "Solvers/TransientSolver.h"
#include "Utilities/AutoptrVector.h"

// forward declerations 
namespace Solver
{
  class NewmarkTransientSolverInfo;
  
  namespace NewmarkSolver
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


#ifndef NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_ID
#define NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_ID 5
#else
#error
#endif

#ifndef NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_NAME
#define NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_NAME "NEWMARK_LINEAR_TRANSIENT_SOLVER"
#else
#error
#endif



#ifndef NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID
#define NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID 6
#else
#error
#endif

#ifndef NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME
#define NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME "NEWMARK_NONLINEAR_TRANSIENT_SOLVER"
#else
#error
#endif






namespace Solver
{
  DeclareEnumName(NEWMARK_LINEAR_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_ID,
                  NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_NAME);
  
  
  DeclareEnumName(NEWMARK_NONLINEAR_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID,
                  NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_NAME);
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class NewmarkLinearTransientSolver : public Solver::LinearTransientSolverBase
    {
public:
      
      /// constructor
      NewmarkLinearTransientSolver(const Solver::LinearTransientSolverInfo& transient_info,
                                   const Solver::LinearSolverInfo& linear_info);
      
      // destructor
      virtual ~NewmarkLinearTransientSolver();
      
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

      std::vector<double> newmark_constants;
      
      /// internal state vectors, arranged in increasing order of 
      /// time derivative. The (n-1)th iterates are stored here
      FESystemUtility::AutoPtrVector<NumericVector<double> > current_state_vectors;

      /// internal state vectors, arranged in increasing order of 
      /// time derivative. The nth iterates are stored here
      FESystemUtility::AutoPtrVector<NumericVector<double> > new_state_vectors;

    };
  
  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class NewmarkNonlinearTransientSolver : public Solver::NonlinearTransientSolverBase
    {
public:
      
      /// constructor
      NewmarkNonlinearTransientSolver(const Solver::NonlinearTransientSolverInfo& transient_info,
                                      const Solver::NonlinearSolverInfo& nonlinear_info,
                                      const Solver::LinearSolverInfo& linear_info);
      
      // destructor
      virtual ~NewmarkNonlinearTransientSolver();
      
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
      
      friend void Solver::NewmarkSolver::getResidual(NumericVector<double>& sol,
                                                     NumericVector<double>& res,
                                                     void* ctx);
      
      friend void Solver::NewmarkSolver::getJacobian(NumericVector<double>& sol,
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
      
      std::vector<double> newmark_constants;
      
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

  
  
  namespace NewmarkSolver
    {
      /// this calculates the constants of the Newmark solver. The order of the 
      /// solver and the input constants are obtained from the solver_info object.
      template <typename NewmarkInfoType>
      void initNewmarkConstants(const NewmarkInfoType& solver_info,
                                std::vector<double>& constants_vector);
      
      /// given all constants, this calculates the LHS matrix of the linear dynamic system 
      template <typename NewmarkInfoType>
        void initNewmarkStiffnessMatrix(const NewmarkInfoType& solver_info,
                                        std::vector<double>& constants_vector,
                                        std::vector<SparseMatrix<double>*>& coeff_matrices,
                                        SparseMatrix<double>& newmark_stiff_matrix);

      /// given all other matrices and vectors, 
      /// this calculated the RHS of the linear dynamic system 
      template <typename NewmarkInfoType>
      void initNewmarkForceVector(const NewmarkInfoType& solver_info,
                                  std::vector<double>& constant_vector,
                                  std::vector<SparseMatrix<double>*>& coeff_matrices,
                                  std::vector<NumericVector<double>*>& state_vectors,
                                  NumericVector<double>& external_force_vec,
                                  NumericVector<double>& scratch_vec1,
                                  NumericVector<double>& scratch_vec2,
                                  NumericVector<double>& newmark_force_vec);

      
      
      /// this initializes the state vector derivatives at the new iterate. This 
      /// method is to be called after the new zeroth order vector has been calculated, 
      /// so that all new derivtives can be calculated. The new zeroth order vector
      /// should be passed as the first element of the new_state_vector. When this function
      /// returns, the previous state vector will be updated with the new state vectors, 
      /// and the new state vectors will all be zeroed
      template <typename NewmarkInfoType>
      void updateStateVectors(const NewmarkInfoType& solver_info,
                              std::vector<double>& constant_vector,
                              std::vector<NumericVector<double>*>& previous_state_vector,
                              std::vector<NumericVector<double>*>& new_state_vector);
      
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


#endif // __fesystem_newmark_transient_solver_h__

