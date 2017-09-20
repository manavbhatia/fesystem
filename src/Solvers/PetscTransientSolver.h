// $Id: PetscTransientSolver.h,v 1.1.2.2 2007-05-08 05:19:13 manav Exp $

#ifndef __fesystem_petsc_transient_solver_h__
#define __fesystem_petsc_transient_solver_h__

// C++ includes

// FESystem includes
#include "Solvers/TransientSolver.h"

// petsc includes
#include "petscts.h"


#ifndef PETSC_BEULER_TRANSIENT_SOLVER_ENUM_ID
#define PETSC_BEULER_TRANSIENT_SOLVER_ENUM_ID 1
#else
#error
#endif

#ifndef PETSC_BEULER_TRANSIENT_SOLVER_ENUM_NAME
#define PETSC_BEULER_TRANSIENT_SOLVER_ENUM_NAME "PETSC_BEULER_TRANSIENT_SOLVER"
#else
#error
#endif



#ifndef PETSC_FEULER_TRANSIENT_SOLVER_ENUM_ID
#define PETSC_FEULER_TRANSIENT_SOLVER_ENUM_ID 2
#else
#error
#endif

#ifndef PETSC_FEULER_TRANSIENT_SOLVER_ENUM_NAME
#define PETSC_FEULER_TRANSIENT_SOLVER_ENUM_NAME "PETSC_FEULER_TRANSIENT_SOLVER"
#else
#error
#endif


#ifndef PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER_ENUM_ID
#define PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER_ENUM_ID 3
#else
#error
#endif

#ifndef PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER_ENUM_NAME
#define PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER_ENUM_NAME "PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER"
#else
#error
#endif



#ifndef PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER_ENUM_ID
#define PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER_ENUM_ID 4
#else
#error
#endif

#ifndef PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER_ENUM_NAME
#define PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER_ENUM_NAME "PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER"
#else
#error
#endif





namespace Solver
{
  DeclareEnumName(PETSC_BEULER_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  PETSC_BEULER_TRANSIENT_SOLVER_ENUM_ID,
                  PETSC_BEULER_TRANSIENT_SOLVER_ENUM_NAME);
  

  DeclareEnumName(PETSC_FEULER_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  PETSC_FEULER_TRANSIENT_SOLVER_ENUM_ID,
                  PETSC_FEULER_TRANSIENT_SOLVER_ENUM_NAME);


  DeclareEnumName(PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER_ENUM_ID,
                  PETSC_RUNGE_KUTTA_TRANSIENT_SOLVER_ENUM_NAME);

  
  DeclareEnumName(PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER, Solver::TransientSolverKindEnum,
                  PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER_ENUM_ID,
                  PETSC_CRANK_NICHOLSON_TRANSIENT_SOLVER_ENUM_NAME);
  

  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class PetscNonlinearTransientSolver : public Solver::NonlinearTransientSolverBase
    {
public:
      
      /// constructor
      PetscNonlinearTransientSolver(const Solver::NonlinearTransientSolverInfo& transient_info,
                                    const Solver::NonlinearSolverInfo& nonlinear_info,
                                    const Solver::LinearSolverInfo& linear_info);
      
      // destructor
      virtual ~PetscNonlinearTransientSolver();
      
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
        
        /// petsc transient solver context
        TS petsc_ts;
        
      /// PETSc nonlinear solver context
      SNES petsc_snes;
      
      /// PETSc KSP solver for the nonlinear solver
      KSP petsc_ksp;
      
      /// a PC context for the linear solver
      PC petsc_pc;
    };
  
}


#endif // __fesystem_petsc_transient_solver_h__
