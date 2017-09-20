// $Id: PetscNonlinearSolver.h,v 1.7 2006-12-07 01:32:27 manav Exp $

#ifndef __fesystem_petsc_nonlinear_solver_h__
#define __fesystem_petsc_nonlinear_solver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/NonlinearSolver.h"


// Petsc includes 
#include "petscsnes.h"

#ifndef NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_ID
#define NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_ID 2
#else
#error
#endif

#ifndef NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_NAME
#define NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_NAME "NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH"
#else
#error
#endif



#ifndef NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_ID
#define NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_ID 3
#else
#error
#endif

#ifndef NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_NAME
#define NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_NAME "NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH"
#else
#error
#endif



#ifndef NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_ID
#define NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_ID 4
#else
#error
#endif

#ifndef NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_NAME
#define NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_NAME "NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH"
#else
#error
#endif

#ifndef NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_ID
#define NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_ID 5
#else
#error
#endif

#ifndef NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_NAME
#define NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_NAME "NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH"
#else
#error
#endif

namespace Solver
{
  // forward declerations
  class PetscLinearSolver;
  
  
  DeclareEnumName(NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH, Solver::NonlinearSolverKindEnum,
                  NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_ID,
                  NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_NAME);
  
  DeclareEnumName(NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH, Solver::NonlinearSolverKindEnum,
                  NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_ID,
                  NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_NAME);

  DeclareEnumName(NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH, Solver::NonlinearSolverKindEnum,
                  NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_ID,
                  NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_NAME);

  DeclareEnumName(NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH, Solver::NonlinearSolverKindEnum,
                  NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_ID,
                  NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_NAME);
  
  /// this class provides an interface to a slover for the solution of a 
  /// nonlinear set of equations. This inherits from the NonlinearSolver
  /// class. It takes a pointer to the analysis driver object that owns this 
  /// object
  class PetscNonlinearSolver:public NonlinearSolver
    {
public:
      
      /// constructor
      PetscNonlinearSolver(const Solver::NonlinearSolverInfo& nonlinear_info,
                           const Solver::LinearSolverInfo& linear_info);
      
      /// destructor
      ~PetscNonlinearSolver();
      
      
      /// method clear the data structures of this object. This should be called each time
      /// the user finishes using this object.
      virtual void clear();
      
      /// function to perform solution steps.
      /// @param solution vector in which the solution will be returned. 
      /// This should be initialized to the 
      /// initial guess
      virtual void solve(NumericVector<double>& solution);
      
      /// this returns the itertion number of the solution process
      virtual inline unsigned int getCurrentIterationNumber() const;
      
      
protected:
        
      /// PETSc nonlinear solver context
      SNES petsc_snes;
     
      /// Petsc KSP context 
      KSP petsc_ksp;
      
      /// Petsc PC context
      PC petsc_pc;
      
      /// PETSc linear solver
      std::auto_ptr<Solver::PetscLinearSolver> linear_solver;
    };
}


inline
unsigned int 
Solver::PetscNonlinearSolver::getCurrentIterationNumber() const
{
  int iter_num = 0;
  
  PetscErrorCode ierr = SNESGetIterationNumber(this->petsc_snes, &iter_num);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  return iter_num;
}




#endif // __fesystem_petsc_nonlinear_solver_h__
