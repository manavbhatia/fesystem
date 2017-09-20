// $Id: NewtonNonlinearSolver.h,v 1.6 2006-12-07 01:32:27 manav Exp $

#ifndef __fesystem_newton_nonlinear_solver_h__
#define __fesystem_newton_nonlinear_solver_h__


// FESystem includes
#include "Solvers/NonlinearSolver.h"

namespace Solver
{
  class LinearSolver;
}

#ifndef FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_ID
#define FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_ID 1
#else
#error
#endif

#ifndef FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_NAME
#define FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_NAME "FESYSTEM_NEWTON_NONLINEAR_SOLVER"
#else
#error
#endif



namespace Solver
{
    
  DeclareEnumName(FESYSTEM_NEWTON_NONLINEAR_SOLVER, Solver::NonlinearSolverKindEnum,
                  FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_ID,
                  FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_NAME);
    
  class NewtonNonlinearSolver : public Solver::NonlinearSolver
    {
public:
      NewtonNonlinearSolver(const Solver::NonlinearSolverInfo& nonlinear_info,
                            const Solver::LinearSolverInfo& linear_info);
      
      ~NewtonNonlinearSolver();
      
      /// method to clear the data structures of this object. This should be called each time
      /// the user completes solution and does not need this object.
      virtual void clear();
      
      /// function to perform solution steps.
      /// @param solution vector in which the solution will be returned. 
      /// This should be initialized to the 
      /// initial guess
      virtual void solve(NumericVector<double>& solution);
      
      /// this returns the itertion number of the solution process
      virtual inline unsigned int getCurrentIterationNumber() const;
      
      
protected:
      
      /// iteration number
      unsigned int iteration_number;

      /// a linear solver is needed to solve the equations at each itertion
      Solver::LinearSolver *linear_solver;
      
    };
  
}





inline 
unsigned int
Solver::NewtonNonlinearSolver::getCurrentIterationNumber() const
{
  return this->iteration_number;
}


#endif// __fesystem_newton_nonlinear_solver_h__
