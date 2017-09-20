// $Id: PetscNonlinearSolver.C,v 1.7.4.4 2008-08-21 00:56:46 manav Exp $

// FESystem includes
#include "Solvers/PetscNonlinearSolver.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/PetscLinearSolver.h"
#include "Utilities/TimeLogs.h"
#include "Utilities/ParallelUtility.h"

// libMesh includes
#include "petsc_vector.h"
#include "petsc_matrix.h"


// forward declerations of functions that handle the jacobian and 
// residual function evaluations
PetscErrorCode 
FESystemNonlinearSolverFunctionEvaluation
(SNES snes, Vec x, Vec f, void *ctx);

PetscErrorCode 
FESystemNonlinearSolverJacobianEvaluation
(SNES snes, Vec x, Mat *A, Mat *B, MatStructure *flag, void *ctx);



Solver::PetscNonlinearSolver::PetscNonlinearSolver
(const Solver::NonlinearSolverInfo& nonlinear_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::NonlinearSolver(nonlinear_info, linear_info)
{
  // create the SNES context and set the KSP and PC options
  PetscErrorCode ierr = 0;
  ierr = SNESCreate(FESystem::COMM_WORLD, &(this->petsc_snes));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // set the options for the SNES
  switch (this->nonlinear_solver_info.getNonlinearSolverKindEnumID())
    {
    case NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_ID:
      {
        // set the solver type to line search, and also set the line
        // search type
        ierr = SNESSetType(this->petsc_snes, SNESLS);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
        
        ierr = SNESLineSearchSet(this->petsc_snes, SNESLineSearchCubic, PETSC_NULL);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
      }
      break;
      
    case NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_ID:
      {
        // set the solver type to line search, and also set the line
        // search type
        ierr = SNESSetType(this->petsc_snes, SNESLS);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
        
        ierr = SNESLineSearchSet(this->petsc_snes, SNESLineSearchQuadratic, PETSC_NULL);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
      }
      break;
      
    case NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_ID:
      {
        // set the solver type to line search, and also set the line
        // search type
        ierr = SNESSetType(this->petsc_snes, SNESLS);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
        
        ierr = SNESLineSearchSet(this->petsc_snes, SNESLineSearchNo, PETSC_NULL);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
      }
      break;
      
    case NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_ID:
      {
        // set the solver type to line search, and also set the line
        // search type
        ierr = SNESSetType(this->petsc_snes, SNESLS);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
        
        ierr = SNESLineSearchSet(this->petsc_snes, SNESLineSearchNoNorms, PETSC_NULL);
        CHKERRABORT(FESystem::COMM_WORLD,ierr);
      }
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  
  const double nonlinear_tol = this->nonlinear_solver_info.getConvergenceTolerance();
  const unsigned int nonlinear_maxits = this->nonlinear_solver_info.getMaximumIterations();
  
  ierr = SNESSetTolerances(this->petsc_snes, nonlinear_tol, nonlinear_tol, nonlinear_tol, 
                           nonlinear_maxits, 1000);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // set the default monitor for the solver
  ierr = SNESMonitorSet(this->petsc_snes,
                        SNESMonitorDefault, 
                        PETSC_NULL, PETSC_NULL);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
  
  
  this->linear_solver.reset(new Solver::PetscLinearSolver(linear_info));
  this->petsc_ksp = linear_solver->getKSP();
  this->petsc_pc = linear_solver->getPC();
  // get the KSP context from the SNES and set its options
  ierr = SNESSetKSP(this->petsc_snes, this->petsc_ksp);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
}




Solver::PetscNonlinearSolver::~PetscNonlinearSolver()
{
  // destroy the nonlinear context
  PetscErrorCode ierr = SNESDestroy(this->petsc_snes);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
}




void
Solver::PetscNonlinearSolver::clear()
{
  Solver::NonlinearSolver::clear();
}



void
Solver::PetscNonlinearSolver::solve(NumericVector<double>& solution )
{
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("PetscNonlinearSolver::solve()", "Solver");

  Assert(this->jacobian != NULL, ExcInternalError());
  Assert(this->residual != NULL, ExcInternalError());

  this->jacobian->close();
  this->residual->close();
  solution.close();
  
  Mat petsc_jac = dynamic_cast<PetscMatrix<double>*>(this->jacobian)->mat();
  Vec petsc_sol = dynamic_cast<PetscVector<double>&>(solution).vec();
  Vec petsc_res = dynamic_cast<PetscVector<double>*>(this->residual)->vec();
  
  // set the SNES function and jacobian routines
  PetscErrorCode ierr = SNESSetFunction(this->petsc_snes, petsc_res, 
                                        FESystemNonlinearSolverFunctionEvaluation,
                                        this);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  ierr = SNESSetJacobian(this->petsc_snes, petsc_jac, petsc_jac,
                         FESystemNonlinearSolverJacobianEvaluation,
                         this);
  //   ierr = SNESSetJacobian(this->petsc_snes, petsc_jac, petsc_jac,
  // 			 SNESDefaultComputeJacobian,
  //   			 this);
  
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // now solve the system
  ierr = SNESSolve(this->petsc_snes, PETSC_NULL, petsc_sol);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  
  // get the convergence reason
  SNESConvergedReason reason; 
  ierr = SNESGetConvergedReason(this->petsc_snes, &reason);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // print out information about convergence
  switch (reason)
    {
    case SNES_CONVERGED_FNORM_ABS: //        =  2, /* F < F_minabs */
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,"SNES Converged due to F < F_minabs");
      }
      break;
      
    case SNES_CONVERGED_FNORM_RELATIVE:  //  =  3, /* F < F_mintol*F_initial */
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Converged due to F < F_mintol * F_initial");
      }
      break;
      
    case SNES_CONVERGED_PNORM_RELATIVE: //    =  4, /* step size small */
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Converged due to Small Step Size");
      }
      break;
      
    case SNES_CONVERGED_TR_DELTA: //        =  7,
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Converged due to TR_DELTA");
      }
      break;
      
      /* diverged */
    case SNES_DIVERGED_FUNCTION_DOMAIN:   // = -1,  
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Diverged in function domain");
      }
      break;
      
    case SNES_DIVERGED_FUNCTION_COUNT: //    = -2,  
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Diverged function count");
      }
      break;
      
    case SNES_DIVERGED_FNORM_NAN: //          = -4, 
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES obtained NAN ||F||");
      }
      break;
      
    case SNES_DIVERGED_MAX_IT: //             = -5,
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Did not converge in maximum iterations");
      }
      break;
      
    case SNES_DIVERGED_LS_FAILURE: //         = -6,
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Line Search Did not converge");
      }
      break;
      
    case SNES_DIVERGED_LOCAL_MIN: //         = -8,  /* || J^T b || is small, implies converged to local minimum of F() */
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Stuck in a local minima");
      }
      break;
      
    case SNES_CONVERGED_ITERATING: //         =  0
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            "SNES Converged Iterating");
      }
      break;
      
    default:
      abort();
      break;
    }

  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
    ("PetscNonlinearSolver::solve()", "Solver");
}




PetscErrorCode 
FESystemNonlinearSolverFunctionEvaluation
(SNES snes, Vec x, Vec f, void *ctx)
{
  // unused parameters
  (void) snes;
  
  // get the pointer to the solver, and ask it to calculate the 
  // residual at the function value given
  Solver::NonlinearSolver *solver = 
  static_cast<Solver::NonlinearSolver*> (ctx);
  
  PetscVector<double> sol_vec(x);
  PetscVector<double> res_vec(f);
  
  solver->evaluateResidual(sol_vec, res_vec);
  
  return 0;
}



PetscErrorCode 
FESystemNonlinearSolverJacobianEvaluation
(SNES snes, Vec x, Mat *A, Mat *B, MatStructure *flag, void *ctx)
{
  // unused parameter
  (void) snes;
  (void) B;

  // get the pointer to the solver, and ask it to calculate the 
  // residual at the function value given
  Solver::NonlinearSolver *solver = 
  static_cast<Solver::NonlinearSolver*> (ctx);
  
  PetscVector<double> sol_vec(x);
  PetscMatrix<double> jac_mat(*A);
  
  solver->evaluateJacobian(sol_vec, jac_mat);
  
  *flag = SAME_NONZERO_PATTERN;
  return 0;
}
