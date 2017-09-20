// $Id: NewtonNonlinearSolver.C,v 1.7.4.6 2008-04-06 04:08:24 manav Exp $

// FESystem includes
#include "Solvers/NewtonNonlinearSolver.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/LinearSolver.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Utilities/TimeLogs.h"
#include "Utilities/ParallelUtility.h"

// libmesh includes
#include "numerics/petsc_matrix.h"
#include "numerics/petsc_vector.h"


Solver::NewtonNonlinearSolver::NewtonNonlinearSolver
(const Solver::NonlinearSolverInfo& nonlinear_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::NonlinearSolver(nonlinear_info, linear_info),
iteration_number(0),
linear_solver(NULL)
{
  this->linear_solver = 
  Solver::createLinearSolver(linear_info).release();
}




Solver::NewtonNonlinearSolver::~NewtonNonlinearSolver()
{
  delete this->linear_solver;
  this->linear_solver = NULL;
}



void
Solver::NewtonNonlinearSolver::clear()
{
  this->iteration_number = FESystemNumbers::InvalidID;
  Solver::NonlinearSolver::clear();
  this->linear_solver->clear();
}




void 
Solver::NewtonNonlinearSolver::solve(NumericVector<double>& solution)
{
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("NewtonNonlinearSolver::solve()", "Solver");

  this->linear_solver->attachAnalysisDriver(this->analysis_driver);
  
  // create a vector to store the solution and delta solution at each iteration
  Vec petsc_delta_sol, petsc_new_sol;
  PetscVector<double>& sol_vec = dynamic_cast<PetscVector<double>&>(solution);
  
  PetscErrorCode ierr = 0;
  ierr = VecDuplicate(sol_vec.vec(), &petsc_delta_sol);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  ierr = VecDuplicate(sol_vec.vec(), &petsc_new_sol);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  PetscVector<double>  delta_solution(petsc_delta_sol), new_sol(petsc_new_sol);
  delta_solution.zero();
  new_sol = solution;
  
  // get the tolerance and max iterations from solver info
  const double nonlinear_tol = this->nonlinear_solver_info.getConvergenceTolerance();
  const unsigned int nonlinear_maxits = this->nonlinear_solver_info.getMaximumIterations();
  
  // next, iterate till the termination criterion is satisfied  
  double old_res_norm = 0.0, new_res_norm = 0.0, delta_sol_norm = 0.0;
  
  for (this->iteration_number = 0 ;this->iteration_number < nonlinear_maxits; 
       this->iteration_number++)
    {
    this->evaluateResidual(new_sol, *(this->residual));
    this->residual->close();
        
    // calculate the norm of rhs. if norm is less than tolerance, break the loop
    new_res_norm = this->residual->l2_norm();
    
      {
        std::ostringstream oss;
        oss << "rhs norm : "<< new_res_norm << std::endl;
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, oss.str());
      }
                                    
    if (this->iteration_number == 0 ||  new_res_norm < old_res_norm)
      {
      	old_res_norm = new_res_norm;
		solution = new_sol;
      
      if (new_res_norm < nonlinear_tol)
        {
	      FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, 
                                            "*** Nonlinear iterations converged through residual norm ***" );
	      break;
        }
      }
    else
      {
      // diverging solution
	  	old_res_norm = new_res_norm;
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, "*** Diverging Nonlinear iterations***");
//      break;
      }
    
    // zero the solution vector before going into the solution
    this->evaluateJacobian(new_sol, *(this->jacobian));
    this->jacobian->close();
    this->residual->close(); delta_solution.close();
    
    //Solve the linear system to calculate the delta_T for this iteration.
    delta_solution.zero();
    this->linear_solver->setSystemMatrix(*(this->jacobian));
    this->linear_solver->setPreconditionerMatrix(*(this->jacobian));
    this->linear_solver->solve(*(this->residual), delta_solution);
    
// {  
//      std::auto_ptr<PetscVector<double> > old_res(new PetscVector<double>(this->residual->size()));
//      std::auto_ptr<PetscVector<double> > new_vec(new PetscVector<double>(this->residual->size()));
//      std::auto_ptr<PetscVector<double> > est_res(new PetscVector<double>(this->residual->size()));
//      std::auto_ptr<PetscVector<double> > frac_delta(new PetscVector<double>(this->residual->size()));
//      
//      *(old_res.get()) = *(this->residual);
//
//      // now calculate res and estimate it at different X and see if the two are
//      // tangent or not
//      unsigned int n_steps = 20;
//      std::vector<double> fract(n_steps+1), estimate(n_steps+1), actual(n_steps+1);
//      
//      double fraction = 0.0, dx = 1.0 / n_steps,
//        lower = -1.0, upper = 1.0;
//      for (unsigned int i=0; i <= n_steps; i++)
//        {
//        // calculate the values of X
//        new_vec->close();
//        new_sol.close();
//        fraction = lower + dx * i * (upper - lower);
//        *new_vec = new_sol;
//        new_vec->close();
//        new_vec->add(-fraction, delta_solution);
//        new_vec->close();
//        // calculate the exact value of residual
//        this->residual->zero();
//        this->residual->close();
//        this->evaluateResidual(*new_vec, *(this->residual));
//        this->residual->close();
//                  
//        // calculate the estimated residual
//        est_res->zero();
//        est_res->close();
//        frac_delta->zero();
//        frac_delta->close();
//        
//        frac_delta->add(-fraction, delta_solution);
//        frac_delta->close();
//        
//        this->jacobian->multiply_vector(*(frac_delta.get()),
//                                        *(est_res.get()));
//        est_res->close();
//        
//        est_res->add(*(old_res.get()));
//        est_res->close();
//        
//        fract[i] = fraction;
//        actual[i] = this->residual->l2_norm();
//        estimate[i] = est_res->l2_norm();
//        }
//
////      for (unsigned int it=0; it < fract.size(); it++)
////        std::cout << fract[it] << "   " 
////          << actual[it] << "   "
////          << estimate[it] << "\n";
//      
//      PetscDraw draw;
//      PetscDrawLG lgdraw;
//      PetscDrawOpenX(FESystem::COMM_WORLD, PETSC_NULL, "residual_plot", PETSC_DECIDE, 
//                     PETSC_DECIDE, PETSC_DRAW_HALF_SIZE, PETSC_DRAW_HALF_SIZE, &draw);
//      PetscDrawSetPause(draw, -1);
//      PetscDrawLGCreate(draw, 2, &lgdraw);
//      double x[2],y[2];
//      for (unsigned int it=0; it < fract.size(); it++)
//        {
//        x[0]= fract[it]; x[1]=fract[it];
//        y[0] = actual[it]; y[1]= estimate[it];
//        PetscDrawLGAddPoint(lgdraw, x,y);
//        }
//      PetscDrawLGDraw(lgdraw);
//      PetscDrawPause(draw);
//      
//      PetscDrawDestroy(draw);
//    }
    
    // calculate the solution vector for next iteration
    delta_solution.close();
    new_sol.add(-1.0, delta_solution);
    new_sol.close(); // this now contains the new solution
    delta_sol_norm = delta_solution.l2_norm();
    
      {
        std::ostringstream oss;
        oss << "delta_sol norm : "<< delta_sol_norm << std::endl;
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, oss.str());
      }
    // calculate the norm of delta_T. if norm is less than tolerance, break the loop
    if (delta_sol_norm < nonlinear_tol)
      {
      solution = new_sol;
      FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, 
                                          "*** Nonlinear iterations converged through dX norm ***" );
      break;
      }
    }
  
  
  if (this->iteration_number == nonlinear_maxits) 
    {
      FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, 
                                          "*** Nonlinear iterations reached maximum iterations limit. ***");
    }
  
  
  // finally, before exiting, destroy the two vectors that were created
  ierr = VecDestroy(petsc_delta_sol);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  ierr = VecDestroy(petsc_new_sol);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
    ("NewtonNonlinearSolver::solve()", "Solver");
}

