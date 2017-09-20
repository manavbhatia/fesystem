// $Id: EulerTransientSolver.C,v 1.1.2.2 2008-08-21 20:44:25 manav Exp $

// FESystem includes
#include "Solvers/EulerTransientSolver.h"
#include "FESystem/FESystemController.h"
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/EulerTransientSolverInfo.h"
#include "Utilities/TimeLogs.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"



Solver::EulerLinearTransientSolver::EulerLinearTransientSolver
(const Solver::LinearTransientSolverInfo& transient_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::LinearTransientSolverBase(transient_info, linear_info),
simulated_time(0.0),
current_time(0.0),
time_step_size(0.0),
simulated_iteration_number(0)
{
}



Solver::EulerLinearTransientSolver::~EulerLinearTransientSolver()
{
  
}




void
Solver::EulerLinearTransientSolver::clear()
{
  this->simulated_time = 0.0;
  this->current_time = 0.0;
  this->time_step_size = 0.0;
  this->simulated_iteration_number = 0;
  this->current_state_vectors.clear();
  Solver::LinearTransientSolverBase::clear();
}





void
Solver::EulerLinearTransientSolver::solve()
{
  this->linear_solver.reset(Solver::createLinearSolver(this->linear_solver_info).release());
  this->linear_solver->attachAnalysisDriver(this->analysis_driver);
  
  const Solver::EulerLinearTransientSolverInfo& newmark_info =
  dynamic_cast<const Solver::EulerLinearTransientSolverInfo&>(this->transient_solver_info);
  
  unsigned int system_order = newmark_info.getSystemOrder();
  // for now, this implementation is only for 1st order system
  Assert(system_order == 1, ExcInternalError());
  
  // make sure that the solver has the necessary information
  Assert(this->coefficient_matrices.size() == (system_order+1), ExcInternalError());
  Assert(this->current_state_vectors.size() == (system_order+1), ExcInternalError());
  
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("EulerLinearTransientSolver::solve()", "Solver");
  
  
  // create the vectors to be used in this solution
  std::auto_ptr<NumericVector<double> > 
  external_force_vector(NumericVector<double>::build().release()),
  scratch_vec1(NumericVector<double>::build().release());
  
  external_force_vector->duplicate_vector(*(this->current_state_vectors[0]), false);
  scratch_vec1->duplicate_vector(*(this->current_state_vectors[0]), false);
    
  // init the iteration data 
  this->simulated_time = this->transient_solver_info.getInitialTime();
  this->current_time = this->simulated_time;
  this->simulated_iteration_number = 0;
  
  // the basic equation to be solved is 
  // dX/dt^(n) = inv(C_(n)) * (F_0 - sum_(i=0..n-1) (C_i * dX/dt^(i)))    
  
  while (this->getSimulatedTime() < this->transient_solver_info.getFinalTime())
    {
      if (this->simulated_iteration_number > 0)
        {
          scratch_vec1->zero();
          scratch_vec1->close();
          scratch_vec1->add(this->transient_solver_info.getInitialTimeStep(),
                            *(this->current_state_vectors[1]));
          scratch_vec1->add(1.0, *(this->current_state_vectors[0]));
          scratch_vec1->close();
          this->current_state_vectors[0]->zero();
          this->current_state_vectors[0]->close();
          this->current_state_vectors[0]->add(1.0, *scratch_vec1);
          this->current_state_vectors[0]->close();
        }
      
      
      // ask the driver to init the matrices. This needs to be done only once for a system where 
      // the matrices are constant with iterations
      if (!this->time_independent_coefficient_matrices || this->simulated_iteration_number == 0)
        {
          dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
          getCoefficientMatrices(this->getCurrentTime(),
                                 this->coefficient_matrices);
          
          this->coefficient_matrices[system_order]->close();
          this->linear_solver->setSystemMatrix(*(this->coefficient_matrices[system_order]));
          this->linear_solver->setPreconditionerMatrix(*(this->coefficient_matrices[system_order]));
          this->linear_solver->useSameMatricesBetweenSolves();
        }
            
      // ask the driver to init the force vector;
      dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
      getForceVector(this->getCurrentTime(),
                     *external_force_vector);
      external_force_vector->close();
      
      // calculate the RHS of the linear equation
      scratch_vec1->zero(); 
      scratch_vec1->close();
      this->coefficient_matrices[0]->close();
      this->coefficient_matrices[0]->multiply_vector(*(this->current_state_vectors[0]),
                                                     *scratch_vec1);
      scratch_vec1->close();
      external_force_vector->add(-1.0, *scratch_vec1);
      external_force_vector->close();

      // now solve for the new iterate
      this->linear_solver->solve(*external_force_vector, 
                                 *(this->current_state_vectors[1]));
            
      this->current_state_vectors[1]->close();
      
      
      
      // increment the current time and the iteration number
      this->current_time += this->transient_solver_info.getInitialTimeStep();
      this->simulated_time = 
      this->getCurrentTime() - this->transient_solver_info.getInitialTimeStep();
      this->simulated_iteration_number++;
      
      // tell the analysis driver to record the iteration increment
      dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
      setNewIteration(this->getSimulatedTime(),
                      this->getSimulatedIterationNumber(),
                      this->current_state_vectors.getReference());
      
    }
  
  this->linear_solver.reset();
  
  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
  ("EulerLinearTransientSolver::solve()", "Solver"); 
}



void
Solver::EulerLinearTransientSolver::setInitialCondition
(std::vector<NumericVector<double>*>& state_vecs)
{
  // make sure that the current solutions have not been initialized
  Assert(this->current_state_vectors.size() == 0, ExcInternalError());
  
  unsigned int system_order = 
  dynamic_cast<const Solver::EulerLinearTransientSolverInfo&>(this->transient_solver_info).getSystemOrder();
  
  // make sure that the number of vectors provided are one less than the system order
  Assert(state_vecs.size() == system_order,
         ExcInternalError());
  
  this->current_state_vectors.resize(system_order+1);
  
  std::auto_ptr<NumericVector<double> > vec;
  
  for (unsigned int i=0; i < system_order; i++)
    {
      vec.reset(NumericVector<double>::build().release());
      vec->duplicate_vector(*(state_vecs[i]), true);
      vec->close();
      this->current_state_vectors.reset(i, vec.release());
    }
  
  // finally, insert the highest order vectors. These vectors are later calculate from 
  // the initial conditions. 
  vec.reset(NumericVector<double>::build().release());
  vec->duplicate_vector(*(state_vecs[0]), false);
  vec->close();
  this->current_state_vectors.reset(system_order, vec.release());
}




double
Solver::EulerLinearTransientSolver::getCurrentTime()
{
  return this->current_time;
}


double
Solver::EulerLinearTransientSolver::getSimulatedTime()
{
  return this->simulated_time;
}


double
Solver::EulerLinearTransientSolver::getCurrentStepSize()
{
  return this->time_step_size;
}




unsigned int
Solver::EulerLinearTransientSolver::getSimulatedIterationNumber()
{
  return this->simulated_iteration_number;
}



unsigned int
Solver::EulerLinearTransientSolver::getCurrentIterationNumber()
{
  return this->simulated_iteration_number + 1;
}




Solver::EulerNonlinearTransientSolver::EulerNonlinearTransientSolver
(const Solver::NonlinearTransientSolverInfo& transient_info,
 const Solver::NonlinearSolverInfo& nonlinear_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::NonlinearTransientSolverBase(transient_info, nonlinear_info, 
                                     linear_info),
simulated_time(0.0),
current_time(0.0),
time_step_size(0.0),
simulated_iteration_number(0)
{
}



Solver::EulerNonlinearTransientSolver::~EulerNonlinearTransientSolver()
{
  
}




void
Solver::EulerNonlinearTransientSolver::clear()
{
  this->simulated_time = 0.0;
  this->current_time = 0.0;
  this->time_step_size = 0.0;
  this->simulated_iteration_number = 0;
  this->last_iter_state_vectors.clear();
  this->current_state_vectors.clear();
  this->new_state_vectors.clear();
  Solver::NonlinearTransientSolverBase::clear();
}





void
Solver::EulerNonlinearTransientSolver::solve()
{
//  this->nonlinear_solver.reset
//  (Solver::createNonlinearSolver(this->nonlinear_solver_info,
//                                 this->linear_solver_info).release());
//  
//  const Solver::EulerNonlinearTransientSolverInfo& newmark_info =
//  dynamic_cast<const Solver::EulerNonlinearTransientSolverInfo&>(this->transient_solver_info);
//  
//  unsigned int system_order = newmark_info.getSystemOrder();
//  
//  // make sure that the solver has the necessary information
//  Assert(this->coefficient_matrices.size() == (system_order+1), ExcInternalError());
//  Assert(this->current_state_vectors.size() == (system_order+1), ExcInternalError());
//  Assert(this->new_state_vectors.size() == (system_order+1), ExcInternalError());
//  Assert(this->jacobian_matrix != NULL, ExcInternalError());
//  
//  this->analysis_driver->getFESystemController().performance_logging->setEvent
//  ("EulerNonlinearTransientSolver::solve()", "Solver");
//  
//  
//  // create the vectors to be used in this solution
//  std::auto_ptr<NumericVector<double> > 
//  res_vec(NumericVector<double>::build().release()),
//  nonlinear_sol_vec(NumericVector<double>::build().release());
//  
//  res_vec->duplicate_vector(*(this->current_state_vectors[0]), false);
//  nonlinear_sol_vec->duplicate_vector(*(this->current_state_vectors[0]), false);
//  
//  // init the constants
//  Solver::EulerSolver::initEulerConstants(newmark_info,
//                                              this->newmark_constants);
//  
//  // init the iteration data 
//  this->simulated_time = this->transient_solver_info.getInitialTime();
//  this->simulated_iteration_number = 0;
//  
//  // calculate the highest order derivative at initial time from the initial conditions
//  // the basic equation to be solved is 
//  // dX/dt^(n) = inv(C_(n)) * (F_0 - sum_(i=0..n-1) (C_i * dX/dt^(i)))
//  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
//  getCoefficientMatrices(this->transient_solver_info.getInitialTime(),
//                         this->current_state_vectors.getReference(),
//                         this->coefficient_matrices);
//  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
//  getForceVector(this->transient_solver_info.getInitialTime(),
//                 this->current_state_vectors.getReference(),
//                 *res_vec);
//  res_vec->close();
//  
//  // calculate the RHS of the linear equation
//  for (unsigned int i=0; i < system_order; i++)
//    {
//      nonlinear_sol_vec->zero(); 
//      nonlinear_sol_vec->close();
//      this->coefficient_matrices[i]->close();
//      this->coefficient_matrices[i]->multiply_vector(*(this->current_state_vectors[i]),
//                                                     *nonlinear_sol_vec);
//      nonlinear_sol_vec->close();
//      res_vec->add(-1.0, *nonlinear_sol_vec);
//    }
//  res_vec->close();
//  
//  // tell the linear solver to factorize this matrix
//  this->coefficient_matrices[system_order]->close();
//  
//  // initialize the linear solver and solve for the derivative
//  std::auto_ptr<Solver::LinearSolver> linear_solver;
//  linear_solver.reset
//  (Solver::createLinearSolver(this->linear_solver_info).release());
//  linear_solver->attachAnalysisDriver(this->analysis_driver);
//  
//  linear_solver->setSystemMatrix(*(this->coefficient_matrices[system_order]),
//                                 Solver::ATTACH_AND_FACTORIZE_MATRIX::num());
//  linear_solver->solve(*res_vec, 
//                       *(this->current_state_vectors[system_order]));
//  this->current_state_vectors[system_order]->close();
//  // now clear the vectors that were used
//  res_vec->zero();
//  res_vec->close();
//  nonlinear_sol_vec->zero();
//  nonlinear_sol_vec->close();
//  // and clear the linear solver
//  linear_solver.reset();
//  
//  
//  // now, increment the iteration number to 1, for which the vectors will be stored
//  this->simulated_iteration_number++;
//  
//  
//  // tell the driver of the 1st iterate, which will be the initial conditions
//  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
//  setNewIteration(this->getSimulatedTime(),
//                  this->getSimulatedIterationNumber(),
//                  this->current_state_vectors.getReference());
//  
//  // also copy the initial conditions to the last iteration vector
//  for (unsigned int i = 0; i <= system_order; i++)
//    {
//      *(this->last_iter_state_vectors[i]) = *(this->current_state_vectors[i]);
//      this->last_iter_state_vectors[i]->close();
//    }
//  
//  
//  this->nonlinear_solver->clear();
//  this->nonlinear_solver->attachAnalysisDriver(this->analysis_driver);
//  this->nonlinear_solver->setFunctions(Solver::EulerSolver::getResidual, 
//                                       Solver::EulerSolver::getJacobian,
//                                       (void *) this);
//  this->nonlinear_solver->attachMatrixAndVector(*jacobian_matrix, *res_vec);
//  
//  
//  
//  while (this->getSimulatedTime() < this->transient_solver_info.getFinalTime())
//    {
//      this->current_time = this->simulated_time + this->transient_solver_info.getInitialTimeStep();
//      
//      nonlinear_sol_vec->zero(); 
//      //    *nonlinear_sol_vec = *(this->last_iter_state_vectors[0]);
//      //    nonlinear_sol_vec->close();
//      this->nonlinear_solver->solve(*nonlinear_sol_vec);
//      
//      // copy the vector
//      *(this->new_state_vectors[0]) = *(this->last_iter_state_vectors[0]);
//      this->new_state_vectors[0]->add(1.0,  *nonlinear_sol_vec);
//      this->new_state_vectors[0]->close();
//      
//      // calculate the new derivatives
//      Solver::EulerSolver::updateStateVectors
//      (newmark_info,
//       this->newmark_constants,
//       this->last_iter_state_vectors.getReference(),
//       this->new_state_vectors.getReference());
//      
//      // increment the current time and the iteration number
//      this->simulated_time = this->getCurrentTime();
//      this->simulated_iteration_number++;
//      
//      // tell the analysis driver to record the iteration increment
//      dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
//      setNewIteration(this->getSimulatedTime(),
//                      this->getSimulatedIterationNumber(),
//                      this->last_iter_state_vectors.getReference());
//    }
//  
//  this->nonlinear_solver.reset();
//  
//  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
//  ("EulerNonlinearTransientSolver::solve()", "Solver"); 
}




void
Solver::EulerSolver::getResidual(NumericVector<double>& sol,
                                   NumericVector<double>& res,
                                   void* ctx)
{
//  Solver::EulerNonlinearTransientSolver* solver = 
//  static_cast<Solver::EulerNonlinearTransientSolver*> (ctx);
//  
//  const Solver::EulerNonlinearTransientSolverInfo& newmark_info =
//  dynamic_cast<const Solver::EulerNonlinearTransientSolverInfo&>(solver->transient_solver_info);
//  
//  unsigned int system_order = newmark_info.getSystemOrder();
//  
//  // copy the previous iteration solutions needef for calculation of new state vectors
//  for (unsigned int i=0; i <= system_order; i++)
//    {
//      *(solver->current_state_vectors[i]) = *(solver->last_iter_state_vectors[i]);
//      solver->current_state_vectors[i]->close();
//    }
//  
//  // now copy the new approximation and calculate the residual
//  sol.close();
//  *(solver->new_state_vectors[0]) = *(solver->last_iter_state_vectors[0]);
//  solver->new_state_vectors[0]->add(1.0, sol);
//  solver->new_state_vectors[0]->close();
//  
//  // calculate the new derivatives
//  Solver::EulerSolver::updateStateVectors
//  (newmark_info,
//   solver->newmark_constants,
//   solver->current_state_vectors.getReference(),
//   solver->new_state_vectors.getReference());
//  
//  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(solver->analysis_driver)->
//  getResidualVector(solver->getCurrentTime(),
//                    solver->current_state_vectors.getReference(),
//                    res);
}




void
Solver::EulerSolver::getJacobian(NumericVector<double>& sol,
                                   SparseMatrix<double>& jac,
                                   void* ctx)
{
//  Solver::EulerNonlinearTransientSolver* solver = 
//  static_cast<Solver::EulerNonlinearTransientSolver*> (ctx);
//  
//  const Solver::EulerNonlinearTransientSolverInfo& newmark_info =
//  dynamic_cast<const Solver::EulerNonlinearTransientSolverInfo&>(solver->transient_solver_info);
//  
//  unsigned int system_order = newmark_info.getSystemOrder();
//  
//  // copy the previous iteration solutions needef for calculation of new state vectors
//  for (unsigned int i=0; i <= system_order; i++)
//    {
//      *(solver->current_state_vectors[i]) = *(solver->last_iter_state_vectors[i]);
//      solver->current_state_vectors[i]->close();
//    }
//  
//  // now copy the new approximation and calculate the residual
//  sol.close();
//  *(solver->new_state_vectors[0]) = *(solver->last_iter_state_vectors[0]);
//  solver->new_state_vectors[0]->add(1.0, sol);
//  solver->new_state_vectors[0]->close();
//  
//  // calculate the new derivatives
//  Solver::EulerSolver::updateStateVectors
//  (newmark_info,
//   solver->newmark_constants,
//   solver->current_state_vectors.getReference(),
//   solver->new_state_vectors.getReference());
//  
//  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(solver->analysis_driver)->
//  getJacobianMatrix(solver->getCurrentTime(),
//                    solver->current_state_vectors.getReference(),
//                    jac);
//  
//  jac.close();
//  
//  // now, get the coefficient matrices to add to the jacobian
//  std::vector<SparseMatrix<double>*> coeff_mats(system_order + 1);
//  
//  switch (system_order)
//  {
//    case 1:
//    {
//      // only the C1 matrix is needed
//      coeff_mats[0] = NULL;
//      coeff_mats[1] = solver->coefficient_matrices[1];
//      
//      dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(solver->analysis_driver)->
//      getCoefficientMatrices(solver->getCurrentTime(),
//                             solver->current_state_vectors.getReference(),
//                             coeff_mats);
//      coeff_mats[1]->close();
//      jac.add(solver->newmark_constants[0], *(coeff_mats[1]));
//      jac.close();
//    }
//      break;
//      
//    default:
//      Assert(false, ExcInternalError());
//  }
}





void
Solver::EulerNonlinearTransientSolver::setInitialCondition
(std::vector<NumericVector<double>*>& state_vecs)
{
//  // make sure that the current solutions have not been initialized
//  Assert(this->current_state_vectors.size() == 0, ExcInternalError());
//  Assert(this->new_state_vectors.size() == 0, ExcInternalError());
//  
//  unsigned int system_order = 
//  dynamic_cast<const Solver::EulerNonlinearTransientSolverInfo&>(this->transient_solver_info).getSystemOrder();
//  
//  // make sure that the number of vectors provided are one less than the system order
//  Assert(state_vecs.size() == system_order,
//         ExcInternalError());
//  
//  this->current_state_vectors.resize(system_order+1);
//  this->last_iter_state_vectors.resize(system_order+1);
//  this->new_state_vectors.resize(system_order+1);
//  
//  std::auto_ptr<NumericVector<double> > vec;
//  
//  for (unsigned int i=0; i < system_order; i++)
//    {
//      vec.reset(NumericVector<double>::build().release());
//      vec->duplicate_vector(*(state_vecs[i]), false);
//      vec->close();
//      this->last_iter_state_vectors.reset(i, vec.release());
//      
//      vec.reset(NumericVector<double>::build().release());
//      vec->duplicate_vector(*(state_vecs[i]), true);
//      vec->close();
//      this->current_state_vectors.reset(i, vec.release());
//      
//      vec.reset(NumericVector<double>::build().release());
//      vec->duplicate_vector(*(state_vecs[i]), false);
//      vec->close();
//      this->new_state_vectors.reset(i, vec.release());
//    }
//  
//  // finally, insert the highest order vectors
//  vec.reset(NumericVector<double>::build().release());
//  vec->duplicate_vector(*(state_vecs[0]), false);
//  vec->close();
//  this->last_iter_state_vectors.reset(system_order, vec.release());
//  
//  vec.reset(NumericVector<double>::build().release());
//  vec->duplicate_vector(*(state_vecs[0]), false);
//  vec->close();
//  this->current_state_vectors.reset(system_order, vec.release());
//  
//  vec.reset(NumericVector<double>::build().release());
//  vec->duplicate_vector(*(state_vecs[0]), false);
//  vec->close();
//  this->new_state_vectors.reset(system_order, vec.release());
}




double
Solver::EulerNonlinearTransientSolver::getCurrentTime()
{
  return this->current_time;
}


double
Solver::EulerNonlinearTransientSolver::getSimulatedTime()
{
  return this->simulated_time;
}


double
Solver::EulerNonlinearTransientSolver::getCurrentStepSize()
{
  return this->time_step_size;
}




unsigned int
Solver::EulerNonlinearTransientSolver::getSimulatedIterationNumber()
{
  return this->simulated_iteration_number;
}



unsigned int
Solver::EulerNonlinearTransientSolver::getCurrentIterationNumber()
{
  return this->simulated_iteration_number + 1;
}

