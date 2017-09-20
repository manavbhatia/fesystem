// $Id: NewmarkTransientSolver.C,v 1.1.2.5 2008-08-21 20:43:10 manav Exp $

// FESystem includes
#include "Solvers/NewmarkTransientSolver.h"
#include "FESystem/FESystemController.h"
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/NewmarkTransientSolverInfo.h"
#include "Utilities/TimeLogs.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"



Solver::NewmarkLinearTransientSolver::NewmarkLinearTransientSolver
(const Solver::LinearTransientSolverInfo& transient_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::LinearTransientSolverBase(transient_info, linear_info),
simulated_time(0.0),
current_time(0.0),
time_step_size(0.0),
simulated_iteration_number(0)
{
}



Solver::NewmarkLinearTransientSolver::~NewmarkLinearTransientSolver()
{
  
}




void
Solver::NewmarkLinearTransientSolver::clear()
{
  this->simulated_time = 0.0;
  this->current_time = 0.0;
  this->time_step_size = 0.0;
  this->simulated_iteration_number = 0;
  this->current_state_vectors.clear();
  this->new_state_vectors.clear();
  Solver::LinearTransientSolverBase::clear();
}





void
Solver::NewmarkLinearTransientSolver::solve()
{
  this->linear_solver.reset(Solver::createLinearSolver(this->linear_solver_info).release());
  this->linear_solver->attachAnalysisDriver(this->analysis_driver);
  
  const Solver::NewmarkLinearTransientSolverInfo& newmark_info =
  dynamic_cast<const Solver::NewmarkLinearTransientSolverInfo&>(this->transient_solver_info);
  
  unsigned int system_order = newmark_info.getSystemOrder();
  
  // make sure that the solver has the necessary information
  Assert(this->coefficient_matrices.size() == (system_order+1), ExcInternalError());
  Assert(this->current_state_vectors.size() == (system_order+1), ExcInternalError());
  Assert(this->new_state_vectors.size() == (system_order+1), ExcInternalError());
  
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("NewmarkLinearTransientSolver::solve()", "Solver");
  
  
  // create the vectors to be used in this solution
  std::auto_ptr<NumericVector<double> > 
  external_force_vector(NumericVector<double>::build().release()),
  scratch_vec1(NumericVector<double>::build().release()),
  scratch_vec2(NumericVector<double>::build().release()),
  newmark_force_vector(NumericVector<double>::build().release());
  
  external_force_vector->duplicate_vector(*(this->current_state_vectors[0]), false);
  scratch_vec1->duplicate_vector(*(this->current_state_vectors[0]), false);
  scratch_vec2->duplicate_vector(*(this->current_state_vectors[0]), false);
  newmark_force_vector->duplicate_vector(*(this->current_state_vectors[0]), false);
  
  // create the temporary matrix for calculation of the Newmark stiffness matrix
  std::auto_ptr<SparseMatrix<double> > newmark_stiffness_matrix(SparseMatrix<double>::build().release());
  this->coefficient_matrices[0]->close();
  newmark_stiffness_matrix->duplicate_matrix(*(this->coefficient_matrices[0]), false);
  
  
  // init the constants
  Solver::NewmarkSolver::initNewmarkConstants(newmark_info,
                                              this->newmark_constants);
  
  // init the iteration data 
  this->simulated_time = this->transient_solver_info.getInitialTime();
  this->simulated_iteration_number = 0;
  
  // calculate the highest order derivative at initial time from the initial conditions
  // the basic equation to be solved is 
  // dX/dt^(n) = inv(C_(n)) * (F_0 - sum_(i=0..n-1) (C_i * dX/dt^(i)))
  dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
  getCoefficientMatrices(this->transient_solver_info.getInitialTime(),
                         this->coefficient_matrices);
  dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
  getForceVector(this->transient_solver_info.getInitialTime(),
                 *external_force_vector);
  external_force_vector->close();
  
  // calculate the RHS of the linear equation
  for (unsigned int i=0; i < system_order; i++)
    {
      scratch_vec1->zero(); 
      scratch_vec1->close();
      this->coefficient_matrices[i]->close();
      this->coefficient_matrices[i]->multiply_vector(*(this->current_state_vectors[i]),
                                                     *scratch_vec1);
      scratch_vec1->close();
      external_force_vector->add(-1.0, *scratch_vec1);
    }
  external_force_vector->close();
  
  // tell the linear solver to factorize this matrix
  this->coefficient_matrices[system_order]->close();
  this->linear_solver->setSystemMatrix(*(this->coefficient_matrices[system_order]));
  this->linear_solver->setPreconditionerMatrix(*(this->coefficient_matrices[system_order]));
  this->linear_solver->solve(*external_force_vector, 
                             *(this->current_state_vectors[system_order]));
  this->current_state_vectors[system_order]->close();

  
  // now, increment the iteration number to 1, for which the vectors will be stored
  this->simulated_iteration_number++;

  // tell the driver of the 1st iterate, which will be the initial conditions
  dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
  setNewIteration(this->getSimulatedTime(),
                  this->getSimulatedIterationNumber(),
                  this->current_state_vectors.getReference());  

  
  while (this->getSimulatedTime() < this->transient_solver_info.getFinalTime())
    {
      this->current_time = this->simulated_time + this->transient_solver_info.getInitialTimeStep();
      
      // ask the driver to init the matrices. This needs to be done only once for a system where 
      // the matrices are constant with iterations
      if (!this->time_independent_coefficient_matrices || this->simulated_iteration_number == 1)
        {
          dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
          getCoefficientMatrices(this->getCurrentTime(),
                                 this->coefficient_matrices);
          
          // calculate the Newmark stiffness matrix
          Solver::NewmarkSolver::initNewmarkStiffnessMatrix
          (newmark_info,
           this->newmark_constants,
           this->coefficient_matrices,
           *newmark_stiffness_matrix);
          
          // tell the linear solver to factorize this matrix
          this->linear_solver->setSystemMatrix(*newmark_stiffness_matrix);
          this->linear_solver->setPreconditionerMatrix(*newmark_stiffness_matrix);
        }
      
      
      // ask the driver to init the force vector;
      dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
      getForceVector(this->getCurrentTime(),
                     *external_force_vector);
      
      
      // calculate the Newmark force vector
      Solver::NewmarkSolver::initNewmarkForceVector
      (newmark_info,
       this->newmark_constants,
       this->coefficient_matrices,
       this->current_state_vectors.getReference(),
       *external_force_vector,
       *scratch_vec1,
       *scratch_vec2,
       *newmark_force_vector);
      
      
      // now solve for the new iterate
      this->linear_solver->solve(*newmark_force_vector, *(this->new_state_vectors[0]));
      
      // calculate the new derivatives
      Solver::NewmarkSolver::updateStateVectors
      (newmark_info,
       this->newmark_constants,
       this->current_state_vectors.getReference(),
       this->new_state_vectors.getReference());
      
      // increment the current time and the iteration number
      this->simulated_time = this->getCurrentTime();
      this->simulated_iteration_number++;
      
      // tell the analysis driver to record the iteration increment
      dynamic_cast<Driver::TransientAnalysisDriver*>(this->analysis_driver)->
      setNewIteration(this->getSimulatedTime(),
                      this->getSimulatedIterationNumber(),
                      this->current_state_vectors.getReference());
      
      
      
    }
  
  this->linear_solver.reset();
  
  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
  ("NewmarkLinearTransientSolver::solve()", "Solver"); 
}



void
Solver::NewmarkLinearTransientSolver::setInitialCondition
(std::vector<NumericVector<double>*>& state_vecs)
{
  // make sure that the current solutions have not been initialized
  Assert(this->current_state_vectors.size() == 0, ExcInternalError());
  Assert(this->new_state_vectors.size() == 0, ExcInternalError());
  
  unsigned int system_order = 
  dynamic_cast<const Solver::NewmarkLinearTransientSolverInfo&>(this->transient_solver_info).getSystemOrder();
  
  // make sure that the number of vectors provided are one less than the system order
  Assert(state_vecs.size() == system_order,
         ExcInternalError());
  
  this->current_state_vectors.resize(system_order+1);
  this->new_state_vectors.resize(system_order+1);
  
  std::auto_ptr<NumericVector<double> > vec;
  
  for (unsigned int i=0; i < system_order; i++)
    {
      vec.reset(NumericVector<double>::build().release());
      vec->duplicate_vector(*(state_vecs[i]), true);
      vec->close();
      this->current_state_vectors.reset(i, vec.release());
      
      vec.reset(NumericVector<double>::build().release());
      vec->duplicate_vector(*(state_vecs[i]), false);
      vec->close();
      this->new_state_vectors.reset(i, vec.release());
    }
  
  // finally, insert the highest order vectors. These vectors are later calculate from 
  // the initial conditions. 
  vec.reset(NumericVector<double>::build().release());
  vec->duplicate_vector(*(state_vecs[0]), false);
  vec->close();
  this->current_state_vectors.reset(system_order, vec.release());
  
  vec.reset(NumericVector<double>::build().release());
  vec->duplicate_vector(*(state_vecs[0]), false);
  vec->close();
  this->new_state_vectors.reset(system_order, vec.release());
}




double
Solver::NewmarkLinearTransientSolver::getCurrentTime()
{
  return this->current_time;
}


double
Solver::NewmarkLinearTransientSolver::getSimulatedTime()
{
  return this->simulated_time;
}


double
Solver::NewmarkLinearTransientSolver::getCurrentStepSize()
{
  return this->time_step_size;
}




unsigned int
Solver::NewmarkLinearTransientSolver::getSimulatedIterationNumber()
{
  return this->simulated_iteration_number;
}



unsigned int
Solver::NewmarkLinearTransientSolver::getCurrentIterationNumber()
{
  return this->simulated_iteration_number + 1;
}




Solver::NewmarkNonlinearTransientSolver::NewmarkNonlinearTransientSolver
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



Solver::NewmarkNonlinearTransientSolver::~NewmarkNonlinearTransientSolver()
{
  
}




void
Solver::NewmarkNonlinearTransientSolver::clear()
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
Solver::NewmarkNonlinearTransientSolver::solve()
{
  this->nonlinear_solver.reset
  (Solver::createNonlinearSolver(this->nonlinear_solver_info,
                                 this->linear_solver_info).release());
  
  const Solver::NewmarkNonlinearTransientSolverInfo& newmark_info =
  dynamic_cast<const Solver::NewmarkNonlinearTransientSolverInfo&>(this->transient_solver_info);
  
  unsigned int system_order = newmark_info.getSystemOrder();
  
  // make sure that the solver has the necessary information
  Assert(this->coefficient_matrices.size() == (system_order+1), ExcInternalError());
  Assert(this->current_state_vectors.size() == (system_order+1), ExcInternalError());
  Assert(this->new_state_vectors.size() == (system_order+1), ExcInternalError());
  Assert(this->jacobian_matrix != NULL, ExcInternalError());
  
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("NewmarkNonlinearTransientSolver::solve()", "Solver");
  
  
  // create the vectors to be used in this solution
  std::auto_ptr<NumericVector<double> > 
  res_vec(NumericVector<double>::build().release()),
  nonlinear_sol_vec(NumericVector<double>::build().release());
  
  res_vec->duplicate_vector(*(this->current_state_vectors[0]), false);
  nonlinear_sol_vec->duplicate_vector(*(this->current_state_vectors[0]), false);
  
  // init the constants
  Solver::NewmarkSolver::initNewmarkConstants(newmark_info,
                                              this->newmark_constants);
  
  // init the iteration data 
  this->simulated_time = this->transient_solver_info.getInitialTime();
  this->simulated_iteration_number = 0;
  
  // calculate the highest order derivative at initial time from the initial conditions
  // the basic equation to be solved is 
  // dX/dt^(n) = inv(C_(n)) * (F_0 - sum_(i=0..n-1) (C_i * dX/dt^(i)))
  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
  getCoefficientMatrices(this->transient_solver_info.getInitialTime(),
                         this->current_state_vectors.getReference(),
                         this->coefficient_matrices);
  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
  getForceVector(this->transient_solver_info.getInitialTime(),
                 this->current_state_vectors.getReference(),
                 *res_vec);
  res_vec->close();
  
  // calculate the RHS of the linear equation
  for (unsigned int i=0; i < system_order; i++)
    {
      nonlinear_sol_vec->zero(); 
      nonlinear_sol_vec->close();
      this->coefficient_matrices[i]->close();
      this->coefficient_matrices[i]->multiply_vector(*(this->current_state_vectors[i]),
                                                     *nonlinear_sol_vec);
      nonlinear_sol_vec->close();
      res_vec->add(-1.0, *nonlinear_sol_vec);
    }
  res_vec->close();
  
  // tell the linear solver to factorize this matrix
  this->coefficient_matrices[system_order]->close();
  
  // initialize the linear solver and solve for the derivative
  std::auto_ptr<Solver::LinearSolver> linear_solver;
  linear_solver.reset
  (Solver::createLinearSolver(this->linear_solver_info).release());
  linear_solver->attachAnalysisDriver(this->analysis_driver);
  
  linear_solver->setSystemMatrix(*(this->coefficient_matrices[system_order]));
  linear_solver->setPreconditionerMatrix(*(this->coefficient_matrices[system_order]));
  linear_solver->solve(*res_vec, 
                       *(this->current_state_vectors[system_order]));
  this->current_state_vectors[system_order]->close();
  // now clear the vectors that were used
  res_vec->zero();
  res_vec->close();
  nonlinear_sol_vec->zero();
  nonlinear_sol_vec->close();
  // and clear the linear solver
  linear_solver.reset();
  
  
  // now, increment the iteration number to 1, for which the vectors will be stored
  this->simulated_iteration_number++;
  
  
  // tell the driver of the 1st iterate, which will be the initial conditions
  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
  setNewIteration(this->getSimulatedTime(),
                  this->getSimulatedIterationNumber(),
                  this->current_state_vectors.getReference());
  
  // also copy the initial conditions to the last iteration vector
  for (unsigned int i = 0; i <= system_order; i++)
    {
      *(this->last_iter_state_vectors[i]) = *(this->current_state_vectors[i]);
      this->last_iter_state_vectors[i]->close();
    }
  
  
  this->nonlinear_solver->clear();
  this->nonlinear_solver->attachAnalysisDriver(this->analysis_driver);
  this->nonlinear_solver->setFunctions(Solver::NewmarkSolver::getResidual, 
                                       Solver::NewmarkSolver::getJacobian,
                                       (void *) this);
  this->nonlinear_solver->attachMatrixAndVector(*jacobian_matrix, *res_vec);
  
  
  
  while (this->getSimulatedTime() < this->transient_solver_info.getFinalTime())
    {
      this->current_time = this->simulated_time + this->transient_solver_info.getInitialTimeStep();
      
      nonlinear_sol_vec->zero(); 
      //    *nonlinear_sol_vec = *(this->last_iter_state_vectors[0]);
      //    nonlinear_sol_vec->close();
      this->nonlinear_solver->solve(*nonlinear_sol_vec);
      
      // copy the vector
      *(this->new_state_vectors[0]) = *(this->last_iter_state_vectors[0]);
      this->new_state_vectors[0]->add(1.0,  *nonlinear_sol_vec);
      this->new_state_vectors[0]->close();
      
      // calculate the new derivatives
      Solver::NewmarkSolver::updateStateVectors
      (newmark_info,
       this->newmark_constants,
       this->last_iter_state_vectors.getReference(),
       this->new_state_vectors.getReference());
      
      // increment the current time and the iteration number
      this->simulated_time = this->getCurrentTime();
      this->simulated_iteration_number++;
      
      // tell the analysis driver to record the iteration increment
      dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver)->
      setNewIteration(this->getSimulatedTime(),
                      this->getSimulatedIterationNumber(),
                      this->last_iter_state_vectors.getReference());
    }
  
  this->nonlinear_solver.reset();
  
  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
  ("NewmarkNonlinearTransientSolver::solve()", "Solver"); 
}




void
Solver::NewmarkSolver::getResidual(NumericVector<double>& sol,
                                   NumericVector<double>& res,
                                   void* ctx)
{
  Solver::NewmarkNonlinearTransientSolver* solver = 
  static_cast<Solver::NewmarkNonlinearTransientSolver*> (ctx);
  
  const Solver::NewmarkNonlinearTransientSolverInfo& newmark_info =
  dynamic_cast<const Solver::NewmarkNonlinearTransientSolverInfo&>(solver->transient_solver_info);
  
  unsigned int system_order = newmark_info.getSystemOrder();
  
  // copy the previous iteration solutions needef for calculation of new state vectors
  for (unsigned int i=0; i <= system_order; i++)
    {
      *(solver->current_state_vectors[i]) = *(solver->last_iter_state_vectors[i]);
      solver->current_state_vectors[i]->close();
    }
  
  // now copy the new approximation and calculate the residual
  sol.close();
  *(solver->new_state_vectors[0]) = *(solver->last_iter_state_vectors[0]);
  solver->new_state_vectors[0]->add(1.0, sol);
  solver->new_state_vectors[0]->close();
  
  // calculate the new derivatives
  Solver::NewmarkSolver::updateStateVectors
  (newmark_info,
   solver->newmark_constants,
   solver->current_state_vectors.getReference(),
   solver->new_state_vectors.getReference());
  
  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(solver->analysis_driver)->
  getResidualVector(solver->getCurrentTime(),
                    solver->current_state_vectors.getReference(),
                    res);
}




void
Solver::NewmarkSolver::getJacobian(NumericVector<double>& sol,
                                   SparseMatrix<double>& jac,
                                   void* ctx)
{
  Solver::NewmarkNonlinearTransientSolver* solver = 
  static_cast<Solver::NewmarkNonlinearTransientSolver*> (ctx);
  
  const Solver::NewmarkNonlinearTransientSolverInfo& newmark_info =
  dynamic_cast<const Solver::NewmarkNonlinearTransientSolverInfo&>(solver->transient_solver_info);
  
  unsigned int system_order = newmark_info.getSystemOrder();
  
  // copy the previous iteration solutions needef for calculation of new state vectors
  for (unsigned int i=0; i <= system_order; i++)
    {
      *(solver->current_state_vectors[i]) = *(solver->last_iter_state_vectors[i]);
      solver->current_state_vectors[i]->close();
    }
  
  // now copy the new approximation and calculate the residual
  sol.close();
  *(solver->new_state_vectors[0]) = *(solver->last_iter_state_vectors[0]);
  solver->new_state_vectors[0]->add(1.0, sol);
  solver->new_state_vectors[0]->close();
  
  // calculate the new derivatives
  Solver::NewmarkSolver::updateStateVectors
  (newmark_info,
   solver->newmark_constants,
   solver->current_state_vectors.getReference(),
   solver->new_state_vectors.getReference());
  
  dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(solver->analysis_driver)->
  getJacobianMatrix(solver->getCurrentTime(),
                    solver->current_state_vectors.getReference(),
                    jac);
  
  jac.close();
  
  // now, get the coefficient matrices to add to the jacobian
  std::vector<SparseMatrix<double>*> coeff_mats(system_order + 1);
  
  switch (system_order)
  {
    case 1:
    {
      // only the C1 matrix is needed
      coeff_mats[0] = NULL;
      coeff_mats[1] = solver->coefficient_matrices[1];
      
      dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(solver->analysis_driver)->
      getCoefficientMatrices(solver->getCurrentTime(),
                             solver->current_state_vectors.getReference(),
                             coeff_mats);
      coeff_mats[1]->close();
      jac.add(solver->newmark_constants[0], *(coeff_mats[1]));
      jac.close();
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
}





void
Solver::NewmarkNonlinearTransientSolver::setInitialCondition
(std::vector<NumericVector<double>*>& state_vecs)
{
  // make sure that the current solutions have not been initialized
  Assert(this->current_state_vectors.size() == 0, ExcInternalError());
  Assert(this->new_state_vectors.size() == 0, ExcInternalError());
  
  unsigned int system_order = 
  dynamic_cast<const Solver::NewmarkNonlinearTransientSolverInfo&>(this->transient_solver_info).getSystemOrder();
  
  // make sure that the number of vectors provided are one less than the system order
  Assert(state_vecs.size() == system_order,
         ExcInternalError());
  
  this->current_state_vectors.resize(system_order+1);
  this->last_iter_state_vectors.resize(system_order+1);
  this->new_state_vectors.resize(system_order+1);
  
  std::auto_ptr<NumericVector<double> > vec;
  
  for (unsigned int i=0; i < system_order; i++)
    {
      vec.reset(NumericVector<double>::build().release());
      vec->duplicate_vector(*(state_vecs[i]), false);
      vec->close();
      this->last_iter_state_vectors.reset(i, vec.release());
      
      vec.reset(NumericVector<double>::build().release());
      vec->duplicate_vector(*(state_vecs[i]), true);
      vec->close();
      this->current_state_vectors.reset(i, vec.release());
      
      vec.reset(NumericVector<double>::build().release());
      vec->duplicate_vector(*(state_vecs[i]), false);
      vec->close();
      this->new_state_vectors.reset(i, vec.release());
    }
  
  // finally, insert the highest order vectors
  vec.reset(NumericVector<double>::build().release());
  vec->duplicate_vector(*(state_vecs[0]), false);
  vec->close();
  this->last_iter_state_vectors.reset(system_order, vec.release());
  
  vec.reset(NumericVector<double>::build().release());
  vec->duplicate_vector(*(state_vecs[0]), false);
  vec->close();
  this->current_state_vectors.reset(system_order, vec.release());
  
  vec.reset(NumericVector<double>::build().release());
  vec->duplicate_vector(*(state_vecs[0]), false);
  vec->close();
  this->new_state_vectors.reset(system_order, vec.release());
}




double
Solver::NewmarkNonlinearTransientSolver::getCurrentTime()
{
  return this->current_time;
}


double
Solver::NewmarkNonlinearTransientSolver::getSimulatedTime()
{
  return this->simulated_time;
}


double
Solver::NewmarkNonlinearTransientSolver::getCurrentStepSize()
{
  return this->time_step_size;
}




unsigned int
Solver::NewmarkNonlinearTransientSolver::getSimulatedIterationNumber()
{
  return this->simulated_iteration_number;
}



unsigned int
Solver::NewmarkNonlinearTransientSolver::getCurrentIterationNumber()
{
  return this->simulated_iteration_number + 1;
}


template <typename InfoType>
void
Solver::NewmarkSolver::initNewmarkConstants(const InfoType& solver_info,
                                            std::vector<double>& constants_vector)
{
  // this should be cleared before the method is called
  Assert(constants_vector.size() == 0, ExcInternalError());
  
  const double dt = solver_info.getInitialTimeStep();
  
  switch (solver_info.getSystemOrder())
  {
    case 0:
    {
      // needs to be implemented
      Assert(false, ExcInternalError());
    }
      break;
      
    case 1:
    {
      const double theta = solver_info.getConstant(0);
      
      constants_vector.resize(2);
      constants_vector[0] = 1.0/ (theta * dt);
      constants_vector[1] = (1.0/theta) - 1.0;
    }
      break;
      
    case 2:
    {
      const double beta = solver_info.getConstant(0);
      const double gamma = solver_info.getConstant(1);
      
      constants_vector.resize(8);
      constants_vector[0] = 1.0/ (beta * dt * dt);
      constants_vector[1] = gamma / (beta * dt);
      constants_vector[2] = 1.0 / (beta * dt);
      constants_vector[3] = 0.5 / beta - 1.0;
      constants_vector[4] = gamma / beta - 1.0;
      constants_vector[5] = 0.5 * dt * (gamma/beta - 2.0);
      constants_vector[6] = dt * (1.0 - gamma);
      constants_vector[7] = dt * gamma;
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
  }
}


template <typename InfoType>
void
Solver::NewmarkSolver::initNewmarkStiffnessMatrix
(const InfoType& solver_info,
 std::vector<double>& constants_vector,
 std::vector<SparseMatrix<double>*>& coeff_matrices,
 SparseMatrix<double>& newmark_stiff_matrix)
{
  Assert(constants_vector.size() != 0, ExcInternalError());
  Assert(coeff_matrices.size() == (solver_info.getSystemOrder()+1), ExcInternalError());
  newmark_stiff_matrix.zero();
  newmark_stiff_matrix.perform_intermediate_matrix_assembly();
  
  switch (solver_info.getSystemOrder())
  {
    case 0:
    {
      // needs to be implemented
      Assert(false, ExcInternalError());
    }
      break;
      
    case 1:
    {
      if (coeff_matrices[0] != NULL)
        newmark_stiff_matrix.add(1.0, *coeff_matrices[0]);
      
      if (coeff_matrices[1] != NULL)
        newmark_stiff_matrix.add(constants_vector[0], *coeff_matrices[1]);
    }
      break;
      
    case 2:
    {
      if (coeff_matrices[0] != NULL)
        newmark_stiff_matrix.add(1.0, *coeff_matrices[0]);
      
      if (coeff_matrices[1] != NULL)
        newmark_stiff_matrix.add(constants_vector[0], *coeff_matrices[1]);
      
      if (coeff_matrices[2] != NULL)
        newmark_stiff_matrix.add(constants_vector[1], *coeff_matrices[2]);
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
  }
  
  newmark_stiff_matrix.close();
}





template <typename InfoType>
void
Solver::NewmarkSolver::initNewmarkForceVector
(const InfoType& solver_info,
 std::vector<double>& constant_vector,
 std::vector<SparseMatrix<double>*>& coeff_matrices,
 std::vector<NumericVector<double>*>& state_vectors,
 NumericVector<double>& external_force_vec,
 NumericVector<double>& scratch_vec1,
 NumericVector<double>& scratch_vec2,
 NumericVector<double>& newmark_force_vec)
{
  Assert(constant_vector.size() != 0, ExcInternalError());
  Assert(coeff_matrices.size() == (solver_info.getSystemOrder()+1), ExcInternalError());
  Assert(state_vectors.size() == (solver_info.getSystemOrder()+1), ExcInternalError());
  newmark_force_vec.zero();
  newmark_force_vec.perform_intermediate_vector_assembly();
  
  switch (solver_info.getSystemOrder())
  {
    case 0:
    {
      // needs to be implemented
      Assert(false, ExcInternalError());
    }
      break;
      
    case 1:
    {
      scratch_vec1.zero();
      scratch_vec1.perform_intermediate_vector_assembly();
      scratch_vec2.zero();
      scratch_vec2.perform_intermediate_vector_assembly();
      
      if (coeff_matrices[1] != NULL)
        {
          scratch_vec1.add(constant_vector[0], *state_vectors[0]);
          scratch_vec1.add(constant_vector[1], *state_vectors[1]);
          scratch_vec1.close();
          coeff_matrices[1]->multiply_vector(scratch_vec1, scratch_vec2);
          scratch_vec2.close();
          newmark_force_vec.add(1.0, scratch_vec2);
        }
    }
      break;
      
    case 2:
    {
      scratch_vec1.zero();
      scratch_vec1.perform_intermediate_vector_assembly();
      scratch_vec2.zero();
      scratch_vec2.perform_intermediate_vector_assembly();
      
      if (coeff_matrices[1] != NULL)
        {
          scratch_vec1.add(constant_vector[1], *state_vectors[0]);
          scratch_vec1.add(constant_vector[4], *state_vectors[1]);
          scratch_vec1.add(constant_vector[5], *state_vectors[2]);
          scratch_vec1.close();
          coeff_matrices[1]->multiply_vector(scratch_vec1, scratch_vec2);
          scratch_vec2.close();
          newmark_force_vec.add(1.0, scratch_vec2);
        }
      
      scratch_vec1.zero();
      scratch_vec1.perform_intermediate_vector_assembly();
      scratch_vec2.zero();
      scratch_vec2.perform_intermediate_vector_assembly();
      
      if (coeff_matrices[2] != NULL)
        {
          scratch_vec1.add(constant_vector[0], *state_vectors[0]);
          scratch_vec1.add(constant_vector[2], *state_vectors[1]);
          scratch_vec1.add(constant_vector[3], *state_vectors[2]);
          scratch_vec1.close();
          coeff_matrices[2]->multiply_vector(scratch_vec1, scratch_vec2);
          scratch_vec2.close();
          newmark_force_vec.add(1.0, scratch_vec2);
        }
      
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
  }
  
  external_force_vec.close();
  newmark_force_vec.add(1.0, external_force_vec);
  newmark_force_vec.close();
}




template <typename InfoType>
void
Solver::NewmarkSolver::updateStateVectors
(const InfoType& solver_info,
 std::vector<double>& constant_vector,
 std::vector<NumericVector<double>*>& previous_state_vector,
 std::vector<NumericVector<double>*>& new_state_vector)
{
  Assert(constant_vector.size() != 0, ExcInternalError());
  Assert(previous_state_vector.size() == (solver_info.getSystemOrder()+1),
         ExcInternalError());
  Assert(previous_state_vector.size() == new_state_vector.size(), ExcInternalError());
  
  switch (solver_info.getSystemOrder())
  {
    case 0:
    {
      // needs to be implemented
      Assert(false, ExcInternalError());
    }
      break;
      
    case 1:
    {
      new_state_vector[0]->close();
      new_state_vector[1]->zero();
      new_state_vector[1]->add(constant_vector[0], *new_state_vector[0]);
      new_state_vector[1]->add(-constant_vector[0], *previous_state_vector[0]);
      new_state_vector[1]->add(-constant_vector[1], *previous_state_vector[1]);
      new_state_vector[1]->close();
      
      // now copy the new vector into the previous state vector for the next iteration
      previous_state_vector[0]->zero();
      *(previous_state_vector[0]) = *(new_state_vector[0]);
      previous_state_vector[0]->close();
      new_state_vector[0]->zero();
      new_state_vector[0]->close();
      
      previous_state_vector[1]->zero();
      *(previous_state_vector[1]) = *(new_state_vector[1]);
      previous_state_vector[1]->close();
      new_state_vector[1]->zero();
      new_state_vector[1]->close();
    }
      break;
      
    case 2:
    {
      new_state_vector[0]->close();
      
      new_state_vector[1]->zero();
      new_state_vector[1]->add(constant_vector[1], *new_state_vector[0]);
      new_state_vector[1]->add(-constant_vector[1], *previous_state_vector[0]);
      new_state_vector[1]->add(-constant_vector[4], *previous_state_vector[1]);
      new_state_vector[1]->add(-constant_vector[5], *previous_state_vector[2]);
      new_state_vector[1]->close();
      
      new_state_vector[2]->zero();
      new_state_vector[2]->add(constant_vector[0], *new_state_vector[0]);
      new_state_vector[2]->add(-constant_vector[0], *previous_state_vector[0]);
      new_state_vector[2]->add(-constant_vector[2], *previous_state_vector[1]);
      new_state_vector[2]->add(-constant_vector[3], *previous_state_vector[2]);
      new_state_vector[2]->close();
      
      // now copy the new vector into the previous state vector for the next iteration
      previous_state_vector[0]->zero();
      *(previous_state_vector[0]) = *(new_state_vector[0]);
      previous_state_vector[0]->close();
      new_state_vector[0]->zero();
      new_state_vector[0]->close();
      
      previous_state_vector[1]->zero();
      *(previous_state_vector[1]) = *(new_state_vector[1]);
      previous_state_vector[1]->close();
      new_state_vector[1]->zero();
      new_state_vector[1]->close();
      
      previous_state_vector[2]->zero();
      *(previous_state_vector[2]) = *(new_state_vector[2]);
      previous_state_vector[2]->close();
      new_state_vector[2]->zero();
      new_state_vector[2]->close();
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
  }
}


template 
void Solver::NewmarkSolver::initNewmarkConstants<Solver::NewmarkLinearTransientSolverInfo>
(const Solver::NewmarkLinearTransientSolverInfo& solver_info,
 std::vector<double>& constants_vector);


template 
void Solver::NewmarkSolver::initNewmarkStiffnessMatrix<Solver::NewmarkLinearTransientSolverInfo>
(const Solver::NewmarkLinearTransientSolverInfo& solver_info,
 std::vector<double>& constants_vector,
 std::vector<SparseMatrix<double>*>& coeff_matrices,
 SparseMatrix<double>& newmark_stiff_matrix);


template 
void Solver::NewmarkSolver::initNewmarkForceVector<Solver::NewmarkLinearTransientSolverInfo>
(const Solver::NewmarkLinearTransientSolverInfo& solver_info,
 std::vector<double>& constant_vector,
 std::vector<SparseMatrix<double>*>& coeff_matrices,
 std::vector<NumericVector<double>*>& state_vectors,
 NumericVector<double>& external_force_vec,
 NumericVector<double>& scratch_vec1,
 NumericVector<double>& scratch_vec2,
 NumericVector<double>& newmark_force_vec);



template 
void Solver::NewmarkSolver::updateStateVectors<Solver::NewmarkLinearTransientSolverInfo>
(const Solver::NewmarkLinearTransientSolverInfo& solver_info,
 std::vector<double>& constant_vector,
 std::vector<NumericVector<double>*>& previous_state_vector,
 std::vector<NumericVector<double>*>& new_state_vector);


template 
void Solver::NewmarkSolver::initNewmarkConstants<Solver::NewmarkNonlinearTransientSolverInfo>
(const Solver::NewmarkNonlinearTransientSolverInfo& solver_info,
 std::vector<double>& constants_vector);


template 
void Solver::NewmarkSolver::initNewmarkStiffnessMatrix<Solver::NewmarkNonlinearTransientSolverInfo>
(const Solver::NewmarkNonlinearTransientSolverInfo& solver_info,
 std::vector<double>& constants_vector,
 std::vector<SparseMatrix<double>*>& coeff_matrices,
 SparseMatrix<double>& newmark_stiff_matrix);


template 
void Solver::NewmarkSolver::initNewmarkForceVector<Solver::NewmarkNonlinearTransientSolverInfo>
(const Solver::NewmarkNonlinearTransientSolverInfo& solver_info,
 std::vector<double>& constant_vector,
 std::vector<SparseMatrix<double>*>& coeff_matrices,
 std::vector<NumericVector<double>*>& state_vectors,
 NumericVector<double>& external_force_vec,
 NumericVector<double>& scratch_vec1,
 NumericVector<double>& scratch_vec2,
 NumericVector<double>& newmark_force_vec);



template 
void Solver::NewmarkSolver::updateStateVectors<Solver::NewmarkNonlinearTransientSolverInfo>
(const Solver::NewmarkNonlinearTransientSolverInfo& solver_info,
 std::vector<double>& constant_vector,
 std::vector<NumericVector<double>*>& previous_state_vector,
 std::vector<NumericVector<double>*>& new_state_vector);

