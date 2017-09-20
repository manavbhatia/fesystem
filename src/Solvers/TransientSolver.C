
// FESystem includes
#include "Solvers/TransientSolver.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "Solvers/TransientSolverInfo.h"
#include "Utilities/TimeLogs.h"


Solver::TransientSolverBase::TransientSolverBase
(const Solver::TransientSolverInfoBase& info):
Solver::FESystemSolverBase(dynamic_cast<const Solver::SolverInfo&>(info)),
transient_solver_info(info)
{
  
}



Solver::TransientSolverBase::~TransientSolverBase()
{
  
}


void
Solver::TransientSolverBase::clear()
{
  this->coefficient_matrices.clear();
  
  Solver::FESystemSolverBase::clear();
}



void
Solver::TransientSolverBase::attachCoefficientMatrices
(std::vector<SparseMatrix<double>* > & matrices)
{
  Assert(this->coefficient_matrices.size() == 0, ExcInternalError());
  Assert(matrices.size() != 0, ExcInternalError());

  this->coefficient_matrices = matrices;
}


double 
Solver::TransientSolverBase::getTimeStepSize() const
{
  // Note that this will change when adaptive time steps are added
  return this->transient_solver_info.getInitialTimeStep();
}


Solver::LinearTransientSolverBase::LinearTransientSolverBase
(const Solver::LinearTransientSolverInfo& transient_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::TransientSolverBase(dynamic_cast<const Solver::TransientSolverInfoBase&>(transient_info)),
linear_solver_info(linear_info),
time_independent_coefficient_matrices(false)
{

}


Solver::LinearTransientSolverBase::~LinearTransientSolverBase()
{
  
}



void
Solver::LinearTransientSolverBase::clear()
{
  this->time_independent_coefficient_matrices = false;
  Solver::TransientSolverBase::clear();
}


void
Solver::LinearTransientSolverBase::setCoefficientMatrixTimeDependence(const bool dependence)
{
  this->time_independent_coefficient_matrices = dependence;
}



Solver::NonlinearTransientSolverBase::NonlinearTransientSolverBase
(const Solver::NonlinearTransientSolverInfo& transient_info,
 const Solver::NonlinearSolverInfo& nonlinear_info,
 const Solver::LinearSolverInfo& linear_info):
Solver::TransientSolverBase(dynamic_cast<const Solver::TransientSolverInfoBase&>(transient_info)),
nonlinear_solver_info(nonlinear_info),
linear_solver_info(linear_info),
jacobian_matrix(NULL)
{
}


Solver::NonlinearTransientSolverBase::~NonlinearTransientSolverBase()
{
  
}


void
Solver::NonlinearTransientSolverBase::clear()
{
  this->jacobian_matrix = NULL;
  Solver::TransientSolverBase::clear();
}



void
Solver::NonlinearTransientSolverBase::attachJacobianMatrix(SparseMatrix<double>& jac_mat)
{
  this->jacobian_matrix = &jac_mat;
}



// void
// Solver::NonlinearTransientSolverBase::evaluateRHSFunction(const double current_time,
//                                                           NumericVector<double>& sol,
//                                                           NumericVector<double>& rhs_func)
// {
//   this->analysis_driver->getFESystemController().performance_logging->setEvent
//   ("NonlinearTransientSolver::evaluateRHSFunction()", "Solver");
  
//   Driver::NonlinearTransientAnalysisDriver* driver = 
//   dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver);
  
//   driver->getRHSFunction(current_time, sol, rhs_func);
  
//   this->analysis_driver->getFESystemController().performance_logging->unsetEvent
//     ("NonlinearSolver::evaluateRHSFunction()", "Solver");  
// }



// void 
// Solver::NonlinearTransientSolverBase::evaluateRHSJacobian(const double current_time,
//                                                           NumericVector<double>& sol,
//                                                           SparseMatrix<double>& jac)
// {
//   this->analysis_driver->getFESystemController().performance_logging->setEvent
//   ("NonlinearTransientSolver::evaluateRHSJacobian()", "Solver");
  
//   Driver::NonlinearTransientAnalysisDriver* driver = 
//   dynamic_cast<Driver::NonlinearTransientAnalysisDriver*>(this->analysis_driver);
  
//   driver->getRHSJacobian(current_time, sol, jac);
  
//   this->analysis_driver->getFESystemController().performance_logging->unsetEvent
//     ("NonlinearSolver::evaluateRHSJacobian()", "Solver");
// }



