// $Id: NonlinearSolver.C,v 1.5.4.4 2007-05-14 16:45:07 manav Exp $

// FESystem includes
#include "Solvers/NonlinearSolver.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Utilities/TimeLogs.h"


Solver::NonlinearSolver::NonlinearSolver(const Solver::NonlinearSolverInfo& nonlinear_info,
                                         const Solver::LinearSolverInfo& linear_info):
Solver::FESystemSolverBase(dynamic_cast<const Solver::SolverInfo&>(nonlinear_info)),
nonlinear_solver_kind_enum_ID(nonlinear_info.getNonlinearSolverKindEnumID()),
nonlinear_solver_info(nonlinear_info),
linear_solver_info(linear_info),
jacobian(NULL),
residual(NULL),
jacobian_function(NULL),
residual_function(NULL),
context(NULL)
{
  
}




Solver::NonlinearSolver::~NonlinearSolver()
{
  
}


void
Solver::NonlinearSolver::clear()
{
  // first clear the data structures for this class and then call the 
  // parent class' method
  this->jacobian = NULL;
  this->residual = NULL;
  this->jacobian_function = NULL;
  this->residual_function = NULL;
  this->context = NULL;
  
  Solver::FESystemSolverBase::clear();
}



void
Solver::NonlinearSolver::setFunctions(void (*res_func)(NumericVector<double>& sol,
                                                       NumericVector<double>& res,
                                                       void* ctx),
                                      void (*jac_func)(NumericVector<double>& sol,
                                                       SparseMatrix<double>& jac,
                                                       void* ctx),
                                      void* ctx)
{
  Assert(this->residual_function == NULL, ExcInternalError());
  Assert(this->jacobian_function == NULL, ExcInternalError());
  Assert(this->context == NULL, ExcInternalError());
  Assert(res_func != NULL, ExcInternalError());
  Assert(jac_func != NULL, ExcInternalError());
  
  
  this->residual_function = res_func;
  this->jacobian_function = jac_func;
  this->context = ctx;
}


void 
Solver::NonlinearSolver::evaluateResidual(NumericVector<double>& sol,
                                          NumericVector<double>& res)
{
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("NonlinearSolver::evaluateResidual()", "Solver");
  
  if (this->jacobian_function == NULL)
    {
    Driver::NonLinearAnalysisDriver* driver = 
    dynamic_cast<Driver::NonLinearAnalysisDriver*>(this->analysis_driver);
    
    driver->getResidualAtVector(sol, res);
    }
  else
    this->residual_function(sol,res, this->context);

  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
    ("NonlinearSolver::evaluateResidual()", "Solver");
}




void Solver::NonlinearSolver::evaluateJacobian(NumericVector<double>& sol,
                                               SparseMatrix<double>& jac)
{
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("NonlinearSolver::evaluateJacobian()", "Solver");

  if (this->jacobian_function == NULL)
    {
    Driver::NonLinearAnalysisDriver* driver = 
    dynamic_cast<Driver::NonLinearAnalysisDriver*>(this->analysis_driver);

    driver->getJacobianAtVector(sol, jac);
    }
  else
    this->jacobian_function(sol, jac, this->context);

  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
    ("NonlinearSolver::evaluateJacobian()", "Solver");
}

