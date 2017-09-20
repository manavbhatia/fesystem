// $Id: LinearSolver.C,v 1.5.6.2 2007-05-11 05:16:54 manav Exp $

// FESystem includes
#include "Solvers/LinearSolver.h"
#include "Solvers/LinearSolverInfo.h"


Solver::LinearSolver::LinearSolver(const Solver::LinearSolverInfo& info,
                                   const unsigned int solver_kind):
Solver::FESystemSolverBase(dynamic_cast<const Solver::SolverInfo&>(info)),
linear_solver_kind_enum_ID(solver_kind),
solver_info(info),
ksp_type_enum_ID(info.getKSPTypeEnumID()),
pc_type_enum_ID(info.getPCTypeEnumID()),
system_matrix(NULL),
new_system_matrix(false),
preconditioner_matrix(NULL),
new_PC_matrix(false),
same_matrices_between_solves(false),
n_solves_after_new_matrices(0)
{
}




Solver::LinearSolver::~LinearSolver()
{
  
}



bool
Solver::LinearSolver::useSameMatricesBetweenSolves()
{
  this->same_matrices_between_solves = true;
}




void 
Solver::LinearSolver::clear()
{
  this->system_matrix = NULL;
  this->preconditioner_matrix = NULL;
  this->new_system_matrix = false;
  this->new_PC_matrix = false;
  this->same_matrices_between_solves = false;
  this->n_solves_after_new_matrices = 0;

  // nothing to be done, just call the parent's method
  Solver::FESystemSolverBase::clear();
}


void
Solver::LinearSolver::setSystemMatrix(SparseMatrix<double>& matrix)
{
  this->system_matrix = &matrix;
  this->new_system_matrix = true;
  this->n_solves_after_new_matrices = 0;
}



void 
Solver::LinearSolver::setPreconditionerMatrix(SparseMatrix<double>& matrix)
{
  this->preconditioner_matrix = &matrix;
  this->new_PC_matrix = true;
  this->n_solves_after_new_matrices = 0;
}

