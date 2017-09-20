// $Id: EigenSolver.C,v 1.13.4.1 2007-03-04 03:37:12 manav Exp $


// FESystem includes
#include "Solvers/EigenSolver.h"
#include "Solvers/EigenSolverInfo.h"
#include "Solvers/LinearSolverInfo.h"



Solver::EigenSolverBase::EigenSolverBase(const Solver::EigenSolverInfo& eigen_info,
                                         const Solver::LinearSolverInfo& linear_info):
Solver::FESystemSolverBase(dynamic_cast<const Solver::SolverInfo&>(eigen_info)),
eigen_solver_info(eigen_info),
linear_solver_info(linear_info)
{
  
}

Solver::EigenSolverBase::~EigenSolverBase()
{
  
}
