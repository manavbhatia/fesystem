// $Id: NonlinearSolverInfo.C,v 1.5 2006-09-05 20:41:35 manav Exp $

// FESystem includes
#include "Solvers/NonlinearSolverInfo.h"
#include "Utilities/InputOutputUtility.h"




Solver::NonlinearSolverInfo::NonlinearSolverInfo():
Solver::SolverInfo(NONLINEAR_SOLVER_INFO::num(),
                   NONLINEAR_SOLVER::num()),
nonlinear_solver_kind_enum_ID(FESystemNumbers::InvalidID),
linear_solver_info_ID(FESystemNumbers::InvalidID),
max_iters(FESystemNumbers::InvalidID),
tolerance(0.0)
{
  
}


Solver::NonlinearSolverInfo::~NonlinearSolverInfo()
{
  
}


std::istream& 
Solver::NonlinearSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, NONLINEAR_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "SOLVER_TYPE", tag);
  
  this->nonlinear_solver_kind_enum_ID = Solver::NonlinearSolverKindEnum::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, "LINEAR_SOLVER_INFO_ID", 
                                   this->linear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "MAXIMUM_ITERATIONS", 
                                   this->max_iters);
  
  FESystemIOUtility::readFromInput(input, "CONVERGENCE_TOLERANCE", 
                                   this->tolerance);

  FESystemIOUtility::readFromInput(input, NONLINEAR_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}



std::istream& 
Solver::operator >> (std::istream& input, Solver::NonlinearSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}
