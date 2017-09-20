// $Id: NewmarkTransientSolverInfo.C,v 1.1.2.1 2007-05-15 20:38:53 manav Exp $

// FESystem includes
#include "Solvers/NewmarkTransientSolverInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/NewmarkTransientSolver.h"


Solver::NewmarkLinearTransientSolverInfo::NewmarkLinearTransientSolverInfo():
Solver::LinearTransientSolverInfo(Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO::num(),
                                  Solver::LINEAR_TRANSIENT_SOLVER::num()),
order(FESystemNumbers::InvalidID)
{
  
}


Solver::NewmarkLinearTransientSolverInfo::~NewmarkLinearTransientSolverInfo()
{
  
}


unsigned int
Solver::NewmarkLinearTransientSolverInfo::getSystemOrder() const
{
  Assert(this->order != FESystemNumbers::InvalidID, ExcInternalError());
  return this->order;
}



double
Solver::NewmarkLinearTransientSolverInfo::getConstant(const unsigned int i) const
{
  Assert(i < this->getSystemOrder() != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->constant_vals.size() == this->getSystemOrder(), ExcInternalError());
  return this->constant_vals[i];
}



std::istream& 
Solver::NewmarkLinearTransientSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  this->transient_solver_kind_enum_ID = Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER::num();
  
  FESystemIOUtility::readFromInput(input, "LINEAR_SOLVER_INFO_ID", this->linear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->max_iterations);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_TIME", this->initial_time);
  
  FESystemIOUtility::readFromInput(input, "TOTAL_DURATION", this->total_duration);
  
  FESystemIOUtility::readFromInput(input, "IF_ADAPTIVE_SIZE", this->adaptive_step_size);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_STEP_SIZE", this->initial_step_size);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_TOLERANCE", this->global_tolerance);
  
  FESystemIOUtility::readFromInput(input, "SYSTEM_ORDER", this->order);

  FESystemIOUtility::readFromInput(input, "CONSTANTS");
  
  this->constant_vals.resize(order);
  for (unsigned int i=0; i < this->order; i++)
    input >> this->constant_vals[i];
  
  FESystemIOUtility::readFromInput(input, Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}



std::istream& 
Solver::operator >> (std::istream& input, Solver::NewmarkLinearTransientSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}



Solver::NewmarkNonlinearTransientSolverInfo::NewmarkNonlinearTransientSolverInfo():
Solver::NonlinearTransientSolverInfo(Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO::num(),
                                     Solver::NONLINEAR_TRANSIENT_SOLVER::num()),
order(FESystemNumbers::InvalidID)
{
  
}


Solver::NewmarkNonlinearTransientSolverInfo::~NewmarkNonlinearTransientSolverInfo()
{
  
}


unsigned int
Solver::NewmarkNonlinearTransientSolverInfo::getSystemOrder() const
{
  Assert(this->order != FESystemNumbers::InvalidID, ExcInternalError());
  return this->order;
}



double
Solver::NewmarkNonlinearTransientSolverInfo::getConstant(const unsigned int i) const
{
  Assert(i < this->getSystemOrder() != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->constant_vals.size() == this->getSystemOrder(), ExcInternalError());
  return this->constant_vals[i];
}


std::istream& 
Solver::NewmarkNonlinearTransientSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  this->transient_solver_kind_enum_ID = Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER::num();
  
  FESystemIOUtility::readFromInput(input, "NONLINEAR_SOLVER_INFO_ID", this->nonlinear_solver_info_ID);
    
  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->max_iterations);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_TIME", this->initial_time);
  
  FESystemIOUtility::readFromInput(input, "TOTAL_DURATION", this->total_duration);
  
  FESystemIOUtility::readFromInput(input, "IF_ADAPTIVE_SIZE", this->adaptive_step_size);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_STEP_SIZE", this->initial_step_size);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_TOLERANCE", this->global_tolerance);
  
  FESystemIOUtility::readFromInput(input, "SYSTEM_ORDER", this->order);
  
  FESystemIOUtility::readFromInput(input, "CONSTANTS");
  
  this->constant_vals.resize(order);
  for (unsigned int i=0; i < this->order; i++)
    input >> this->constant_vals[i];
  
  FESystemIOUtility::readFromInput(input, Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}




std::istream& 
Solver::operator >> (std::istream& input, Solver::NewmarkNonlinearTransientSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}
