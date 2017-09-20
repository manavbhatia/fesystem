// $Id: EulerTransientSolverInfo.C,v 1.1.2.1 2008-02-25 04:31:21 manav Exp $

// FESystem includes
#include "Solvers/EulerTransientSolverInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/EulerTransientSolver.h"


Solver::EulerLinearTransientSolverInfo::EulerLinearTransientSolverInfo():
Solver::LinearTransientSolverInfo(Solver::EULER_LINEAR_TRANSIENT_SOLVER_INFO::num(),
                                  Solver::LINEAR_TRANSIENT_SOLVER::num()),
order(FESystemNumbers::InvalidID)
{
  
}


Solver::EulerLinearTransientSolverInfo::~EulerLinearTransientSolverInfo()
{
  
}


unsigned int
Solver::EulerLinearTransientSolverInfo::getSystemOrder() const
{
  Assert(this->order != FESystemNumbers::InvalidID, ExcInternalError());
  return this->order;
}




std::istream& 
Solver::EulerLinearTransientSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, Solver::EULER_LINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  this->transient_solver_kind_enum_ID = Solver::EULER_LINEAR_TRANSIENT_SOLVER::num();
  
  FESystemIOUtility::readFromInput(input, "LINEAR_SOLVER_INFO_ID", this->linear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->max_iterations);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_TIME", this->initial_time);
  
  FESystemIOUtility::readFromInput(input, "TOTAL_DURATION", this->total_duration);
  
  FESystemIOUtility::readFromInput(input, "IF_ADAPTIVE_SIZE", this->adaptive_step_size);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_STEP_SIZE", this->initial_step_size);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_TOLERANCE", this->global_tolerance);
  
  FESystemIOUtility::readFromInput(input, "SYSTEM_ORDER", this->order);
    
  FESystemIOUtility::readFromInput(input, Solver::EULER_LINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}



std::istream& 
Solver::operator >> (std::istream& input, Solver::EulerLinearTransientSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}



Solver::EulerNonlinearTransientSolverInfo::EulerNonlinearTransientSolverInfo():
Solver::NonlinearTransientSolverInfo(Solver::EULER_NONLINEAR_TRANSIENT_SOLVER_INFO::num(),
                                     Solver::NONLINEAR_TRANSIENT_SOLVER::num()),
order(FESystemNumbers::InvalidID)
{
  
}


Solver::EulerNonlinearTransientSolverInfo::~EulerNonlinearTransientSolverInfo()
{
  
}


unsigned int
Solver::EulerNonlinearTransientSolverInfo::getSystemOrder() const
{
  Assert(this->order != FESystemNumbers::InvalidID, ExcInternalError());
  return this->order;
}




std::istream& 
Solver::EulerNonlinearTransientSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, Solver::EULER_NONLINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  this->transient_solver_kind_enum_ID = Solver::EULER_NONLINEAR_TRANSIENT_SOLVER::num();
  
  FESystemIOUtility::readFromInput(input, "NONLINEAR_SOLVER_INFO_ID", this->nonlinear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->max_iterations);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_TIME", this->initial_time);
  
  FESystemIOUtility::readFromInput(input, "TOTAL_DURATION", this->total_duration);
  
  FESystemIOUtility::readFromInput(input, "IF_ADAPTIVE_SIZE", this->adaptive_step_size);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_STEP_SIZE", this->initial_step_size);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_TOLERANCE", this->global_tolerance);
  
  FESystemIOUtility::readFromInput(input, "SYSTEM_ORDER", this->order);
  
  FESystemIOUtility::readFromInput(input, Solver::EULER_NONLINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}




std::istream& 
Solver::operator >> (std::istream& input, Solver::EulerNonlinearTransientSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}
