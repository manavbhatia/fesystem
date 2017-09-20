// $Id: TransientSolverInfo.C,v 1.1.2.4 2007-05-11 05:16:54 manav Exp $

// FESystem includes
#include "Solvers/TransientSolverInfo.h"
#include "Utilities/InputOutputUtility.h"




Solver::TransientSolverInfoBase::TransientSolverInfoBase(const unsigned int info_enum_ID, 
                                                         const unsigned int solver_class):
Solver::SolverInfo(info_enum_ID, solver_class),
transient_solver_kind_enum_ID(FESystemNumbers::InvalidID),
max_iterations(FESystemNumbers::InvalidID),
initial_time(0.0),
total_duration(0.0),
adaptive_step_size(false),
initial_step_size(0.0),
global_tolerance(0.0)
{
  
}



Solver::TransientSolverInfoBase::~TransientSolverInfoBase()
{
  
}


Solver::LinearTransientSolverInfo::LinearTransientSolverInfo():
Solver::TransientSolverInfoBase(Solver::LINEAR_TRANSIENT_SOLVER_INFO::num(),
				Solver::LINEAR_TRANSIENT_SOLVER::num()),
linear_solver_info_ID(FESystemNumbers::InvalidID)
{
  
}



Solver::LinearTransientSolverInfo::LinearTransientSolverInfo(const unsigned int info_enum_ID,
                                                             const unsigned int solver_class):
Solver::TransientSolverInfoBase(info_enum_ID,
                                solver_class),
linear_solver_info_ID(FESystemNumbers::InvalidID)
{
  
}



Solver::LinearTransientSolverInfo::~LinearTransientSolverInfo()
{
  
}


std::istream& 
Solver::LinearTransientSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, Solver::LINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "TRANSIENT_SOLVER_TYPE", tag);
  
  this->transient_solver_kind_enum_ID = Solver::TransientSolverKindEnum::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, "LINEAR_SOLVER_INFO_ID", this->linear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->max_iterations);

  FESystemIOUtility::readFromInput(input, "INITIAL_TIME", this->initial_time);

  FESystemIOUtility::readFromInput(input, "TOTAL_DURATION", this->total_duration);
  
  FESystemIOUtility::readFromInput(input, "IF_ADAPTIVE_SIZE", this->adaptive_step_size);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_STEP_SIZE", this->initial_step_size);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_TOLERANCE", this->global_tolerance);
  
  FESystemIOUtility::readFromInput(input, Solver::LINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}



std::istream& 
Solver::operator >> (std::istream& input, Solver::LinearTransientSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}



Solver::NonlinearTransientSolverInfo::NonlinearTransientSolverInfo():
Solver::TransientSolverInfoBase(Solver::NONLINEAR_TRANSIENT_SOLVER_INFO::num(),
                                Solver::NONLINEAR_TRANSIENT_SOLVER::num()),
nonlinear_solver_info_ID(FESystemNumbers::InvalidID)
{
  
}



Solver::NonlinearTransientSolverInfo::NonlinearTransientSolverInfo(const unsigned int info_enum_ID,
                                                                   const unsigned int solver_class):
Solver::TransientSolverInfoBase(info_enum_ID,
                                solver_class)
{
  
}



Solver::NonlinearTransientSolverInfo::~NonlinearTransientSolverInfo()
{
  
}



std::istream& 
Solver::NonlinearTransientSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, Solver::NONLINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "TRANSIENT_SOLVER_TYPE", tag);
  
  this->transient_solver_kind_enum_ID = Solver::TransientSolverKindEnum::enumID(tag);

  FESystemIOUtility::readFromInput(input, "NONLINEAR_SOLVER_INFO_ID", this->nonlinear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->max_iterations);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_TIME", this->initial_time);
  
  FESystemIOUtility::readFromInput(input, "TOTAL_DURATION", this->total_duration);
  
  FESystemIOUtility::readFromInput(input, "IF_ADAPTIVE_SIZE", this->adaptive_step_size);
  
  FESystemIOUtility::readFromInput(input, "INITIAL_STEP_SIZE", this->initial_step_size);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_TOLERANCE", this->global_tolerance);
  
  FESystemIOUtility::readFromInput(input, Solver::NONLINEAR_TRANSIENT_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}




std::istream& 
Solver::operator >> (std::istream& input, Solver::NonlinearTransientSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}


