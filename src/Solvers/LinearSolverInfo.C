// $Id: LinearSolverInfo.C,v 1.3 2006-09-05 20:41:35 manav Exp $

// FESystem includes
#include "Solvers/LinearSolverInfo.h"
#include "Utilities/InputOutputUtility.h"




Solver::LinearSolverInfo::LinearSolverInfo():
Solver::SolverInfo(LINEAR_SOLVER_INFO::num(),
                   LINEAR_SOLVER::num()),
pc_type_enum_ID(FESystemNumbers::InvalidID),
ksp_type_enum_ID(FESystemNumbers::InvalidID)
{
  
}


Solver::LinearSolverInfo::~LinearSolverInfo()
{
  
}


std::istream& 
Solver::LinearSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, LINEAR_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);

  tag.clear();
  FESystemIOUtility::readFromInput(input, "SOLVER_PACKAGE", tag);
  this->linear_solver_enum_id = Solver::LinearSolverKindEnum::enumID(tag);

  tag.clear();
  FESystemIOUtility::readFromInput(input, "KSP_TYPE", tag);
  
  this->ksp_type_enum_ID = Solver::KSPTypeEnum::enumID(tag);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "PC_TYPE", tag);

  this->pc_type_enum_ID = Solver::PCTypeEnum::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, LINEAR_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}



std::istream& 
Solver::operator >> (std::istream& input, Solver::LinearSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}
