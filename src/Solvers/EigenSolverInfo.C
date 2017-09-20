// $Id: EigenSolverInfo.C,v 1.7.4.2 2007-06-13 15:01:46 manav Exp $

// FESystem includes
#include "Solvers/EigenSolverInfo.h"
#include "Utilities/InputOutputUtility.h"




Solver::EigenSolverInfo::EigenSolverInfo():
Solver::SolverInfo(EIGEN_SOLVER_INFO::num(),
                   EIGEN_SOLVER::num()),
problem_kind_enum_ID(FESystemNumbers::InvalidID),
eigen_solver_kind_enum_ID(FESystemNumbers::InvalidID),
linear_solver_info_ID(0),
n_eigen_pairs(0),
calculate_eigen_vectors(false),
solver_shift_kind_enum_ID(FESystemNumbers::InvalidID),
solver_shift_value(0.0),
eigen_spectrum_end(FESystemNumbers::InvalidID),
tolerance_value(0.0),
maximum_iterations(FESystemNumbers::InvalidID),
swtich_A_and_B_matrices(false)
{
  
}


Solver::EigenSolverInfo::~EigenSolverInfo()
{
  
}


std::istream& 
Solver::EigenSolverInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, EIGEN_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "EIGEN_SOLVER_TYPE", tag);
  
  this->eigen_solver_kind_enum_ID = Solver::EigenSolverKindEnum::enumID(tag);

  tag.clear();
  FESystemIOUtility::readFromInput(input, "EIGEN_PROBLEM_TYPE", tag);
  
  this->problem_kind_enum_ID = Solver::EigenProblemKindEnum::enumID(tag);

  FESystemIOUtility::readFromInput(input, "SWAP_A_AND_B_MATRICES", 
                                   this->swtich_A_and_B_matrices);
  
  FESystemIOUtility::readFromInput(input, "LINEAR_SOLVER_INFO_ID", this->linear_solver_info_ID);
  
  FESystemIOUtility::readFromInput(input, "N_EIGEN_PAIRS", this->n_eigen_pairs);

  FESystemIOUtility::readFromInput(input, "CALCULATE_EIGEN_VECTORS", this->calculate_eigen_vectors);
  

  tag.clear();
  FESystemIOUtility::readFromInput(input, "EIGEN_SHIFT_TYPE", tag);
  this->solver_shift_kind_enum_ID = Solver::EigenSolutionShiftKindEnum::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, "EIGEN_SHIFT_VALUE", this->solver_shift_value);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "EIGEN_SPECTRUM_END_TO_CALCULATE", tag);
  
  this->eigen_spectrum_end = Solver::EigenSolutionSpectrumKindEnum::enumID(tag);

  FESystemIOUtility::readFromInput(input, "TOLERANCE", this->tolerance_value);

  FESystemIOUtility::readFromInput(input, "MAX_ITERATIONS", this->maximum_iterations);
    
  FESystemIOUtility::readFromInput(input, EIGEN_SOLVER_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  this->initialized = true;
  
  return input;
}



std::istream& 
Solver::operator >> (std::istream& input, Solver::EigenSolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}
