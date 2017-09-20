// $Id: EigenSolverInfo.h,v 1.8.4.2 2007-06-13 15:01:46 manav Exp $

#ifndef __fesystem_eigen_solver_info_h__
#define __fesystem_eigen_solver_info_h__


// FESystem includes
#include "Solvers/SolverInfo.h"
#include "Solvers/EigenSolver.h"
#include "Solvers/SlepcEigenSolver.h"


#ifndef EIGEN_SOLVER_INFO_ENUM_ID
#define EIGEN_SOLVER_INFO_ENUM_ID 3
#else
#error
#endif

#ifndef EIGEN_SOLVER_INFO_ENUM_NAME
#define EIGEN_SOLVER_INFO_ENUM_NAME "EIGEN_SOLVER_INFO"
#else
#error
#endif


namespace Solver
{
  
  
  DeclareEnumName(EIGEN_SOLVER_INFO, Solver::SolverInfoEnum,
                  EIGEN_SOLVER_INFO_ENUM_ID,
                  EIGEN_SOLVER_INFO_ENUM_NAME);
  
  
  /// stores the information about an eigen solver. This information is used to configure
  /// the eigen solver
  class EigenSolverInfo: public SolverInfo
    {
public:
      /// constructor
      EigenSolverInfo();
      
      /// destructor
      virtual ~EigenSolverInfo();
      
      /// @returns problem kind enum ID
      inline unsigned int getEigenProblemKindEnumID() const;

      /// @returns problem kind enum name
      inline std::string getEigenProblemKindEnumName() const;
      
      /// @returns eigen solver kind enum ID
      inline unsigned int getEigenSolverKindEnumID() const;

      /// @returns eigen solver kind enum name
      inline std::string getEigenSolverKindEnumName() const;
      
      /// @returns linear solver info ID for this solver
      unsigned int getLinearSolverInfoID() const;
      
      /// @returns number of eigen pairs to solve
      inline unsigned int getNEigenPairsToSolve() const;
      
      /// @returns the type of shift to apply to the eigen solver
      inline unsigned int getEigenSolverShiftKind() const;
      
      /// @returns the value of the eigen solver shift
      inline double getEigenSolverShiftValue() const;
      
      /// @returns true/false if the eigen vectors should be calculated
      inline bool calculateEigenVectors() const;
      
      /// @returns the eigen spectrum end to compute
      inline unsigned int getEigenSpectrumToCompute() const;
      
      /// @returns the tolerance of the solver
      inline double getTolerance() const;

      /// @returns the tolerance of the solver
      inline unsigned int getMaxAllowableIterations() const;

      //// @returns boolean stating whether the A and B matrices should be switched
      /// before using the other parameters that the user has set for the
      /// solver
      bool ifSwapMatrices() const;
      
      /// reads from input
      virtual std::istream& readFromInputStream(std::istream& input);
      
      /// overloaded input stream operator
      friend std::istream& operator >> (std::istream& input, 
                                        Solver::EigenSolverInfo& info);
      
protected:
      
      /// kind of eigen problem
      unsigned int problem_kind_enum_ID;

      /// kind of eigen solver
      unsigned int eigen_solver_kind_enum_ID;

      /// ID of linear solver info
      unsigned int linear_solver_info_ID;
      
      /// number of eigen pairs to compute
      unsigned int n_eigen_pairs;

      /// if eigen vectors should also be computed
      bool calculate_eigen_vectors;
      
      /// type of shift to apply to the solver
      unsigned int solver_shift_kind_enum_ID;
      
      /// value of shift
      double solver_shift_value;
      
      /// spectrum of eigen values to compute
      unsigned int eigen_spectrum_end;
      
      /// value of tolerance
      double tolerance_value;

      /// value of tolerance
      unsigned int maximum_iterations;
      
      /// whether to switch the A and B matrices or not
      bool swtich_A_and_B_matrices;
};
  
}


inline 
unsigned int 
Solver::EigenSolverInfo::getEigenProblemKindEnumID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->problem_kind_enum_ID;
}



inline 
std::string
Solver::EigenSolverInfo::getEigenProblemKindEnumName() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return Solver::EigenProblemKindEnum::enumName(this->problem_kind_enum_ID);
}


inline 
unsigned int 
Solver::EigenSolverInfo::getEigenSolverKindEnumID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->eigen_solver_kind_enum_ID;
}



inline 
std::string
Solver::EigenSolverInfo::getEigenSolverKindEnumName() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return Solver::EigenSolverKindEnum::enumName(this->eigen_solver_kind_enum_ID);
}


inline 
unsigned int 
Solver::EigenSolverInfo::getLinearSolverInfoID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->linear_solver_info_ID;
}



inline 
unsigned int 
Solver::EigenSolverInfo::getNEigenPairsToSolve() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->n_eigen_pairs;
}


inline 
bool
Solver::EigenSolverInfo::calculateEigenVectors() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->calculate_eigen_vectors;
}



inline 
unsigned int 
Solver::EigenSolverInfo::getEigenSpectrumToCompute() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->eigen_spectrum_end;
}



inline
unsigned int
Solver::EigenSolverInfo::getEigenSolverShiftKind() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->solver_shift_kind_enum_ID;
}



inline
double
Solver::EigenSolverInfo::getEigenSolverShiftValue() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->solver_shift_value;
}


inline
double
Solver::EigenSolverInfo::getTolerance() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->tolerance_value;
}


inline
unsigned int 
Solver::EigenSolverInfo::getMaxAllowableIterations() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->maximum_iterations;
}


inline
bool
Solver::EigenSolverInfo::ifSwapMatrices() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->swtich_A_and_B_matrices;
}


#endif // __fesystem_eigen_solver_info_h__
