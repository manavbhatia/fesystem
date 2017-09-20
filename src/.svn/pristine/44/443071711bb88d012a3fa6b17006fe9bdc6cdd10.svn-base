// $Id: NonlinearSolverInfo.h,v 1.5 2006-09-05 20:41:35 manav Exp $

#ifndef __fesystem_nonlinear_solver_info_h__
#define __fesystem_nonlinear_solver_info_h__


// FESystem includes
#include "Solvers/SolverInfo.h"
#include "Solvers/NonlinearSolver.h"


#ifndef NONLINEAR_SOLVER_INFO_ENUM_ID
#define NONLINEAR_SOLVER_INFO_ENUM_ID 2
#else
#error
#endif

#ifndef NONLINEAR_SOLVER_INFO_ENUM_NAME
#define NONLINEAR_SOLVER_INFO_ENUM_NAME "NONLINEAR_SOLVER_INFO"
#else
#error
#endif

namespace Solver
{
  
  DeclareEnumName(NONLINEAR_SOLVER_INFO, SolverInfoEnum,
                  NONLINEAR_SOLVER_INFO_ENUM_ID,
                  NONLINEAR_SOLVER_INFO_ENUM_NAME);
  
  /// this class stores the information for initialization of a nonlinear solver.
  class NonlinearSolverInfo: public SolverInfo
    {
public:
      
      /// constructor
      NonlinearSolverInfo();
      
      /// destructor
      virtual ~NonlinearSolverInfo();
      
      /// @returns enum ID of the nonlinear solver kind
      virtual inline unsigned int getNonlinearSolverKindEnumID() const;

      /// @returns enum name of the nonlinear solver kind
      virtual inline const std::string getNonlinearSolverKindEnumName() const;
      
      /// @returns the solver info ID for LinearSolver to be used for this
      /// NonlinearSolver
      inline unsigned int getLinearSolverInfoID() const;

      /// @returns the maximum number of iterations for this solver
      inline unsigned int getMaximumIterations() const;

      /// @returns the maximum tolerance to define convergence of this solver
      inline double getConvergenceTolerance() const;

      /// reads from the input stream
      virtual std::istream& readFromInputStream(std::istream& input);
      
      /// overloaded stream input operator for this class
      friend std::istream& operator >> (std::istream& input, 
                                        Solver::NonlinearSolverInfo& info);
      
protected:

        /// solver kind enum ID
      unsigned int nonlinear_solver_kind_enum_ID;
      
      /// LinearSolverInfo ID
      unsigned int linear_solver_info_ID;
      
      /// max iterations
      unsigned int max_iters;
      
      /// tolerance
      double tolerance;
    };
}

inline 
unsigned int 
Solver::NonlinearSolverInfo::getNonlinearSolverKindEnumID() const
{
  return this->nonlinear_solver_kind_enum_ID;
}



inline 
const std::string
Solver::NonlinearSolverInfo::getNonlinearSolverKindEnumName() const
{
  return Solver::NonlinearSolverKindEnum::enumName(this->nonlinear_solver_kind_enum_ID);
}


inline 
unsigned int 
Solver::NonlinearSolverInfo::getLinearSolverInfoID() const
{
  return this->linear_solver_info_ID;
}


inline
unsigned int
Solver::NonlinearSolverInfo::getMaximumIterations() const
{
  return this->max_iters;
}


inline
double 
Solver::NonlinearSolverInfo::getConvergenceTolerance() const
{
  return this->tolerance;
}



#endif // __fesystem_nonlinear_solver_info_h__
