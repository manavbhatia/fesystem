// $Id: NewmarkTransientSolverInfo.h,v 1.1.2.2 2008-02-25 04:32:21 manav Exp $

#ifndef __fesystem_newmark_transient_solver_info_h__
#define __fesystem_newmark_transient_solver_info_h__


// FESystem includes
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/SolverInfo.h"
#include "Solvers/TransientSolver.h"


#ifndef NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID
#define NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID 6
#else
#error
#endif

#ifndef NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME
#define NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME "NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO"
#else
#error
#endif


#ifndef NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID
#define NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID 7
#else
#error
#endif

#ifndef NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME
#define NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME "NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO"
#else
#error
#endif



namespace Solver
{
  
  
  DeclareEnumName(NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO, Solver::SolverInfoEnum,
                  NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID,
                  NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME);
  
  DeclareEnumName(NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO, Solver::SolverInfoEnum,
                  NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID,
                  NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME);

  

/// this is a class specifically for linear transient solvers
class NewmarkLinearTransientSolverInfo: public Solver::LinearTransientSolverInfo
{
public: 
  /// constructor
  NewmarkLinearTransientSolverInfo();
  
  /// destructor
  ~NewmarkLinearTransientSolverInfo();

  /// @returns the order of time derivative of this system
  unsigned int getSystemOrder() const;
  
  /// @returns the constant for the given order of time derivative
  double getConstant(const unsigned int order) const;
  
  /// reads from input
  virtual std::istream& readFromInputStream(std::istream& input);
  
  /// overloaded input stream operator
  friend std::istream& operator >> (std::istream& input, 
                                    Solver::NewmarkLinearTransientSolverInfo& info);
  
protected:

    /// order of the system
    unsigned int order;
    
  /// vector of constants
  std::vector<double> constant_vals;
};


/// this is a class specifically for linear transient solvers
class NewmarkNonlinearTransientSolverInfo: public Solver::NonlinearTransientSolverInfo
{
public: 
  /// constructor
  NewmarkNonlinearTransientSolverInfo();
  
  /// destructor
  ~NewmarkNonlinearTransientSolverInfo();
  
  /// @returns the order of time derivative of this system
  unsigned int getSystemOrder() const;
  
  /// @returns the constant for the given order of time derivative
  double getConstant(const unsigned int order) const;
  
  /// reads from input
  virtual std::istream& readFromInputStream(std::istream& input);
  
  /// overloaded input stream operator
  friend std::istream& operator >> (std::istream& input, 
                                    Solver::NewmarkNonlinearTransientSolverInfo& info);

protected:
    
    /// order of the system
    unsigned int order;
  
  /// vector of constants
  std::vector<double> constant_vals;
  
};

}



#endif // __fesystem_newmark_transient_solver_info_h__
