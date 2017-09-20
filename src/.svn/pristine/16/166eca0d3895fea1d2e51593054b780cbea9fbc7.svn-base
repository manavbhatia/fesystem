// $Id: EulerTransientSolverInfo.h,v 1.1.2.1 2008-02-25 04:31:21 manav Exp $

#ifndef __fesystem_euler_transient_solver_info_h__
#define __fesystem_euler_transient_solver_info_h__


// FESystem includes
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/SolverInfo.h"
#include "Solvers/TransientSolver.h"


#ifndef EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID
#define EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID 8
#else
#error
#endif

#ifndef EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME
#define EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME "EULER_LINEAR_TRANSIENT_SOLVER_INFO"
#else
#error
#endif


#ifndef EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID
#define EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID 9
#else
#error
#endif

#ifndef EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME
#define EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME "EULER_NONLINEAR_TRANSIENT_SOLVER_INFO"
#else
#error
#endif



namespace Solver
{
  
  
  DeclareEnumName(EULER_LINEAR_TRANSIENT_SOLVER_INFO, Solver::SolverInfoEnum,
                  EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID,
                  EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME);
  
  DeclareEnumName(EULER_NONLINEAR_TRANSIENT_SOLVER_INFO, Solver::SolverInfoEnum,
                  EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID,
                  EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME);
  
  
  
  /// this is a class specifically for linear transient solvers
  class EulerLinearTransientSolverInfo: public Solver::LinearTransientSolverInfo
  {
  public: 
    /// constructor
    EulerLinearTransientSolverInfo();
    
    /// destructor
    ~EulerLinearTransientSolverInfo();
    
    /// @returns the order of time derivative of this system
    unsigned int getSystemOrder() const;
        
    /// reads from input
    virtual std::istream& readFromInputStream(std::istream& input);
    
    /// overloaded input stream operator
    friend std::istream& operator >> (std::istream& input, 
                                      Solver::EulerLinearTransientSolverInfo& info);
    
  protected:
    
    /// order of the system
    unsigned int order;
  };
  
  
  /// this is a class specifically for linear transient solvers
  class EulerNonlinearTransientSolverInfo: public Solver::NonlinearTransientSolverInfo
  {
  public: 
    /// constructor
    EulerNonlinearTransientSolverInfo();
    
    /// destructor
    ~EulerNonlinearTransientSolverInfo();
    
    /// @returns the order of time derivative of this system
    unsigned int getSystemOrder() const;
        
    /// reads from input
    virtual std::istream& readFromInputStream(std::istream& input);
    
    /// overloaded input stream operator
    friend std::istream& operator >> (std::istream& input, 
                                      Solver::EulerNonlinearTransientSolverInfo& info);
    
  protected:
    
    /// order of the system
    unsigned int order;
  };
  
}



#endif // __fesystem_euler_transient_solver_info_h__
