// $Id: TransientSolverInfo.h,v 1.1.2.3 2007-05-11 05:16:54 manav Exp $

#ifndef __fesystem_transient_solver_info_h__
#define __fesystem_transient_solver_info_h__


// FESystem includes
#include "Solvers/SolverInfo.h"
#include "Solvers/TransientSolver.h"


#ifndef LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID
#define LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID 4
#else
#error
#endif

#ifndef LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME
#define LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME "LINEAR_TRANSIENT_SOLVER_INFO"
#else
#error
#endif


#ifndef NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID
#define NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID 5
#else
#error
#endif

#ifndef NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME
#define NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME "NONLINEAR_TRANSIENT_SOLVER_INFO"
#else
#error
#endif



namespace Solver
{
  
  
  DeclareEnumName(LINEAR_TRANSIENT_SOLVER_INFO, Solver::SolverInfoEnum,
                  LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID,
                  LINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME);
  
  DeclareEnumName(NONLINEAR_TRANSIENT_SOLVER_INFO, Solver::SolverInfoEnum,
                  NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID,
                  NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_NAME);
  
  
  /// stores the information about an eigen solver. This information is used to configure
  /// the eigen solver
  class TransientSolverInfoBase: public SolverInfo
    {
public:
      /// constructor
      TransientSolverInfoBase(const unsigned int info_enum_ID, 
                              const unsigned int solver_class);
      
      /// destructor
      virtual ~TransientSolverInfoBase();
      
      /// @returns eigen solver kind enum ID
      inline unsigned int getTransientSolverKindEnumID() const;
      
      /// @returns eigen solver kind enum name
      virtual std::string getTransientSolverKindEnumName() const;
      
      /// @returns total duration of time integration
      inline double getInitialTime() const;
      
      /// @returns total duration of time integration
      inline double getDuration() const;
      
      /// @returns the final time of the transient simulation
      double getFinalTime() const;

      /// @returns the maximum number of time iterations allowed
      inline unsigned int getMaxTimeIterations() const;
      
      /// @returns true/false based on whether the user has asked for 
      /// adaptive calculation of time steps. If this is set to true, 
      /// then the initial time step would be ignored, and the global 
      /// tolerance level will be specified.
      inline bool ifAdaptiveTimeSteps() const;
      
      /// @returns the value of the initial time step. This is used 
      /// when adaptive time steps is not used.
      inline double getInitialTimeStep() const;
      
      /// @returns global tolerance level. This is used only when adaptive 
      /// time steps is turned on.
      inline double getGlobalTolerance() const;
      
      /// reads from input
      virtual std::istream& readFromInputStream(std::istream& input) = 0;
      
      
protected:
        
        /// kind of eigen solver
        unsigned int transient_solver_kind_enum_ID;
      
      /// number of eigen pairs to compute
      unsigned int max_iterations;
      
      /// initial time to start the iteration at
      double initial_time;
      
      /// total duration of simulation
      double total_duration;
      
      /// if adaptive step size should be used or not
      bool adaptive_step_size;
      
      /// value of initial step size
      double initial_step_size;
      
      /// global tolerance used for adaptive time calculation
      double global_tolerance;
    };
  
  
  /// this is a class specifically for linear transient solvers
  class LinearTransientSolverInfo: public Solver::TransientSolverInfoBase
    {
public: 
      /// constructor
      LinearTransientSolverInfo();
      
      /// constructor for cases where another class derives from this class
      LinearTransientSolverInfo(const unsigned int info_enum_ID, 
                                const unsigned int solver_class_enum_ID);
      
      /// destructor
      ~LinearTransientSolverInfo();
      
      /// reads from input
      virtual std::istream& readFromInputStream(std::istream& input);
      
      /// @returns linear solver info ID for this solver
      unsigned int getLinearSolverInfoID() const;
            
      /// overloaded input stream operator
      friend std::istream& operator >> (std::istream& input, 
                                        Solver::LinearTransientSolverInfo& info);
      
protected:
        
        /// ID of linear solver info
        unsigned int linear_solver_info_ID;      

    };
  
  /// this is a class specifically for linear transient solvers
  class NonlinearTransientSolverInfo: public Solver::TransientSolverInfoBase
    {
public: 
      /// constructor
      NonlinearTransientSolverInfo();
      
      /// constructor for cases where another class derives from this class
      NonlinearTransientSolverInfo(const unsigned int info_enum_ID, 
                                   const unsigned int solver_class_enum_ID);
      
      /// destructor
      ~NonlinearTransientSolverInfo();
      
      /// @returns nonlinear solver info ID for this solver
      unsigned int getNonlinearSolverInfoID() const;
      
      /// reads from input
      virtual std::istream& readFromInputStream(std::istream& input);
      
      /// overloaded input stream operator
      friend std::istream& operator >> (std::istream& input, 
                                        Solver::NonlinearTransientSolverInfo& info);
      
protected:
        
        /// ID of nonlinear solver info
        unsigned int nonlinear_solver_info_ID;
      
    };
  
}



inline 
unsigned int 
Solver::TransientSolverInfoBase::getTransientSolverKindEnumID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->transient_solver_kind_enum_ID;
}




inline 
std::string
Solver::TransientSolverInfoBase::getTransientSolverKindEnumName() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return Solver::TransientSolverKindEnum::enumName(this->transient_solver_kind_enum_ID);
}





inline 
double
Solver::TransientSolverInfoBase::getInitialTime() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->initial_time;
}


inline 
double
Solver::TransientSolverInfoBase::getDuration() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->total_duration;
}


inline 
double
Solver::TransientSolverInfoBase::getFinalTime() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return (this->total_duration + this->initial_time);
}



inline 
unsigned int 
Solver::TransientSolverInfoBase::getMaxTimeIterations() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->max_iterations;
}



inline
bool
Solver::TransientSolverInfoBase::ifAdaptiveTimeSteps() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->adaptive_step_size;
}



inline
double
Solver::TransientSolverInfoBase::getInitialTimeStep() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->initial_step_size;
}


inline
double
Solver::TransientSolverInfoBase::getGlobalTolerance() const
{
  Assert(this->initialized,
         ExcInvalidState());
  
  // only if adaptive step size is enabled
  AssertThrow(this->adaptive_step_size, ExcInternalError());
  
  return this->global_tolerance;
}



inline 
unsigned int 
Solver::LinearTransientSolverInfo::getLinearSolverInfoID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->linear_solver_info_ID;
}




inline 
unsigned int 
Solver::NonlinearTransientSolverInfo::getNonlinearSolverInfoID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->nonlinear_solver_info_ID;
}


#endif // __fesystem_transient_solver_info_h__
