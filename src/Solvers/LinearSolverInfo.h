// $Id: LinearSolverInfo.h,v 1.3 2006-09-05 20:41:35 manav Exp $

#ifndef __fesystem_linear_solver_info_h__
#define __fesystem_linear_solver_info_h__


// FESystem includes
#include "Solvers/SolverInfo.h"
#include "Solvers/LinearSolver.h"


#ifndef LINEAR_SOLVER_INFO_ENUM_ID
#define LINEAR_SOLVER_INFO_ENUM_ID 1
#else
#error
#endif

#ifndef LINEAR_SOLVER_INFO_ENUM_NAME
#define LINEAR_SOLVER_INFO_ENUM_NAME "LINEAR_SOLVER_INFO"
#else
#error
#endif


namespace Solver
{
  
  
  DeclareEnumName(LINEAR_SOLVER_INFO, Solver::SolverInfoEnum,
                  LINEAR_SOLVER_INFO_ENUM_ID,
                  LINEAR_SOLVER_INFO_ENUM_NAME);
  
  
  
  class LinearSolverInfo: public SolverInfo
    {
public:
      
      LinearSolverInfo();
      virtual ~LinearSolverInfo();
      
      inline unsigned int getLinearSolverEnumID() const;
      
      inline unsigned int getPCTypeEnumID() const;
      inline std::string getPCTypeEnumName() const;

      inline unsigned int getKSPTypeEnumID() const;
      inline std::string getKSPTypeEnumName() const;

      virtual std::istream& readFromInputStream(std::istream& input);
      
      friend std::istream& operator >> (std::istream& input, 
                                        Solver::LinearSolverInfo& info);
      
protected:
        
      unsigned int linear_solver_enum_id;
      unsigned int pc_type_enum_ID;
      unsigned int ksp_type_enum_ID;
    };
  
}


inline 
unsigned int 
Solver::LinearSolverInfo::getLinearSolverEnumID() const
{
  Assert(this->initialized,ExcInvalidState());
  return this->linear_solver_enum_id;
}


inline 
unsigned int 
Solver::LinearSolverInfo::getPCTypeEnumID() const
{
  Assert(this->initialized,ExcInvalidState());
  return this->pc_type_enum_ID;
}



inline 
std::string
Solver::LinearSolverInfo::getPCTypeEnumName() const
{
  Assert(this->initialized,ExcInvalidState());
  return Solver::PCTypeEnum::enumName(this->pc_type_enum_ID);
}


inline 
unsigned int 
Solver::LinearSolverInfo::getKSPTypeEnumID() const
{
  Assert(this->initialized,ExcInvalidState());
  return this->ksp_type_enum_ID;
}



inline 
std::string
Solver::LinearSolverInfo::getKSPTypeEnumName() const
{
  Assert(this->initialized,ExcInvalidState());
  return Solver::KSPTypeEnum::enumName(this->ksp_type_enum_ID);
}


#endif // __fesystem_linear_solver_info_h__
