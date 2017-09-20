// $Id: SolverInfo.h,v 1.4.4.1 2007-03-04 03:37:13 manav Exp $

#ifndef __fesystem_solver_info_h__
#define __fesystem_solver_info_h__

// C++ includes
#include <iostream>
#include <string>
#include <memory>


// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"
#include "Utilities/NameEnumHandler.h"
#include "Solvers/FESystemSolverBase.h"

namespace Solver
{
  DeclareEnumClass(SolverInfoEnum);
  
  class SolverInfo
    {
public:
      SolverInfo(const unsigned int info_enum_ID,
                 const unsigned int solver_class);
      
      virtual  ~SolverInfo();
      
      inline unsigned int getID() const;
      virtual inline unsigned int getSolverInfoEnumID() const;
      virtual inline const std::string getSolverInfoEnumName() const;
      
      inline unsigned int getSolverClassEnumID() const;
      inline std::string getSolverClassEnumName() const;
      
      virtual std::istream& readFromInputStream(std::istream& input) = 0;
  
      /// overloaded input operator
      friend std::istream& operator>>(std::istream& input, Solver::SolverInfo& info);
      
protected:
        
      unsigned int ID;
      
      const unsigned int solver_info_enum_ID;
      
      unsigned int solver_class_enum_ID;
      
      bool initialized;
    };
  
  /// @returns a pointer to a newly created SolverInfo
  std::auto_ptr<Solver::SolverInfo> createSolverInfo(const unsigned int enum_ID);
}




inline
unsigned int
Solver::SolverInfo::getID() const
{
  Assert(this->ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(this->ID));
  return this->ID;
}



inline
unsigned int
Solver::SolverInfo::getSolverInfoEnumID() const
{
  return this->solver_info_enum_ID;
}



inline
const std::string
Solver::SolverInfo::getSolverInfoEnumName() const
{
  return SolverInfoEnum::enumName(this->solver_info_enum_ID);
}



inline
unsigned int 
Solver::SolverInfo::getSolverClassEnumID() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return this->solver_class_enum_ID;
}



inline
std::string
Solver::SolverInfo::getSolverClassEnumName() const
{
  Assert(this->initialized,
         ExcInvalidState());
  return Solver::SolverClassEnum::enumName(this->solver_class_enum_ID);
}


#endif // __fesystem_solver_info_h__
