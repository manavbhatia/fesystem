// $Id: ThermalSolution.h,v 1.1.6.1 2007-03-14 22:05:03 manav Exp $

#ifndef __fesystem_thermal_solution_h__
#define __fesystem_thermal_solution_h__


// C++ includes
#include <iostream>
#include <string>

// FESystem includes
#include "Solutions/SolutionBase.h"

// forward declerations
namespace FESystem
{
  class FESystemController;
}

#ifndef THERMAL_SOLUTION_ENUM_ID 
#define THERMAL_SOLUTION_ENUM_ID 1
#else 
#error
#endif

#ifndef THERMAL_SOLUTION_ENUM_NAME 
#define THERMAL_SOLUTION_ENUM_NAME "THERMAL_SOLUTION"
#else 
#error
#endif



namespace Solution
{
  
  DeclareEnumName(THERMAL_SOLUTION, Solution::SolutionEnum,
                  THERMAL_SOLUTION_ENUM_ID, THERMAL_SOLUTION_ENUM_NAME);
  
  
  /// this class provides the base functionality for all solutions
  class ThermalSolution: public Solution::SolutionBase
    {
public:
      /// constructor
      ThermalSolution();
      
      /// destructor
      virtual ~ThermalSolution();
      
      
      /// solves the problem
      virtual void solve();
      
    /// @returns the transient nature of the given discipline
    virtual unsigned int getDisciplineTransientNatureEnumID(const unsigned int discipline) const ;

    /// @returns the vector of SolNameInfo structures for the nodal solutions, and for the
    /// specified disciplines
    virtual FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
      getTimeIndependentNodalSolutionDataInfo(const unsigned int discipline_enum_ID) ;
    
    /// @returns the vector of SolNameInfo structures for the nodal solutions, and for the
    /// specified disciplines
    virtual FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
      getTimeDependentNodalSolutionDataInfo(const unsigned int discipline_enum_ID) ;

      /// @this will read from the input stream
      virtual std::istream& readFromInputStream(std::istream& input);

protected:
        
    };
}


#endif // __fesystem_thermal_solution_h__
