// $Id: LinearStressSolution.h,v 1.4.6.3 2007-06-13 14:59:40 manav Exp $

#ifndef __fesystem_linear_stress_solution_h__
#define __fesystem_linear_stress_solution_h__


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

#ifndef LINEAR_STRESS_SOLUTION_ENUM_ID 
#define LINEAR_STRESS_SOLUTION_ENUM_ID 2
#else 
#error
#endif

#ifndef LINEAR_STRESS_SOLUTION_ENUM_NAME 
#define LINEAR_STRESS_SOLUTION_ENUM_NAME "LINEAR_STRESS_SOLUTION"
#else 
#error
#endif



namespace Solution
{
  
  DeclareEnumName(LINEAR_STRESS_SOLUTION, Solution::SolutionEnum,
                  LINEAR_STRESS_SOLUTION_ENUM_ID, LINEAR_STRESS_SOLUTION_ENUM_NAME);
  
  
  /// this class provides the base functionality for all solutions
  class LinearStressSolution: public Solution::SolutionBase
{
public:
      /// constructor
      LinearStressSolution();
      
      /// destructor
      virtual ~LinearStressSolution();

      
      /// solves the system
      virtual void solve();
      
    /// @returns the transient nature of the given discipline
    virtual unsigned int getDisciplineTransientNatureEnumID
      (const unsigned int discipline) const;
    
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
        
        /// boolean to check if coupled thermoelastic solution is needed, which will 
        /// need solution of the thermal discipline
        bool coupled_thermoelastic;
      
      /// ID of the interpolation case if a coupled thermal structural analysis is requested
      unsigned int interpolation_case_ID;
    };
}


#endif // __fesystem_linear_stress_solution_h__
