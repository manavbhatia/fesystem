// $Id: TransientStructuralSolution.h,v 1.1.2.1 2007-09-06 07:14:15 manav Exp $

#ifndef __fesystem_transient_structural_solution_h__
#define __fesystem_transient_structural_solution_h__


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

#ifndef TRANSIENT_STRUCTURAL_SOLUTION_ENUM_ID 
#define TRANSIENT_STRUCTURAL_SOLUTION_ENUM_ID 8
#else 
#error
#endif

#ifndef TRANSIENT_STRUCTURAL_SOLUTION_ENUM_NAME 
#define TRANSIENT_STRUCTURAL_SOLUTION_ENUM_NAME "TRANSIENT_STRUCTURAL_SOLUTION"
#else 
#error
#endif



namespace Solution
{
  
  DeclareEnumName(TRANSIENT_STRUCTURAL_SOLUTION, Solution::SolutionEnum,
                  TRANSIENT_STRUCTURAL_SOLUTION_ENUM_ID, TRANSIENT_STRUCTURAL_SOLUTION_ENUM_NAME);
  
  
  /// this class provides the base functionality for all solutions
  class TransientStructuralSolution: public Solution::SolutionBase
    {
public:
      /// constructor
      TransientStructuralSolution();
      
      /// destructor
      virtual ~TransientStructuralSolution();
      
      
      /// solves the problem
      virtual void solve();
      
      /// @returns the transient nature of the given discipline
      virtual unsigned int getDisciplineTransientNatureEnumID(const unsigned int discipline) const ;
      
      /// @returns the vector of SolNameInfo structures for the nodal solutions, and for the
      /// specified disciplines
      virtual FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
        getTimeIndependentNodalSolutionDataInfo(const unsigned int discipline_enum_ID);
      
      /// @returns the vector of SolNameInfo structures for the nodal solutions, and for the
      /// specified disciplines
      virtual FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
        getTimeDependentNodalSolutionDataInfo(const unsigned int discipline_enum_ID);
      
      /// @this will read from the input stream
      virtual std::istream& readFromInputStream(std::istream& input);
      
protected:
        
    };
}


#endif // __fesystem_transient_structural_solution_h__
