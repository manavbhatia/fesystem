// $Id: LinearizedBucklingEigenSolution.h,v 1.3.6.2 2007-06-13 14:59:40 manav Exp $

#ifndef __fesystem_structural_linearized_buckling_eigen_solution_h__
#define __fesystem_structural_linearized_buckling_eigen_solution_h__


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

#ifndef LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID 
#define LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID 5
#else 
#error
#endif

#ifndef LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_NAME 
#define LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_NAME "LINEARIZED_BUCKLING_EIGEN_SOLUTION"
#else 
#error
#endif



namespace Solution
{
  
  DeclareEnumName(LINEARIZED_BUCKLING_EIGEN_SOLUTION, Solution::SolutionEnum,
                  LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_ID, 
                  LINEARIZED_BUCKLING_EIGEN_SOLUTION_ENUM_NAME);
  
  
  /// this class provides the base functionality for all solutions
  class LinearizedBucklingEigenSolution: public Solution::SolutionBase
    {
public:
      /// constructor
      LinearizedBucklingEigenSolution();
      
      /// destructor
      virtual ~LinearizedBucklingEigenSolution();
      
      /// this method will solve for all the load cases and sensitivity analysis
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

        /// boolean to check if coupled thermoelastic solution is needed, which will 
        /// need solution of the thermal discipline
        bool coupled_thermoelastic;
      
      /// ID of the interpolation case if a coupled thermal structural analysis is requested
      unsigned int interpolation_case_ID;
      
    };
}


#endif // __fesystem_structural_linearized_buckling_eigen_solution_h__
