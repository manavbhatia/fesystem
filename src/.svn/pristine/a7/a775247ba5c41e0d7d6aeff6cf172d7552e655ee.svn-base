// $Id: StructuralVibrationEigenSolution.h,v 1.6.6.2 2008-08-25 04:39:11 manav Exp $

#ifndef __fesystem_structural_vibration_eigen_solution_h__
#define __fesystem_structural_vibration_eigen_solution_h__


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

#ifndef STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID 
#define STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID 4
#else 
#error
#endif

#ifndef STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_NAME 
#define STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_NAME "STRUCTURAL_VIBRATION_EIGEN_SOLUTION"
#else 
#error
#endif



namespace Solution
{
  
  DeclareEnumName(STRUCTURAL_VIBRATION_EIGEN_SOLUTION, Solution::SolutionEnum,
                  STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_ID, 
                  STRUCTURAL_VIBRATION_EIGEN_SOLUTION_ENUM_NAME);
  
  
  /// this class provides the base functionality for all solutions
  class StructuralVibrationEigenSolution: public Solution::SolutionBase
    {
public:
      /// constructor
      StructuralVibrationEigenSolution();
      
      /// destructor
      virtual ~StructuralVibrationEigenSolution();
      
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
      
      /// @returns true if the geometric effects are to be included
      bool ifIncludeGeometricEffects() const;
            
protected:
      
//      /// solution w/o geometric effects
//      void solveWithOutGeometricEffects();
//      
//      /// solution w/o geometric effects
//      void solveWithGeometricEffects();
      
      /// if the geometric effects should be included in the solution
      bool include_geometric_effects;

      /// boolean to check if coupled thermoelastic solution is needed, which will 
      /// need solution of the thermal discipline
      bool coupled_thermoelastic;
      
      /// ID of the interpolation case if a coupled thermal structural analysis is requested
      unsigned int interpolation_case_ID;
      
    };
}


inline
bool 
Solution::StructuralVibrationEigenSolution::ifIncludeGeometricEffects() const
{
  return this->include_geometric_effects;
}


#endif // __fesystem_structural_vibration_eigen_solution_h__
