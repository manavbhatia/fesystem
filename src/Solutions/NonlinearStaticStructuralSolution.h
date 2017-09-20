// $Id: LinearStressSolution.C,v 1.4.6.3 2007-06-13 14:59:40 manav Exp $

/*
 *  NonlinearStaticStructuralSolution.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 11/21/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

#ifndef __fesystem_nonlinear_static_structural_solution_h__
#define __fesystem_nonlinear_static_structural_solution_h__


// C++ includes 


// FESystem includes


// libmesh includes

// forward declerations
namespace FESystem
{
  class FESystemController;
}

#ifndef NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_ID 
#define NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_ID 10 
#else  
#error 
#endif 

#ifndef NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_NAME 
#define NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_NAME "NONLINEAR_STATIC_STRUCTURAL_SOLUTION"
#else 
#error
#endif



namespace Solution
{
  
  DeclareEnumName(NONLINEAR_STATIC_STRUCTURAL_SOLUTION, Solution::SolutionEnum,
                  NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_ID, 
                  NONLINEAR_STATIC_STRUCTURAL_SOLUTION_ENUM_NAME);
  
  
  /// this class provides the base functionality for all solutions
  class NonlinearStaticStructuralSolution: public Solution::SolutionBase
  {
  public:
    /// constructor
    NonlinearStaticStructuralSolution();
    
    /// destructor
    virtual ~NonlinearStaticStructuralSolution();
    
    
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





#endif // __fesystem_nonlinear_static_structural_solution_h__

