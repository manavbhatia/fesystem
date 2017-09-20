/*
 *  AeroelasticitySolution.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// $Id:$

#ifndef __fesystem_aeroelasticity_solution_h__
#define __fesystem_aeroelasticity_solution_h__

// C++ includes
#include "Solutions/SolutionBase.h"
#include "geom/point.h" 

// forward declerations
namespace FESystem
{
  class FESystemController;
}

namespace Loads {
  class FlightCondition;
}

#ifndef AEROELASTICITY_SOLUTION_ENUM_ID 
#define AEROELASTICITY_SOLUTION_ENUM_ID 9
#else 
#error
#endif

#ifndef AEROELASTICITY_SOLUTION_ENUM_NAME 
#define AEROELASTICITY_SOLUTION_ENUM_NAME "AEROELASTICITY_SOLUTION"
#else 
#error
#endif


namespace Solution {
  
  
  DeclareEnumName(AEROELASTICITY_SOLUTION, Solution::SolutionEnum,
                  AEROELASTICITY_SOLUTION_ENUM_ID, 
                  AEROELASTICITY_SOLUTION_ENUM_NAME);
  
  /// This is the solution class of the aeroelasticity solution using the panel methods.
  class AeroelasticitySolution : public SolutionBase {
  public:
    
    /// constructor
    AeroelasticitySolution();
    
    
    /// destructor
    virtual ~AeroelasticitySolution();
    
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
        
    /// returns the analysis driver enumeration ID to be used in the solution
    unsigned int getSolutionMethodEnumID() const;
    
    /// returns the analysis driver enumeration name to be used in the solution
    unsigned int getSolutionMethodEnumName() const;

    /// @returns the number of lag terms for this solution
    virtual unsigned int numberOfAerodynamicLagTerms() const;
    
    /// @returns the enumeration ID of the kind of aeroelastic analysis driver
    unsigned int getAeroelasticDriverEnumID() const;
    
    /// @returns the enumeration ID of the kind of aerodynamic discipline modeling to be used
    unsigned int getAerodynamicDisciplineEnumID() const;

    /// @returns a reference to the vector of flight conditions for which aeroelastic analysis 
    /// needs to be performed 
    const std::vector<Loads::FlightCondition>& getFlightConditionVector() const; 

  protected:
    
    /// the solution sequeunce for a generaic aeroelastic eigenvalue analyis
    void aeroelasticSolve();
    
    /// whether modal order reduction is to be used
    bool use_modal_order_reduction;
    
    /// whether coupled thermal structural analysis is to be performed (unidirectional coupling)
    bool coupled_thermoelastic;
    
    /// whether thermal geometric effects should be included
    bool include_geometric_effects;
    
    /// linear solution completed
    bool completed_linear_solve;
    
    /// interpolation case ID for coupled thermal structural solution
    unsigned int interpolation_case_ID;
    
    /// number of lag terms
    unsigned int n_lag_terms;
    
    /// aerodynamic discipline enumeration ID
    unsigned int aerodynamic_discipline_enum_ID;

    /// aeroelastic method enumeration ID
    unsigned int aeroelastic_driver_enum_ID;
    
    /// vector of dynamic pressures
    std::vector<Loads::FlightCondition> flight_condition_vector;
  };
}




#endif // __fesystem_aeroelasticity_solution_h__

