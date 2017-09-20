// $Id:$
/*
 *  AeroelasticityAnalysisDriver.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/20/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */


#ifndef __fesystem_aeroelasticity_analysis_driver_base_h__
#define __fesystem_aeroelasticity_analysis_driver_base_h__


// C++ includes


// FESystem includes
#include "AnalysisDriver/AnalysisDriver.h"


namespace FESystem
{
  class FESystemController;
}


namespace Discipline
{
  class AnalysisDisciplineBase;
}

namespace Loads {
  class FlightCondition;
}

namespace Driver {

  
  class AeroelasticityAnalysisDriverBase: public AnalysisDriver
    {
    public:
      
      /// constructor
      AeroelasticityAnalysisDriverBase(const unsigned int ID,
                                       FESystem::FESystemController& controller,
                                       unsigned int driver_enum_ID);

      /// destructor
      virtual ~AeroelasticityAnalysisDriverBase() = 0;
            
      /// @returns structural discipline
      Discipline::AnalysisDisciplineBase& getStructuralDiscipline() const;
      
      /// @returns aerodynamic discipline
      Discipline::AnalysisDisciplineBase& getAerodynamicDiscipline() const;
      
      /// @returns the current flight condition
      const Loads::FlightCondition& getCurrentFlightCondition() const;
      
    protected:
            
      /// current flight condition
      const Loads::FlightCondition* current_flight_condition;
      
    };
}



#endif // __fesystem_aeroelasticity_analysis_driver_base_h__

