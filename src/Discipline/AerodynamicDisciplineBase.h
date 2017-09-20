/*
 *  AerodynamicDisciplineBase.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */


#ifndef __fesystem_aerodynamic_discipline_base_h__
#define __fesystem_aerodynamic_discipline_base_h__


// C++ includes


// FESystem includes
#include "Discipline/AnalysisDisciplineBase.h"
#include "Loads/load.h"

#ifndef AERODYNAMIC_OPERATION_CONDITION_ENUM_ID
#define AERODYNAMIC_OPERATION_CONDITION_ENUM_ID 12
#else
#error
#endif

#ifndef AERODYNAMIC_OPERATION_CONDITION_ENUM_NAME
#define AERODYNAMIC_OPERATION_CONDITION_ENUM_NAME "AERODYNAMIC_OPERATION_CONDITION"
#else
#error
#endif



DeclareEnumName(AERODYNAMIC_OPERATION_CONDITION, LoadNameEnum,
                AERODYNAMIC_OPERATION_CONDITION_ENUM_ID,
                AERODYNAMIC_OPERATION_CONDITION_ENUM_NAME);



namespace Discipline {

  // forward decleration
  class DisciplineInfo;
  
  /// base class of the panel methods
  class AerodynamicDisciplineBase: public AnalysisDisciplineBase
    {
    public:      
      /// constructor
      AerodynamicDisciplineBase(FESystem::FESystemController& controller,
                                const Discipline::DisciplineInfo& info);
      
      /// destructor
      virtual ~AerodynamicDisciplineBase();
      
      
    protected:
      
      
    };
  
}



#endif // __fesystem_aerodynamic_discipline_base_h__
