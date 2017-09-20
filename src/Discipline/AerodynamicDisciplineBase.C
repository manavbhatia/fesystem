/*
 *  AerodynamicDisciplineBase.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// C++ includes 


// FESystem includes
#include "AerodynamicDisciplineBase.h"




Discipline::AerodynamicDisciplineBase::AerodynamicDisciplineBase(FESystem::FESystemController& controller,
                                                                 const Discipline::DisciplineInfo& info):
Discipline::AnalysisDisciplineBase(controller, info)
{
  
}



Discipline::AerodynamicDisciplineBase::~AerodynamicDisciplineBase()
{
  
}
