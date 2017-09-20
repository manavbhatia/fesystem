/*
 *  AeroelasticityAnalysisDriver.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/20/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */


// C++ includes


// FESystem includes
#include "AeroelasticityAnalysisDriverBase.h"
#include "Discipline/StructuralAnalysis.h"
#include "Solutions/AeroelasticitySolution.h"



Driver::AeroelasticityAnalysisDriverBase::AeroelasticityAnalysisDriverBase
(const unsigned int ID, FESystem::FESystemController& controller,
 unsigned int driver_enum_ID):
Driver::AnalysisDriver(ID, controller, driver_enum_ID),
current_flight_condition(NULL)
{
  
}


Driver::AeroelasticityAnalysisDriverBase::~AeroelasticityAnalysisDriverBase()
{
  
}


Discipline::AnalysisDisciplineBase& 
Driver::AeroelasticityAnalysisDriverBase::getStructuralDiscipline() const
{
  Assert(this->analysis_discipline_map.size()==2, ExcInternalError());
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::const_iterator it, end;
  it = this->analysis_discipline_map.find(Discipline::STRUCTURAL_DISCIPLINE::num());
  end = this->analysis_discipline_map.end();
  Assert(it!=end, ExcInternalError());
  
  return *(it->second);
}




Discipline::AnalysisDisciplineBase& 
Driver::AeroelasticityAnalysisDriverBase::getAerodynamicDiscipline() const
{
  Assert(this->analysis_discipline_map.size()==2, ExcInternalError());
  Assert(this->solution_base != NULL, ExcInternalError());
  
  Solution::AeroelasticitySolution* aeroelasticity_sol = 
  dynamic_cast<Solution::AeroelasticitySolution*> (this->solution_base);
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::const_iterator it, end;
  it = this->analysis_discipline_map.find(aeroelasticity_sol->getAerodynamicDisciplineEnumID());
  end = this->analysis_discipline_map.end();
  Assert(it!=end, ExcInternalError());
  
  return *(it->second);
}


const Loads::FlightCondition& 
Driver::AeroelasticityAnalysisDriverBase::getCurrentFlightCondition() const
{
  Assert(this->current_flight_condition != NULL, ExcInternalError());
  return *(this->current_flight_condition);
}

