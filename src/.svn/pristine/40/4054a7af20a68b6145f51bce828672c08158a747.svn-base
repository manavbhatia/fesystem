// $Id: ThermalDisciplineInfo.C,v 1.5 2006-11-13 10:11:57 manav Exp $

// FESystem includes
#include "Discipline/ThermalDisciplineInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/FESystemSolverBase.h"
#include "Properties/PropertyCardParameter.h"


Discipline::ThermalDisciplineInfo::ThermalDisciplineInfo():
DisciplineInfo()
{
  
}



Discipline::ThermalDisciplineInfo::~ThermalDisciplineInfo()
{
}



std::istream& 
Discipline::ThermalDisciplineInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0, rad_id = 0;
  
  FESystemIOUtility::readFromInput(input, Discipline::THERMAL_DISCIPLINE_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  FESystemIOUtility::readFromInput(input, "MESH_ID", this->mesh_ID);
  tag.clear();
  FESystemIOUtility::readFromInput(input, "ANALYSIS_TYPE", tag);
  this->analysis_type = Discipline::DisciplineAnalysisTypeEnum::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, "N_RADIATION_CAVITIES", num);
  this->radiation_cavity_IDs.resize(num);
  
  for (unsigned int i=0; i < num; i++)
    {
    input >> rad_id;
    (this->radiation_cavity_IDs)[i] = rad_id;
    }

  FESystemIOUtility::readFromInput(input, "N_LOCAL_PARAMETERS", num);
  unsigned int param_id = 0; 
  bool insert_success = false;
  for (unsigned int i=0; i < num; i++)
    {
    tag.clear();
    input >> tag;
    param_id = Property::LocalParameterType::enumID(tag);
    
    insert_success = this->local_parameters.insert(param_id).second;
    
    AssertThrow(insert_success, FESystemExceptions::ExcDuplicateID("Local Parameter", param_id));
    }
  
  FESystemIOUtility::readFromInput(input, Discipline::THERMAL_DISCIPLINE_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  
  return input;
}




std::istream& 
Discipline::operator >> (std::istream& input, 
                         Discipline::ThermalDisciplineInfo& info)
{
  info.readFromInputStream(input);
  return input;
}

