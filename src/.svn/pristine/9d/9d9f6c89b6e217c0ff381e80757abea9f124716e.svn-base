// $Id: FluidDisciplineInfo.C,v 1.1.6.1 2008-02-25 04:21:42 manav Exp $

// FESystem includes
#include "Discipline/StructuralDisciplineInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/FESystemSolverBase.h"
#include "Properties/PropertyCardParameter.h"


Discipline::StructuralDisciplineInfo::StructuralDisciplineInfo():
DisciplineInfo()
{
  
}



Discipline::StructuralDisciplineInfo::~StructuralDisciplineInfo()
{
}



std::istream& 
Discipline::StructuralDisciplineInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num=0;
  
  FESystemIOUtility::readFromInput(input, Discipline::STRUCTURAL_DISCIPLINE_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  FESystemIOUtility::readFromInput(input, "MESH_ID", this->mesh_ID);
  tag.clear();
  FESystemIOUtility::readFromInput(input, "ANALYSIS_TYPE", tag);
  this->analysis_type = Discipline::DisciplineAnalysisTypeEnum::enumID(tag);
  
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
  
  FESystemIOUtility::readFromInput(input, Discipline::STRUCTURAL_DISCIPLINE_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  
  return input;
}




std::istream& 
Discipline::operator >> (std::istream& input, 
                         Discipline::StructuralDisciplineInfo& info)
{
  info.readFromInputStream(input);
  return input;
}

