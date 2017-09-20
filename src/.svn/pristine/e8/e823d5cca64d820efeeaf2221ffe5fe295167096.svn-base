// $Id: ShapeParameter.C,v 1.3 2007-01-15 18:59:43 manav Exp $

// FESystem includes
#include "DesignData/ShapeParameter.h"
#include "Utilities/InputOutputUtility.h"
#include "Discipline/AnalysisDisciplineBase.h"


DesignData::ShapeParameter::ShapeParameter():
DesignData::DesignParameter()
{}


DesignData::ShapeParameter::~ShapeParameter()
{}



std::istream&
DesignData::ShapeParameter::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, DesignData::SHAPE_PARAMETER::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  FESystemIOUtility::readFromInput(input, "NAME", this->name);
  FESystemIOUtility::readFromInput(input, "VALUE", this->value);
  FESystemIOUtility::readFromInput(input, "SENSITIVITY_METHOD");
  FESystemIOUtility::readFromInput(input, "EULER_FD_SENSITIVITY");
  
  this->sensitivity_method = DesignData::EULER_FD_SENSITIVITY::num();
  FESystemIOUtility::readFromInput(input, "PERTURBATION_STEP_SIZE",
                                   this->finite_difference_step_size);
  unsigned int n_mesh = 0, mesh_ID =0, discipline_ID = 0;
  bool insert_success = false;
  FESystemIOUtility::readFromInput(input, "PERTURBED_MESH", n_mesh);
  for (unsigned int i=0; i < n_mesh; i++)
    {
    tag.clear();
    input >> tag;
    discipline_ID = Discipline::AnalysisDisciplineEnum::enumID(tag);
    input >> mesh_ID;
    
    insert_success = this->perturbed_mesh_ID_map.insert
      (DesignData::ShapeParameter::IDMap::value_type(discipline_ID, mesh_ID)).second;
    
    Assert(insert_success, ExcInternalError());
    }
  
  FESystemIOUtility::readFromInput(input, DesignData::SHAPE_PARAMETER::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}


std::istream& 
DesignData::operator>>(std::istream& input, DesignData::ShapeParameter& param)
{
  param.readFromInputStream(input);
  
  return input;
}


