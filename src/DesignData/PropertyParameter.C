// $Id: PropertyParameter.C,v 1.2.6.1 2008-08-21 20:42:16 manav Exp $


// FESystem includes
#include "DesignData/PropertyParameter.h"
#include "Utilities/InputOutputUtility.h"


DesignData::PropertyParameter::PropertyParameter():
DesignData::DesignParameter()
{}


DesignData::PropertyParameter::~PropertyParameter()
{}


std::istream&
DesignData::PropertyParameter::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, DesignData::PROPERTY_PARAMETER::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  FESystemIOUtility::readFromInput(input, "NAME", this->name);
  FESystemIOUtility::readFromInput(input, "VALUE", this->value);
  FESystemIOUtility::readFromInput(input, "SENSITIVITY_METHOD", tag);
  
  this->sensitivity_method = DesignData::SensitivityMethodEnum::enumID(tag);
  
  if (this->sensitivity_method == DesignData::EULER_FD_SENSITIVITY::num())
    {
    FESystemIOUtility::readFromInput(input, "PERTURBATION_STEP_SIZE", 
                                     this->finite_difference_step_size);
    }
  
  FESystemIOUtility::readFromInput(input, DesignData::PROPERTY_PARAMETER::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}


