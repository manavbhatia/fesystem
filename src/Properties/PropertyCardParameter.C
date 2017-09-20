// $Id: PropertyCardParameter.C,v 1.4 2006-09-05 20:41:51 manav Exp $

// FESystem includes
#include "Properties/PropertyCardParameter.h"



Property::ParameterBase::ParameterBase(): 
ID(FESystemNumbers::InvalidID),
name(),
ref_value(0.0),
current_value(0.0)
{} 


Property::ParameterBase::~ParameterBase()
{
  
}




Property::LocalParameter::LocalParameter():
ParameterBase()
{
  
}


Property::LocalParameter::~LocalParameter()
{
  
}




Property::GlobalParameter::GlobalParameter():
ParameterBase()
{
  
}



Property::GlobalParameter::~GlobalParameter()
{
  
}



std::istream& 
Property::LocalParameter::readFromInputStream(std::istream& input)
{
  FESystemIOUtility::readFromInput(input,
                                   Property::LOCAL_PROPERTY_PARAMETER::name());
  
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  FESystemIOUtility::readFromInput(input, "NAME", this->name);
  
  std::string tag;
  FESystemIOUtility::readFromInput(input, "PARAMETER", tag);
  this->param_enum_ID = LocalParameterType::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, "REFERENCE_VALUE", this->ref_value);
  
  FESystemIOUtility::readFromInput(input,
                                   Property::LOCAL_PROPERTY_PARAMETER::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}



 
std::istream& 
Property::GlobalParameter::readFromInputStream(std::istream& input)
{
  FESystemIOUtility::readFromInput(input,
                                   Property::GLOBAL_PROPERTY_PARAMETER::name());
  
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  FESystemIOUtility::readFromInput(input, "NAME", this->name);
  
  FESystemIOUtility::readFromInput(input, "GLOBAL_PARAM_ID", 
                                   this->global_parameter_ID);
  
  FESystemIOUtility::readFromInput(input, "REFERENCE_VALUE", this->ref_value);
  
  FESystemIOUtility::readFromInput(input,
                                   Property::GLOBAL_PROPERTY_PARAMETER::name());

  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}





std::istream& 
Property::operator>>(std::istream& input, LocalParameter& param)
{
  param.readFromInputStream(input);
  
  return input;
}




std::istream& 
Property::operator>>(std::istream& input, GlobalParameter& param)
{
  param.readFromInputStream(input);
  
  return input;
}
  

