// $Id: PropertyCardBase.C,v 1.2 2006-09-05 20:41:51 manav Exp $


// FESystem includes
#include "Properties/PropertyCardBase.h"
#include "Utilities/InputOutputUtility.h"


Property::PropertyCardBase::PropertyCardBase():
cardID(FESystemNumbers::InvalidID),
local_params_initialized(false),
global_params_initialized(false),
function_database(NULL),
parameter_function_ID_table(NULL)
{
  this->parameter_function_ID_table = new Table<2,unsigned int>;
  this->parameter_function_value_table = new Table<3, double>;
}





Property::PropertyCardBase::~PropertyCardBase()
{
  // delete the function ID table
  delete this->parameter_function_ID_table; 
  this->parameter_function_ID_table = NULL;
  
  delete this->parameter_function_value_table;
  this->parameter_function_value_table = NULL;
}




void 
Property::PropertyCardBase::addLocalParameter(Property::LocalParameter* param)
{
  Assert(param != NULL,
         ExcEmptyObject());
  
  unsigned int param_ID = param->getID(),
    enum_ID = param->getLocalParameterEnumID(),
    local_ID = this->parameter_map.size();
  
  bool insert_success = false;
  
  // now insert the param in the param_map and the local_parameter map
  insert_success = this->parameter_map.insert
    (Property::PropertyCardBase::ParameterMap::value_type(local_ID, param)).second;
  
  // an error thrown if the insert failed, which should never happen, since the ID is 
  // always chosen to be unique
  Assert(insert_success, ExcInternalError());

  // now add this to the local parameter map
  insert_success = this->local_parameter_map.insert
    (Property::PropertyCardBase::LocalParameterMap::value_type(enum_ID, param)).second;
  
  // if this is true, then the local parameter has been specified a 2nd time.
  Assert(insert_success,
         Property::PropertyCardBase::ExcDuplicateParameterSpecified(enum_ID));

  // once this is done, insert the ID entries in the ID map. 
  insert_success = this->parameter_external_to_internal_ID_map.insert
    (Property::PropertyCardBase::IDMap::value_type(param_ID, local_ID)).second;
  
  Assert(insert_success,
         Property::PropertyCardBase::ExcDuplicateParameterSpecified(param_ID));
}




void 
Property::PropertyCardBase::addGlobalParameter(Property::GlobalParameter* param)
{
  Assert(param != NULL,
         ExcEmptyObject());
  
  unsigned int param_ID = param->getID(),
    global_ID = param->getGlobalParameterID(),
    local_ID = this->parameter_map.size();
  
  bool insert_success = false;
  
  // now insert the param in the param_map and the local_parameter map
  insert_success = this->parameter_map.insert
    (Property::PropertyCardBase::ParameterMap::value_type(local_ID, param)).second;
  
  // an error thrown if the insert failed, which should never happen, since the ID is 
  // always chosen to be unique
  Assert(insert_success, ExcInternalError());
  
  // now add this to the local parameter map
  insert_success = this->global_parameter_map.insert
    (Property::PropertyCardBase::GlobalParameterMap::value_type(global_ID, param)).second;
  
  // if this is true, then the local parameter has been specified a 2nd time.
  Assert(insert_success,
         Property::PropertyCardBase::ExcDuplicateParameterSpecified(global_ID));
  
  // once this is done, insert the ID entries in the ID map. 
  insert_success = this->parameter_external_to_internal_ID_map.insert
    (Property::PropertyCardBase::IDMap::value_type(param_ID, local_ID)).second;
  
  Assert(insert_success,
         Property::PropertyCardBase::ExcDuplicateParameterSpecified(param_ID));
  
}





std::istream& 
Property::PropertyCardBase::readParametersFromInputStream(std::istream& input)
{
  
  // read the local parameters
  {
    unsigned int n_params = 0;
    
    FESystemIOUtility::readFromInput(input, "N_LOCAL_PARAMETERS", n_params);
    
    for (unsigned int i=0; i<n_params; i++)
      {
      Property::LocalParameter *param = new Property::LocalParameter;
      input >> (*param);
      
      // insert the parameter in the map
      this->addLocalParameter(param);
      }
  }
  
  
  // now read the global parameters
  {
    unsigned int n_params = 0;
    
    FESystemIOUtility::readFromInput(input, "N_GLOBAL_PARAMETERS", n_params);
    
    for (unsigned int i=0; i<n_params; i++)
      {
      Property::GlobalParameter *param = new Property::GlobalParameter;
      input >> (*param);
      
      // insert the parameter in the map
      this->addGlobalParameter(param);
      }
  }
  
  return input;
}


