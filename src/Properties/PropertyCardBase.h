// $Id: PropertyCardBase.h,v 1.6 2006-09-05 20:41:51 manav Exp $

#ifndef __fesystem_property_card_base_h__
#define __fesystem_property_card_base_h__

// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"
#include "Database/FunctionDatabase.h"
#include "Properties/PropertyCardParameter.h"

// deal.II includes
#include "base/table.h"




namespace Property
{
  
  /// The \p PropertyCard class that provides the necessary data structure 
  /// for storage of properties in general. This is a base class,
  /// and different material/structural property cards can derive from this. Each property card
  /// can have dependence on some parameter values, which is defined through the input. If a 
  /// parameter dependence has been defined, then before the parameter values can be accessed, 
  /// the card needs to be initialized at some parameter values. During reinitialization, if
  /// no value has been specified for a parameter, a default value for the same is assumed.
  /// Between any two consecutive initializations, the card needs to be cleared of any previous
  /// initialization. This is needed to make sure that the card is cleared of any remnant value
  /// initializations to avoid errors. Once initialized, the value of the property can be accessed from
  /// the card.
  
  
  class PropertyCardBase 
  {
public:	
    
    /// default constructor
    PropertyCardBase();
    
    /// destructor
    virtual ~PropertyCardBase();
    
    /// @returns the id of the card
    inline unsigned int getID() const;
    
    virtual unsigned int getNProperties() const = 0;
    

    inline unsigned int getNParameters() const;
    

    inline unsigned int getNLocalParameters() const;
    
    
    inline unsigned int getNGlobalParameters() const;
    
    /// @returns the property card kind enumeration ID
    virtual unsigned int getPropertyCardKindEnumID() const = 0;
    
    /// @returns the property card kind enumeration name
    virtual const std::string getPropertyCardKindEnumName() const = 0;
    
    /// @returns true if the card depends on the parameter
    bool checkLocalParameterDependence(const unsigned int enumID) const;
    
    /// @returns true if the card depends on the parameter
    bool checkGlobalParameterDependence(const unsigned int param_global_ID) const;
    
    /// attaches the function database to this card
    inline void attachFunctionDatabase(FunctionDatabase &func_database);
    
    
    /// this method initializes the property card at the values specified for each 
    /// parameter. The input to this method is a map of the parameter and its value. 
    /// For parameters whose values have not been specified in the input map, 
    /// default values are used.
    void reinitLocalParameters(const std::map<unsigned int, double>* value_map = NULL);
    
    /// this method initializes the property card at the values specified for each 
    /// parameter. This method will update all values from the input values in the parameter
    void reinitGlobalParameters(const std::map<unsigned int, double>& value_map);
    
    
    /// this method clears any previous initialization. This needs to be called before any 
    /// reinitialization. Before re-intialization of any parameter, its value should be cleared
    /// using this function. 
    inline void clearLocalParameterInitialization();
    
    /// this method clears the function value initialization for local params specified 
    /// in the input. No other params are cleared.
    inline void
      partialClearLocalParameterInitialization(std::vector<unsigned int>& params_to_clear);
      
    
    /// this method clears any previous initialization. This needs to be called before any 
    /// reinitialization. Before re-intialization of any parameter, its value should be cleared
    /// using this function. If no input is provided, then all the values are cleared. Otherwise, 
    /// only the values listed in this vector are cleared.
    inline void clearGlobalParameterInitialization();
    
    
    
protected:
      
      
    /// @returns parameter internal ID for the external ID specifies in the input. 
    /// This external 
    /// ID is the ID of the parameter as specified in the input for the property card
    inline unsigned int getParameterInternalID(const unsigned int exernal_ID) const;
    
    /// @returns parameter internal ID for a global parameter. 
    /// The ID specified in the input
    /// is the global ID of the parameter in the design database. 
    inline unsigned int 
      getGlobalParameterInternalIDForGlobalID(const unsigned int external_ID) const;

    /// @returns parameter internal ID for a local parameter. 
    /// The ID specified in the input is the enum ID of the parameter 
    inline unsigned int 
      getLocalParameterInternalIDForEnumID(const unsigned int enum_ID) const;
    
    
    /// @returns property internal ID for the external ID specifies in the input
    inline unsigned int getPropertyInternalID(const unsigned int enum_ID) const;

    /// reinitializes the function values and the function derivatives for the parameter IDs specified
    /// in the input
    void reinitParameterFunctionValues(const std::vector<unsigned int>& param_ids);
    
    /// adds a global parameter to the card
    void addLocalParameter(LocalParameter* param);
    
    /// adds a local parameter to the card
    void addGlobalParameter(GlobalParameter* param);    
    
    /// reads parameters from the input
    std::istream& readParametersFromInputStream(std::istream& input);
    
    ///property card ID
    unsigned int cardID;
    
    /// if the local params are initialized or not
    bool local_params_initialized;
    
    /// if the global params are initialized or not
    bool global_params_initialized;
    
    /// a pointer to the function database that stores the functions defining the
    /// property-parameter relations
    FunctionDatabase *function_database;
    
    /// a table to store the IDs of the functions for definition of parameter dependence
    /// of each property. After data has been read into this card, if a paricular entry is still
    /// equal to the invalidID, then that implies that the property is not dependent on that
    /// parameter. Else, it is dependent, and the entry specifies the function defining the
    /// dependence
    Table<2, unsigned int>* parameter_function_ID_table;
    
    /// a table to store the function values and sensitivity for the parameters and property
    Table<3, double>* parameter_function_value_table;
    
    /// an internal error to define non existence of property
    DeclException1(ExcPropertyDoesNotExist, std::string, 
                   << "Requested Property: " << arg1 
                   << "  does not exist in Card.");
    
    /// an internal error to define non existence of property
    DeclException1(ExcDuplicatePropertySpecified, std::string, 
                   << "Specified Property: " << arg1 
                   << "  already exists in Card.");
    
    /// an internal error if the specified parameter does not exist in the 
    /// card
    DeclException1(ExcParameterDoesNotExist, std::string, 
                   << "Specified Parameter: " << arg1 
                   << "  does not exist in Card.");
    
    /// an internal error if the specified parameter does not exist in the 
    /// card
    DeclException1(ExcDuplicateParameterSpecified, unsigned int, 
                   << "Specified Parameter ID: " << arg1 
                   << "  already exist in Card.");
    
    
    
    /// an internal error if the specified parameter does not exist in the 
    /// card
    DeclException1(ExcParameterIDDoesNotExist, unsigned int, 
                   << "Specified Parameter ID: " << arg1 
                   << "  does not exist in Card.");
    
    /// an internal error if the global parameter is not specified in the input
    DeclException1(ExcGlobalParameterNotSpecified, unsigned int, 
                   << "Global Parameter: " << arg1 
                   << "  not specified in input.");
    
    
    
    /// this is a local type definition for a more readable code
    typedef std::map<unsigned int, Property::ParameterBase*>  ParameterMap;
    typedef std::map<unsigned int, Property::LocalParameter*>  LocalParameterMap;
    typedef std::map<unsigned int, Property::GlobalParameter*>  GlobalParameterMap;
    typedef std::map<unsigned int, unsigned int> IDMap;
    
    
    /// The map for storing the pointers to parameters. The key used is the internal ID
    /// of the parameter. The internal IDs start from 0.
    ParameterMap parameter_map;

    /// map for storing the pointers to local parameters. These pointers, are stored for 
    /// convenience. This map uses the parameter enum ID as key
    LocalParameterMap local_parameter_map;
    
    /// map for storing the pointers to local parameters. These pointers, are stored for 
    /// convenience. This map uses the global ID of the parameter as key
    GlobalParameterMap global_parameter_map;

    /// external to internal ID map for the parameters
    IDMap parameter_external_to_internal_ID_map;
    
    /// enumID to internal ID map for the properties
    IDMap property_enumID_to_internal_ID_map;
    
    /// set of local param IDs that have been set to their default values
    std::set<unsigned int> default_local_params;
  };
  
}




//
// ****************************  function definitions for the base class **************************
//



inline unsigned int 
Property::PropertyCardBase::getID() const
{
  Assert(this->cardID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(this->cardID));
	
  return this->cardID;
}



inline 
unsigned int 
Property::PropertyCardBase::getNParameters() const
{
  return this->parameter_map.size();
}


inline
unsigned int 
Property::PropertyCardBase::getNLocalParameters() const
{
  return this->local_parameter_map.size();
}



inline
unsigned int 
Property::PropertyCardBase::getNGlobalParameters() const
{
  return this->global_parameter_map.size();
}




inline void
Property::PropertyCardBase::
attachFunctionDatabase(FunctionDatabase &func_database)
{
  this->function_database = &func_database;
}





inline void
Property::PropertyCardBase::
reinitLocalParameters(const std::map<unsigned int, double>* value_map)
{
  // before calling this function, the user should have cleared previous initializations
  Assert(this->local_params_initialized == false,
         ExcInvalidState());
  
  
  static std::vector<unsigned int> reinitialized_param_ids;
  reinitialized_param_ids.clear();
  
  // set the values of the parameters, and calculate and store the function values 
  // and their derivatives. 
  // a list of all function values that are set to their default values will have be 
  // stored so that there is no need to initialize them or clear them when a partial 
  // clear is performed 
  Property::PropertyCardBase::LocalParameterMap::iterator 
  param_it, param_end;
  param_it = this->local_parameter_map.begin();
  param_end = this->local_parameter_map.end();
  
  double value = 0.0;
  
  bool insert_success = false;
  // if the map is null, then all params will have their default values
  if (value_map == NULL)
    {
    for (; param_it != param_end; param_it++)
      {
      // if the local parameter has not already been initialized to its local value, 
      // initialize it
      if (this->default_local_params.count(param_it->second->getID()) == 0)
        {
        value = param_it->second->getReferenceValue();
        param_it->second->setCurrentValue(value);
        
        // also add this param to the set of params that have been set to their default 
        // values.
        insert_success = this->default_local_params.insert(param_it->second->getID()).second;
        Assert(insert_success, FESystemExceptions::ExcDuplicateID("Local Parameter", 
                                                                  param_it->second->getID()));
        
        reinitialized_param_ids.push_back(param_it->second->getID());
        }
      }
    }
  else
    {
    
    // now set the values of the parameters by syncing the parameters that have been 
    // specified in the input
    unsigned int param_count_in_local_default_param_map = 0, 
    param_enum_id = 0;
    
    std::map<unsigned int, double>::const_iterator input_param_it, input_param_end;
    input_param_end = value_map->end();
    
    for (; param_it != param_end; param_it++)
      {
      param_enum_id = param_it->first;
      input_param_it = value_map->find(param_enum_id);
      
      // check if this has been set to its default value
      param_count_in_local_default_param_map = 
        this->default_local_params.count(param_it->second->getID());
      
      if (input_param_it == input_param_end)
        {
        // if this param has not already been initialized to its default value, 
        // then do the same, and add this ID to the set.
        // otherwise, if it has already been initialized to its local value, then
        // no need to do anything
        if (param_count_in_local_default_param_map == 0)
          {
          value = param_it->second->getReferenceValue();

          // also add this param to the set of params that have been set to their default 
          // values.
          insert_success = this->default_local_params.insert(param_it->second->getID()).second;
          Assert(insert_success, FESystemExceptions::ExcDuplicateID("Local Parameter", 
                                                                    param_it->second->getID()));
          }
        else
          continue;
        }
      else
        {
        // we need to set the value of this param to the one specified in the input. 
        // for this, we need to check if the param has been previously set to its 
        // default value. If it has been, then its ID needs to be cleared from the 
        // local default param set
        if (param_count_in_local_default_param_map != 0)
          this->default_local_params.erase(param_it->second->getID());
        
        value = input_param_it->second;
        }
      
      param_it->second->setCurrentValue(value);
      reinitialized_param_ids.push_back(param_it->second->getID());
      }
    
    }
  
  
  // once the values have been set, the functions need to be initialized
  // call the function value initialization for the parameter IDs
  // that have just been initialized
  this->reinitParameterFunctionValues(reinitialized_param_ids);
  
  // finally, set the tag for initialization of the local params
  this->local_params_initialized = true;
}




inline void
Property::PropertyCardBase::
reinitGlobalParameters(const std::map<unsigned int, double>& value_map)
{
  // before calling this function, the user should have cleared previous initializations
  Assert(this->global_params_initialized == false,
         ExcInvalidState());
  
  // iterate over all the params in this card and set the current values of the params
  // that are specified in the map
  PropertyCardBase::GlobalParameterMap::iterator 
    param_it, param_end;
  param_it = this->global_parameter_map.begin();
  param_end = this->global_parameter_map.end();
  
  double value = 0.0;
  unsigned int param_id = 0;
  
  std::map<unsigned int, double>::const_iterator input_param_it, input_param_end;
  input_param_end = value_map.end();
  
  for (; param_it != param_end; param_it++)
    {
    param_id = param_it->first;
    input_param_it = value_map.find(param_id);

    Assert(input_param_it != input_param_end,
           ExcGlobalParameterNotSpecified(param_id));

    value = input_param_it->second;
    param_it->second->setCurrentValue(value);
    }
  
  // once the values have been set, the functions need to be initialized
  // call the function function value initialization for the parameter IDs
  // that have just been initialized
  static std::vector<unsigned int> 
    reinitialized_param_ids(this->global_parameter_map.size());
  param_it = this->global_parameter_map.begin();
  
  unsigned int n_param = 0;
  for (; param_it != param_end; param_it++)
    {
    reinitialized_param_ids[n_param] = param_it->second->getID();
    n_param++;
    }
  
  this->reinitParameterFunctionValues(reinitialized_param_ids);
  
  // finally, set the tag for initialization of the local params
  this->global_params_initialized = true;
}





inline void
Property::PropertyCardBase::
reinitParameterFunctionValues(const std::vector<unsigned int>& param_ids)
{
  Assert(this->function_database != NULL,
         ExcInvalidState());
  
  
  std::vector<unsigned int>::const_iterator param_it, param_end;
  param_it = param_ids.begin();
  param_end = param_ids.end();
  
  unsigned int func_id = 0, param_internal_id = 0, prop_it=0;
  const unsigned int  n_prop = this->getNProperties();
  double param_value = 0.0, func_value = 0.0, func_derivative = 0.0;
  
  for ( ; param_it != param_end; param_it++)
    {
    param_internal_id = this->getParameterInternalID(*param_it);
    param_value = this->parameter_map[param_internal_id]->getCurrentValue();
    
    for (prop_it=0; prop_it < n_prop; prop_it++)
      {
      // get the function ID and then the function values and its derivative wrt the 
      // parameter
      func_id = (*this->parameter_function_ID_table)(prop_it, param_internal_id);

      if (func_id == FESystemNumbers::InvalidID)
        continue;
      
      this->function_database->getScalarFunctionValue(func_id, param_value, func_value);
      this->function_database->getScalarFunctionDerivative(func_id, param_value, func_derivative);
      
      // not enter the function values in the table
      (*this->parameter_function_value_table)(prop_it, param_internal_id, 0) = func_value;
      (*this->parameter_function_value_table)(prop_it, param_internal_id, 1) = func_derivative;
      }
    }
}





inline void
Property::PropertyCardBase::
partialClearLocalParameterInitialization(std::vector<unsigned int>& params_to_clear)
{
  // iterate over all the params specified and clear their initialization
  std::vector<unsigned int>::const_iterator it, end;
  it = params_to_clear.begin();
  end = params_to_clear.end();
  
  for (; it != end; it++)
    this->default_local_params.erase(*it);
  
  this->local_params_initialized = false;
}




inline void
Property::PropertyCardBase::
clearLocalParameterInitialization()
{
  this->local_params_initialized = false;
}
  




inline void
Property::PropertyCardBase::
clearGlobalParameterInitialization()
{
  this->global_params_initialized = false;
}






inline 
unsigned int 
Property::PropertyCardBase::
getParameterInternalID(const unsigned int external_ID) const
{
  // make sure that this given ID is valid
  Assert(external_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(external_ID));
  
  // now get the iterators to the ID, and make sure that the ID exists in the card
  
  Property::PropertyCardBase::IDMap::const_iterator it, end;
  it = this->parameter_external_to_internal_ID_map.find(external_ID);
  end = this->parameter_external_to_internal_ID_map.end();
  
  Assert(it != end,
         Property::PropertyCardBase::ExcParameterIDDoesNotExist(external_ID));
  
  return it->second;
}




inline 
unsigned int 
Property::PropertyCardBase::
getGlobalParameterInternalIDForGlobalID(const unsigned int global_ID) const
{
  // make sure that this given ID is valid
  Assert(global_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(global_ID));
  
  // now get the iterators to the ID, and make sure that the ID exists in the card
  
  Property::PropertyCardBase::GlobalParameterMap::const_iterator it, end;
  it = this->global_parameter_map.find(global_ID);
  end = this->global_parameter_map.end();
  
  Assert(it != end,
         Property::PropertyCardBase::ExcParameterIDDoesNotExist(global_ID));
  
  unsigned int external_ID = it->second->getID();
  
  return this->getParameterInternalID(external_ID);
}






inline 
unsigned int 
Property::PropertyCardBase::
getLocalParameterInternalIDForEnumID(const unsigned int enum_ID) const
{
  // make sure that this given ID is valid
  Assert(enum_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(enum_ID));
  
  // now get the iterators to the ID, and make sure that the ID exists in the card
  Property::PropertyCardBase::LocalParameterMap::const_iterator it, end;
  it = this->local_parameter_map.find(enum_ID);
  end = this->local_parameter_map.end();
  
  Assert(it != end,
         Property::PropertyCardBase::ExcParameterIDDoesNotExist(enum_ID));

  unsigned int external_ID = it->second->getID();
  
  return this->getParameterInternalID(external_ID);
}






inline 
unsigned int 
Property::PropertyCardBase::
getPropertyInternalID(const unsigned int external_ID) const
{
  // make sure that this given ID is valid
  Assert(external_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(external_ID));
  
  // now get the iterators to the ID, and make sure that the ID exists in the card
  Property::PropertyCardBase::IDMap::const_iterator it, end;
  it = this->property_enumID_to_internal_ID_map.find(external_ID);
  end = this->property_enumID_to_internal_ID_map.end();
  
  Assert(it != end,
         Property::PropertyCardBase::ExcParameterIDDoesNotExist(external_ID));
  
  return it->second;
}


inline bool
Property::PropertyCardBase::checkLocalParameterDependence(const unsigned int enumID) const
{
  // to check for dependence on the specified local param, look in the map
  if (this->local_parameter_map.count(enumID) ==0)
    return false;
  else 
    return true;
}



inline bool
Property::PropertyCardBase::checkGlobalParameterDependence(const unsigned int param_global_ID) const
{
  // to check for dependence on the specified global param, look in the map
  if (this->global_parameter_map.count(param_global_ID) ==0)
    return false;
  else 
    return true;
}




#endif // #define __fesystem_property_card_base_h__
