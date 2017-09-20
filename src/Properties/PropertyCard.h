// $Id: PropertyCard.h,v 1.4.6.1 2007-06-13 14:58:48 manav Exp $

#ifndef __fesystem_property_card_h__
#define __fesystem_property_card_h__

// C++ includes
#include <iostream>
#include <string>
#include <map>

// FESystem includes
#include "Properties/PropertyCardBase.h"


namespace Property
{
  
  /// this is a general prameter, which has no implementation. It should not be used, 
  /// but only the specializations should be used.
  template<typename DataType>
  class PropertyCard: public Property::PropertyCardBase
  {
  };
  
  
  /// this is a specialization of the class for property values that are stored as double
  template < >
    class PropertyCard<double>: 
    public Property::PropertyCardBase
    {
public:
      inline PropertyCard();
      
      inline ~PropertyCard();
      
      virtual inline unsigned int getNProperties() const;
      
      
      /// this method returns the value of a property 
      virtual inline void getPropertyValue(const unsigned int enum_id,
                                           double& value) const;
      
      
      /// this method returns the derivative of a property value with respect to a specified 
      /// global parameter, whose global ID in the design database is specified
      /// the parameter values at which the derivatives should be calculated are obtained from the 
      /// initialization of the property card.
      virtual inline void getPropertyDerivativeForGlobalParameter(const unsigned int enum_ID,
                                                                  const unsigned int param_ID,
                                                                  double& value) const;
      
      
      /// this method returns the derivative of a property value with respect to a specified 
      /// local parameter, whose enum ID is specified
      /// the parameter values at which the derivatives should be calculated are obtained from the 
      /// initialization of the property card.
      virtual inline void getPropertyDerivativeForLocalParameter(const unsigned int enum_ID,
                                                                 const unsigned int param_ID,
                                                                 double& value) const;
      
protected:
        
        /// this method returns the derivative of a property value with respect to a 
        /// specified parameter.
        /// the parameter values at which the derivatives should be calculated is
        /// obtained from the 
        /// initialization of the property card. The parameter ID is the external
        /// ID of the parameter
        /// in the card. This method should not be used by the user, since the external IDs are 
        /// particular only to the card. The user can access the public methods for global and 
        /// local parameters
        //      template<typename ValueType>
        virtual inline void getPropertyDerivative(const unsigned int enum_ID,
                                                  const unsigned int param_ID,
                                                  double& value) const;
      
      
      
      
      /// @returns the property name from enumeration ID of the property
      virtual std::string getPropertyEnumName(const unsigned int enum_ID) const = 0;
      
      /// @returns the property ID from enumeration name of the property
      virtual const unsigned int getPropertyEnumID(const std::string& enum_name) const = 0;
      
      /// this method reads from the input stream
      inline std::istream& readFromInputStream(std::istream& input);
      
      typedef std::map<unsigned int, double> PropertyMap;
      
      PropertyMap property_reference_values;
    };
  
}



inline
Property::PropertyCard<double>::PropertyCard():
PropertyCardBase()
{
  
}



inline
Property::PropertyCard<double>::~PropertyCard()
{
  // delete the parameters since they were created using the new operator
  Property::PropertyCardBase::ParameterMap::iterator param_it, param_end;
  param_it = this->parameter_map.begin();
  param_end = this->parameter_map.end();
  
  for (; param_it != param_end; param_it++)
    {
    delete param_it->second;
    param_it->second = NULL;
    }
}




inline 
unsigned int 
Property::PropertyCard<double>::getNProperties() const
{
  return this->property_reference_values.size();
}





inline
void
Property::PropertyCard<double>::
getPropertyValue(const unsigned int enumID, double& value) const
{
  // check if the param has been initialized
  Assert (this->local_params_initialized == true,
          ExcInvalidState());
  
  Assert (this->global_params_initialized == true,
          ExcInvalidState());
  
  // make sure that the property name asked for, exists, and if it
  // does, get the value for same
  Property::PropertyCard<double>::PropertyMap::const_iterator 
    property_it = this->property_reference_values.find(enumID);
  
  static std::string name;
  name.clear(); name = this->getPropertyEnumName(enumID);
  
  AssertThrow(property_it != this->property_reference_values.end(), 
         Property::PropertyCardBase::ExcPropertyDoesNotExist
         (name));
  
  // now get the value, and add the contributions from the parameter functions
  //  static DataType prop_value;
  value = property_it->second;
  
  const unsigned int n_params = this->parameter_map.size(),
    prop_int_ID = this->getPropertyInternalID(enumID);
  
  double func_value = 0.0;
  
  for (unsigned int param_it = 0; param_it < n_params; param_it++)
    {
    func_value = (*this->parameter_function_value_table)(prop_int_ID, param_it, 0);
    value += func_value;
    }
}



inline
void
Property::PropertyCard<double>::getPropertyDerivative(const unsigned int enum_ID,
                                                      const unsigned int param_ID,
                                                      double& value) const
{
  // check if the data has been initialized
  Assert (this->local_params_initialized == true,
          ExcInvalidState());
  
  Assert (this->global_params_initialized == true,
          ExcInvalidState());
  
  // make sure that the property name asked for, exists, and if it
  // does, get the value of derivative for same
  Property::PropertyCard<double>::PropertyMap::const_iterator property_it = 
    this->property_reference_values.find(enum_ID);
  
  AssertThrow(property_it != this->property_reference_values.end(), 
         Property::PropertyCardBase::ExcPropertyDoesNotExist
         (this->getPropertyEnumName(enum_ID)));
  
  value = 0.0;
  
  const unsigned int param_internal_ID = this->getParameterInternalID(param_ID),
    prop_int_ID = this->getPropertyInternalID(enum_ID);
  
  double func_value = 0.0;
  
  func_value = (*this->parameter_function_value_table)(prop_int_ID, param_internal_ID, 1);
  value += func_value;
  
}



inline
void
Property::PropertyCard<double>::
getPropertyDerivativeForGlobalParameter(const unsigned int enum_ID ,
                                        const unsigned int param_ID,
                                        double& value) const
{
  // make sure that this parameter is a global param
  Property::PropertyCardBase::GlobalParameterMap::const_iterator it, end;
  it = this->global_parameter_map.find(param_ID);
  end = this->global_parameter_map.end();
  
  // if the property card is not dependent on this parameter, then the sensitivity
  // should be zero
  if(it == end) 
    {
    value = 0.0;
    return;
    }
  
  unsigned int param_external_ID = it->second->getID();
  
  this->getPropertyDerivative(enum_ID, param_external_ID, value);
}




inline
void
Property::PropertyCard<double>::
getPropertyDerivativeForLocalParameter(const unsigned int enum_ID ,
                                       const unsigned int param_ID,
                                       double& value) const
{
  // make sure that this parameter is a global param
  Property::PropertyCardBase::LocalParameterMap::const_iterator it, end;
  it = this->local_parameter_map.find(param_ID);
  end = this->local_parameter_map.end();
  
  // if the property card is not dependent on this parameter, then the sensitivity
  // should be zero
  if(it == end) 
    {
    value = 0.0;
    return;
    }
  
  unsigned int param_external_ID = it->second->getID();
  
  this->getPropertyDerivative(enum_ID, param_external_ID, value);
}




inline
std::istream& 
Property::PropertyCard<double>::
readFromInputStream(std::istream& input)
{
  std::string tag;
  
  // read the beginning data of the card
  tag = this->getPropertyCardKindEnumName();
  FESystemIOUtility::readFromInput(input, tag);
  
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->cardID);
  
  // now read in the property base values
  {
    unsigned int n_prop = 0, prop_ID = 0;
    double value = 0.0;
    bool insert_success = false;
    FESystemIOUtility::readFromInput(input, "N_PROPERTIES", n_prop);
    
    Property::PropertyCard<double>::PropertyMap::iterator prop_it;
    for (unsigned int i=0; i<n_prop; i++)
      {
      tag.clear(); value = 0.0;
      input >> tag;
      prop_ID = this->getPropertyEnumID(tag);
      input >> value;
      // now insert the value in the map, and also in the external to 
      // internal ID map
      insert_success = this->property_reference_values.insert
        (Property::PropertyCard<double>::PropertyMap::value_type
         (prop_ID, value)).second;

      AssertThrow(insert_success, 
             Property::PropertyCardBase::ExcDuplicatePropertySpecified(tag));
      
      insert_success = this->property_enumID_to_internal_ID_map.insert
        (Property::PropertyCardBase::IDMap::value_type(prop_ID, i)).second;
      
      // if the previous assert was true, this one should also be true.
      AssertThrow(insert_success, 
             Property::PropertyCardBase::ExcDuplicatePropertySpecified(tag));
      }
  }
  
  // now read the parameters
  this->readParametersFromInputStream(input);
  
  // now read in the property-parameter function ID table
  {
    // first initialize the table
    unsigned int n_prop = this->getNProperties(),
    n_param = this->getNParameters();
    
    this->parameter_function_ID_table->reinit(n_prop, n_param);
    this->parameter_function_ID_table->reset_values(FESystemNumbers::InvalidID);
    
    TableIndices<3> indices(n_prop, n_param, 2);
    this->parameter_function_value_table->reinit(indices);
    this->parameter_function_value_table->reset_values(0.0);
    
    unsigned int n_func_IDs, prop_ID = 0, param_ID = 0, func_ID = 0,
      prop_internal_ID = 0, param_internal_ID = 0;
    
    FESystemIOUtility::readFromInput(input, "N_FUNCTION_IDS", n_func_IDs);
    
    for (unsigned int i=0; i<n_func_IDs; i++)
      {
      tag.clear(); 
      input >> tag;
      prop_ID = this->getPropertyEnumID(tag);
      input >> param_ID;
      input >> func_ID;
      
      // get the internal IDS
      prop_internal_ID = this->getPropertyInternalID(prop_ID);
      param_internal_ID = this->getParameterInternalID(param_ID);
      
      // and, insert the value in the table
      // first make sure that this property hasn't already been set. 
      AssertThrow((*this->parameter_function_ID_table)(prop_internal_ID, param_internal_ID)  ==
                  FESystemNumbers::InvalidID,
                  FESystemExceptions::ExcDuplicateID("Property Function",prop_internal_ID));
      
      (*this->parameter_function_ID_table)(prop_internal_ID, param_internal_ID) = func_ID;
      }
  }
  
  
  // the card should end with the END tag
  tag = this->getPropertyCardKindEnumName();
  FESystemIOUtility::readFromInput(input, tag);

  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}

#endif //__fesystem_property_card_h__

