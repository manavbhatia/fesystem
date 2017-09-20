// $Id: IsotropicElemDataCard.h,v 1.6.6.1 2007-03-14 22:05:03 manav Exp $

#ifndef __fesystem_isotropic_elem_data_card_h__
#define __fesystem_isotropic_elem_data_card_h__

// FESystem includes
#include "Properties/ElemDataCard.h"
#include "Properties/PropertyDatabase.h"


class IsotropicElemDataCard: public ElemDataCard
{
public:
  inline IsotropicElemDataCard();
  
  virtual inline ~IsotropicElemDataCard();
  
  /// returns the material card ID for this card
  inline unsigned int getMaterialCardID() const;
  
  /// returns the material card for this card
  inline Property::PropertyCard<double>& getMaterialCard();

  
  /// @returns true if either the element data card, or the material property card are
  /// dependent on the specified global parameter
  virtual inline bool checkElemAndMaterialCardGlobalParameterDependence
    (const unsigned int global_param_ID);
  
  /// @returns true if either the element data card, or the material property card are
  /// dependent on the specified local parameter
  virtual inline bool checkElemAndMaterialCardLocalParameterDependence
    (const unsigned int local_param_ID);

  /// initializes the card and the material card for the given local parameters
  virtual inline void reinitElemAndMaterialCardForLocalParameters
    (std::map<unsigned int, double> *value_map = NULL);
  
  /// initializes the card and the material card for the given global parameters
  virtual inline void reinitElemAndMaterialCardForGlobalParameters
    (std::map<unsigned int, double> &value_map);
  
  /// clears initialization of the card and the material card for local parameters
  virtual inline void clearElemAndMaterialCardLocalParameterInitialization();
  
  /// clears initialization of the card and the material card for local parameters
  virtual inline void 
    partialClearElemAndMaterialCardLocalParameterInitialization
    (std::vector<unsigned int >& param_ids);

  /// clears initialization of the card and the material card for global parameters
  virtual inline void clearElemAndMaterialCardGlobalParameterInitialization();
  
  /// returns value from material card
  virtual inline void 
    getPropertyValueFromMaterialCard(const unsigned int prop_enum_ID, double& value,
                                     const unsigned int layer = 0);
  
  /// returns value of the derivative of a property from the material card
  virtual inline void
    getPropertyValueDerivativeForGlobalParameterFromMaterialCard
    (const unsigned int prop_enum_ID, 
     const unsigned int global_param_ID,
     double& value, 
     const unsigned int layer = 0);

  
  /// returns value of the derivative of a property from the material card
  virtual inline void
    getPropertyValueDerivativeForLocalParameterFromMaterialCard
    (const unsigned int prop_enum_ID, 
     const unsigned int local_param_ID,
     double& value,
     const unsigned int layer = 0);
  
  virtual inline std::istream& readFromInputStream(std::istream& input);

protected:
    
    /// material card ID for this card
    unsigned int material_ID;

  /// this is a pointer to the material property card for this element, and is stored for 
    /// convenience
    Property::PropertyCard<double> *material_card;
  

};




inline
IsotropicElemDataCard::IsotropicElemDataCard():
ElemDataCard(),
material_ID(FESystemNumbers::InvalidID),
material_card(NULL)
{
  
}



inline
IsotropicElemDataCard::~IsotropicElemDataCard()
{
  
}




inline 
void
IsotropicElemDataCard::reinitElemAndMaterialCardForLocalParameters
(std::map<unsigned int, double> *value_map)
{
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  this->reinitLocalParameters(value_map);
  property.reinitLocalParameters(value_map);
}






inline
void
IsotropicElemDataCard::reinitElemAndMaterialCardForGlobalParameters
(std::map<unsigned int, double> &value_map)
{
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  this->reinitGlobalParameters(value_map);
  property.reinitGlobalParameters(value_map);
}






inline
void
IsotropicElemDataCard::clearElemAndMaterialCardLocalParameterInitialization()
{
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  this->clearLocalParameterInitialization();
  property.clearLocalParameterInitialization();
}



inline
void
IsotropicElemDataCard::
partialClearElemAndMaterialCardLocalParameterInitialization(std::vector<unsigned int>& param_ids)
{
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  this->partialClearLocalParameterInitialization(param_ids);
  property.partialClearLocalParameterInitialization(param_ids);
}



inline
void
IsotropicElemDataCard::clearElemAndMaterialCardGlobalParameterInitialization()
{
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  this->clearGlobalParameterInitialization();
  property.clearGlobalParameterInitialization();
}




inline unsigned int 
IsotropicElemDataCard::getMaterialCardID() const
{
  Assert(this->material_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(this->material_ID));
  
  return this->material_ID;
}



inline
Property::PropertyCard<double>&
IsotropicElemDataCard::getMaterialCard()
{
  if (this->material_card == NULL)
    this->material_card =
      &(this->property_database->getMaterialPropertyCardFromID<double>(this->material_ID));
  
  
  return *(this->material_card);
}




inline void 
IsotropicElemDataCard::getPropertyValueFromMaterialCard(const unsigned int prop_enum_ID, 
							double& value, 
                                                        const unsigned int layer)
{
  // layer not used for an Isotropic card
  (void) layer;
  
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  property.getPropertyValue(prop_enum_ID, value);
}



inline void
IsotropicElemDataCard::getPropertyValueDerivativeForGlobalParameterFromMaterialCard
(const unsigned int prop_enum_ID, 
 const unsigned int global_param_ID,
 double& value, 
 const unsigned int layer)
{
  // layer not used for an Isotropic card
  (void) layer;
  
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  property.getPropertyDerivativeForGlobalParameter(prop_enum_ID,
                                                   global_param_ID,
                                                   value);
}



inline void
IsotropicElemDataCard::getPropertyValueDerivativeForLocalParameterFromMaterialCard
(const unsigned int prop_enum_ID, 
 const unsigned int local_param_ID,
 double& value, 
 const unsigned int layer)
{
  // layer not used for an Isotropic card
  (void) layer;
  
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  property.getPropertyDerivativeForLocalParameter(prop_enum_ID,
                                                  local_param_ID,
                                                  value);
}


inline
std::istream& 
IsotropicElemDataCard::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  // read the beginning data of the card
  tag = this->getPropertyCardKindEnumName();
  FESystemIOUtility::readFromInput(input, tag);
  
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->cardID);
  
  FESystemIOUtility::readFromInput(input, "MATERIAL_PROPERTY_CARD_ID", this->material_ID);
  
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
      // now insert the value in the map
      insert_success = this->property_reference_values.insert
        (Property::PropertyCard<double>::PropertyMap::value_type
         (prop_ID, value)).second;
      
      Assert(insert_success, 
             Property::PropertyCardBase::ExcDuplicatePropertySpecified(tag));
      
      insert_success = this->property_enumID_to_internal_ID_map.insert
        (Property::PropertyCardBase::IDMap::value_type(prop_ID, i)).second;
      
      // if the previous assert was true, this one should also be true.
      Assert(insert_success, 
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
      Assert((*this->parameter_function_ID_table)(prop_internal_ID, param_internal_ID)  ==
             FESystemNumbers::InvalidID,
             FESystemExceptions::ExcDuplicateID("Property Function",prop_internal_ID));
      
      // and, insert the value in the table
      (*this->parameter_function_ID_table)(prop_internal_ID, param_internal_ID) = func_ID;
      }
  }
  
  // the card should end with the END tag
  tag = this->getPropertyCardKindEnumName();
  FESystemIOUtility::readFromInput(input, tag);
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
  
}



inline
bool
IsotropicElemDataCard::checkElemAndMaterialCardGlobalParameterDependence
(const unsigned int global_param_ID)
{
  if (this->checkGlobalParameterDependence(global_param_ID) &&
      this->getMaterialCard().checkGlobalParameterDependence(global_param_ID))
    return true;
  else
    return false;
}




inline
bool
IsotropicElemDataCard::checkElemAndMaterialCardLocalParameterDependence
(const unsigned int local_param_ID)
{
  if (this->checkLocalParameterDependence(local_param_ID) &&
      this->getMaterialCard().checkLocalParameterDependence(local_param_ID))
    return true;
  else
    return false;
}

#endif // __fesystem_isotropic_elem_data_card_h__
