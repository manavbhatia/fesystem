// $Id: OrthotropicMaterialPropertyCard.h,v 1.2 2006-09-05 20:41:51 manav Exp $

#ifndef __fesystem_orthotropic_material_property_card_h__
#define __fesystem_orthotropic_material_property_card_h__

// C++ includes 
#include <iostream>
#include <map>

// FESystem includes
#include "Properties/MaterialPropertyNameEnums.h"
#include "Properties/PropertyCard.h"

// Forward declerations
class TensorBase;
class MaterialPropertyNameBase;



#define ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_ID 1
#define ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_NAME "ORTHOTROPIC_MATERIAL_PROPERTY_CARD"
DeclareEnumName(ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM, MaterialPropertyCardType, 
                ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_ID, 
                ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_NAME);


/// The \p PropertyCard class that provides the necessary data structure 
/// for storage of material and structural properties. This is a base class, 
/// and different material cards can derive from this. Each property card
/// can have dependence on some parameter values, which is defined through the input. If a 
/// parameter dependence has been defined, then before the parameter values can be accessed, 
/// the card needs to be initialized at some parameter values. During reinitialization, if
/// no value has been specified for a parameter, a default value for the same is assumed.
/// Between any two consecutive initializations, the card needs to be cleared of any previous
/// initialization. This is needed to make sure that the card is cleared of any remnant value
/// initializations to avoid errors. Once initialized, the value of the property can be accessed from
/// the card.
//
template <unsigned int dim>
class OrthotropicMaterialPropertyCard : 
public Property::PropertyCard<TensorBase*, MaterialPropertyNameBase>
{
public:	
  
  /// default constructor
  OrthotropicMaterialPropertyCard();
  
  /// destructor
  virtual ~OrthotropicMaterialPropertyCard();
  
  /// @returns the property card kind enumeration ID
  virtual const unsigned int getPropertyCardKindEnumID() const;
  
  /// @returns the property card kind enumeration name
  virtual const std::string getPropertyCardKindEnumName() const;
  
  
  /// @returns the dimensionality of properties 
  virtual const unsigned int getDimension() const;
  
  /// overloaded input operator 
  template <unsigned int dim_>
    friend std::istream& operator>> (std::istream& input, OrthotropicMaterialPropertyCard<dim_>& card);
  
  /// overloaded output stream operator
  template <unsigned int dim_>
    friend std::ostream& operator<< (std::ostream& os , const OrthotropicMaterialPropertyCard<dim_>& card);
  
protected:
    
    /// this method reads the information for this card from the input stream provided. This
    /// has been templetized based on a dimension parameter, so that it can be used for isotropic 
    /// as well as orthotropic cards
    std::istream& readFromInputStream(std::istream& input);
};



template <unsigned int dim>
const unsigned int 
OrthotropicMaterialPropertyCard<dim>::getPropertyCardKindEnumID() const
{
  return ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM::num();
}



template <unsigned int dim>
const std::string 
OrthotropicMaterialPropertyCard<dim>::getPropertyCardKindEnumName() const
{
  return ORTHOTROPIC_MATERIAL_PROPERTY_CARD_ENUM::name();
}



template <unsigned int dim>
const unsigned int 
OrthotropicMaterialPropertyCard<dim>::getDimension() const
{
  return dim;
}





template <unsigned int dim>
std::istream&
OrthotropicMaterialPropertyCard<dim>::readFromInputStream(std::istream& input)
{
  //   std::string tag;
  //   unsigned int enum_id = 0;
  //   double value = 0.0;
  //   MaterialPropertyCard::PropertyName name = MaterialPropertyCard::INVALID_PROPERTY_NAME;
  
  //   tag.clear();
  //   is >> tag;
  //   assert (tag == "BEGIN_PROPERTY_CARD");
	
  //   tag.clear();
  //   is >> tag;
  //   assert (tag == "ID");
  //   is >> card.cardID;
  
  //   tag.clear();
  //   is >> tag;
  //   while (tag != "END_PROPERTY_CARD")
  //     {
  //       enum_id = EnumBase::getEnumerationID(tag);
		
  //       //make sure that the property does not already exist
  //       assert (card.property_values.find(enum_id) == card._property_values.end());
		
  //       is >> value;
  //       std::pair<MaterialPropertyCard::PropertyMapType::iterator, bool> 
  // 	insert_return_value = card._property_values.insert
  // 	(MaterialPropertyCard::PropertyMapType::value_type(enum_id,value));
		
  //       assert (insert_return_value.second == true);
		
  //       tag.clear(); value = 0.0; enum_id = 0;
  //       is >> tag;
  //     }
	
  return input;
}

typedef OrthotropicMaterialPropertyCard<2> OrthotropicMaterialPropertyCard2D;
typedef OrthotropicMaterialPropertyCard<3> OrthotropicMaterialPropertyCard3D;

#endif // #define __fesystem_orthotropic_material_property_card_h__
