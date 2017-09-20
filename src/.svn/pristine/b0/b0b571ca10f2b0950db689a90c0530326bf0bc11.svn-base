// $Id: IsotropicMaterialPropertyCard.h,v 1.4 2006-09-05 20:41:51 manav Exp $

#ifndef __fesystem_isotropic_material_property_card_h__
#define __fesystem_isotropic_material_property_card_h__

// C++ includes 
#include <iostream>

// FESystem includes
#include "Properties/PropertyCard.h"
#include "Properties/MaterialPropertyNameEnums.h"


// define the enumeration class for this card type
#ifndef ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_ID
#define ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_ID 1
#else
#error
#endif

#ifndef ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_NAME
#define ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_NAME "ISOTROPIC_MATERIAL_PROPERTY_CARD"
#else 
#error
#endif 

DeclareEnumName(ISOTROPIC_MATERIAL_PROPERTY_CARD, MaterialPropertyCardEnum, 
                ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_ID, 
                ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_NAME);



/// The \p IsotropicMaterialPropertyCard class provides the necessary functionality 
/// for storage of material properties.
//
class IsotropicMaterialPropertyCard : 
public Property::PropertyCard<double>
{
public:	
  
  /// default constructor
  IsotropicMaterialPropertyCard();
  
  /// destructor
  virtual ~IsotropicMaterialPropertyCard();
  
  /// @returns the property card kind enumeration ID
  virtual inline unsigned int getPropertyCardKindEnumID() const;
  
  /// @returns the property card kind enumeration name
  virtual inline const std::string getPropertyCardKindEnumName() const;
  
  /// @returns the property name from enumeration ID of the property
  virtual inline std::string getPropertyEnumName(const unsigned int enum_ID) const;
  
  /// @returns the property ID from enumeration name of the property
  virtual inline const unsigned int getPropertyEnumID(const std::string& enum_name) const;
  
  /// @returns the dimensionality of properties 
  virtual inline unsigned int getDimension() const;
  
  /// overloaded input operator 
  friend std::istream& operator>> (std::istream& input, 
                                   IsotropicMaterialPropertyCard& card);
  
  /// overloaded output stream operator
  friend std::ostream& operator<< (std::ostream& os , 
                                   const IsotropicMaterialPropertyCard& card);

    
protected:

};




inline
unsigned int 
IsotropicMaterialPropertyCard::getPropertyCardKindEnumID() const
{
  return ISOTROPIC_MATERIAL_PROPERTY_CARD::num();
}



inline
const std::string 
IsotropicMaterialPropertyCard::getPropertyCardKindEnumName() const
{
  return ISOTROPIC_MATERIAL_PROPERTY_CARD::name();
}



inline
unsigned int 
IsotropicMaterialPropertyCard::getDimension() const
{
  return 0;
}



inline
std::string
IsotropicMaterialPropertyCard::getPropertyEnumName(const unsigned int enum_ID) const
{
  return MaterialPropertyNameBase::enumName(enum_ID);
}



inline
const unsigned int
IsotropicMaterialPropertyCard::getPropertyEnumID(const std::string& enum_name) const
{
  return MaterialPropertyNameBase::enumID(enum_name); 
}



#endif // #define __fesystem_isotropic_material_property_card_h__
