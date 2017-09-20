// $Id: MaterialPropertyNameEnums.h,v 1.6.6.1 2007-05-08 05:19:06 manav Exp $

#ifndef __fesystem_material_property_name_h__
#define __fesystem_material_property_name_h__

// C/C++ includes
#include <string>
#include <map>


// FESystem includes
#include "Utilities/NameEnumHandler.h"

class MaterialPropertyNameEnumerationHandler: public NameEnumerationHandler
{
public: 
  inline void registerEnumData(const unsigned int ID, 
                               const std::string& name,
                               const bool if_isotropic);
  
  
  inline bool ifOnlyIsotropic(const unsigned int ID) const;
protected:
    
    /// map for storing extra information about the properties, whether 
    /// they are isotropic in nature or not
    typedef std::map<unsigned int, bool> PropertyKindMapType;
  static PropertyKindMapType property_kind_map;
  
};

inline
void
MaterialPropertyNameEnumerationHandler::registerEnumData(const unsigned int ID, 
                                                         const std::string& name,
                                                         const bool if_isotropic)
{
  this->registerNameAndID(ID, name);
  
  // find the name in the map. If it does not exist, add it, else return the 
  // id for this name
  MaterialPropertyNameEnumerationHandler::PropertyKindMapType::const_iterator it, end;
  it = this->property_kind_map.find(ID);
  end = this->property_kind_map.end();
  
  // make sure that the id and name are unique
  Assert(it == end, ExcFESystemNameAlreadyExists(ID, name));
  
  // now add the ID to the map
  std::pair<MaterialPropertyNameEnumerationHandler::PropertyKindMapType::iterator, bool> 
    return_pair = this->property_kind_map.insert
    (MaterialPropertyNameEnumerationHandler::PropertyKindMapType::value_type(ID, if_isotropic));
  
  assert (return_pair.second == true); 
}


inline
bool
MaterialPropertyNameEnumerationHandler::ifOnlyIsotropic(const unsigned int ID) const
{
  // get the iterator from the map and return the ID
  MaterialPropertyNameEnumerationHandler::PropertyKindMapType::const_iterator it, end;
  it = this->property_kind_map.find(ID);
  end = this->property_kind_map.end();
  
  Assert(it != end, ExcFESystemIDDoesNotExist(ID));
  
  return it->second;
}



/// this class provides a way to define property names anywhere in the code, and still
/// be used by the PropertyCard to store values in a map. This is an abstract class,
/// hence the property names will be derived from this class. 
class MaterialPropertyNameBase
{
public: 
  /// constructor. 
  /// @param ID  the enumeration ID of this class
  /// @param name the enumeration name of this class
  inline MaterialPropertyNameBase(unsigned int ID,
                                  const std::string& name,
                                  bool only_isotropic);
  
  /// this provides a method to obtain the string name of the enumeration
  /// with the enum id is provided
  /// @param enum_id integer ID of the enumeration
  static inline const std::string& enumName(const unsigned int enum_id);
  
  /// this provides a method to obtain the enum id of the enumeration
  /// when the string name is provided
  static inline unsigned int enumID(const std::string& name);
  
  /// returns whether the enumeration is only an isotropic quantity
  static inline bool ifOnlyIsotropic(const unsigned int ID);
  
  /// returns whether the enumeration is only an isotropic quantity
  static inline bool ifOnlyIsotropic(const std::string& name);
  
  
  /// @returns the enum id of this object
  inline unsigned int enumID() const;
  
  /// @returns the enum name of this object
  inline const std::string& enumName() const;
  
  /// @returns if the quantity is isotropic or not
  inline bool ifOnlyIsotropic() const;
  
protected:
    
    /// enumeration ID of this object
    const  unsigned int enum_id;
  
  /// enumeration name of this object
  const std::string enum_name;
  
  /// whether this property is a pure isotropic quantity, in any dimensional space
  const bool if_only_isotropic;
  
  /// this stores the enumeration constants for all derived classes
  static MaterialPropertyNameEnumerationHandler enum_handler;
  
  
  /// this defines an object that will serve the purpose of registering 
  /// the object to the object handler. It is important that a static object
  /// of this class be added inside every derived enum class, since it will register the
  /// enum class with this base.
  template <typename EnumClass>
    struct EnumRegistration
    { EnumRegistration()
      {
      MaterialPropertyNameBase::enum_handler.registerEnumData(EnumClass::num(),
                                                              EnumClass::name(),
                                                              EnumClass::only_isotropic);
      }
    };
};


inline 
MaterialPropertyNameBase::MaterialPropertyNameBase(const unsigned int ID, 
                                                   const std::string& name,
                                                   const bool isotropic):
enum_id(ID),
enum_name(name),
if_only_isotropic(isotropic)
{}





inline 
const std::string& 
MaterialPropertyNameBase::enumName(const unsigned int enum_id)
{return MaterialPropertyNameBase::enum_handler.getEnumerationName(enum_id);}





inline
unsigned int 
MaterialPropertyNameBase::enumID(const std::string& name)
{return MaterialPropertyNameBase::enum_handler.getEnumerationID(name);}




inline
bool 
MaterialPropertyNameBase::ifOnlyIsotropic(const unsigned int ID) 
{return MaterialPropertyNameBase::enum_handler.ifOnlyIsotropic(ID);}



inline
bool 
MaterialPropertyNameBase::ifOnlyIsotropic(const std::string& name)
{
  unsigned int ID = MaterialPropertyNameBase::enumID(name);
  return MaterialPropertyNameBase::enum_handler.ifOnlyIsotropic(ID);
}





inline unsigned int 
MaterialPropertyNameBase::enumID() const
{return this->enum_id;}




inline const std::string& 
MaterialPropertyNameBase::enumName() const
{return this->enum_name;}





inline bool
MaterialPropertyNameBase::ifOnlyIsotropic() const
{return this->if_only_isotropic;}

//***************************************************************************************//



//***************************************************************************************//

/// define a macro that will instantiate the enumeration name class and 
/// also register it
#define DeclareMaterialEnumName(ClassName, id, string_name,    \
                                if_isotropic)              \
class ClassName: public MaterialPropertyNameBase           \
{                                                          \
public:                                                   \
  ClassName():MaterialPropertyNameBase(id , string_name, if_isotropic)\
  {}                                                         \
  \
  static std::string name()                                \
  {return string_name;}                                   \
  \
  static unsigned int num()                                \
  {return id;}      \
  \
  static const bool only_isotropic = if_isotropic;             \
    \
protected:                                                 \
    static MaterialPropertyNameBase::EnumRegistration<ClassName> registration; \
}


#ifndef DENSITY_ENUM_ID
#define DENSITY_ENUM_ID  1
#else
#error
#endif

#ifndef DENSITY_ENUM_NAME
#define DENSITY_ENUM_NAME "DENSITY"
#else
#error
#endif


#ifndef SPECIFIC_HEAT_ENUM_ID
#define SPECIFIC_HEAT_ENUM_ID 2
#else
#error
#endif

#ifndef SPECIFIC_HEAT_ENUM_NAME
#define SPECIFIC_HEAT_ENUM_NAME  "SPECIFIC_HEAT"
#else
#error
#endif


#ifndef ALPHA_EXPANSION_ENUM_ID
#define ALPHA_EXPANSION_ENUM_ID 3
#else
#error
#endif

#ifndef ALPHA_EXPANSION_ENUM_NAME
#define ALPHA_EXPANSION_ENUM_NAME "ALPHA_EXPANSION"
#else
#error
#endif


#ifndef THERMAL_CONDUCTIVITY_ENUM_ID
#define THERMAL_CONDUCTIVITY_ENUM_ID 4
#else
#error
#endif

#ifndef THERMAL_CONDUCTIVITY_ENUM_NAME
#define THERMAL_CONDUCTIVITY_ENUM_NAME  "THERMAL_CONDUCTIVITY"
#else
#error
#endif


#ifndef YOUNGS_MODULUS_ENUM_ID
#define YOUNGS_MODULUS_ENUM_ID 5
#else
#error
#endif


#ifndef YOUNGS_MODULUS_ENUM_NAME
#define YOUNGS_MODULUS_ENUM_NAME "YOUNGS_MODULUS"
#else
#error
#endif


#ifndef POISSONS_RATIO_ENUM_ID
#define POISSONS_RATIO_ENUM_ID 6
#else
#error
#endif

#ifndef POISSONS_RATIO_ENUM_NAME
#define POISSONS_RATIO_ENUM_NAME "POISSONS_RATIO"
#else
#error
#endif


#ifndef SPRING_CONSTANT_ENUM_ID
#define SPRING_CONSTANT_ENUM_ID 7
#else
#error
#endif

#ifndef SPRING_CONSTANT_ENUM_NAME
#define SPRING_CONSTANT_ENUM_NAME "SPRING_CONSTANT"
#else
#error
#endif


#ifndef EMISSIVITY_ENUM_ID
#define EMISSIVITY_ENUM_ID 8
#else
#error
#endif

#ifndef EMISSIVITY_ENUM_NAME
#define EMISSIVITY_ENUM_NAME "EMISSIVITY"
#else
#error
#endif



DeclareMaterialEnumName(DENSITY, DENSITY_ENUM_ID, DENSITY_ENUM_NAME, true);
DeclareMaterialEnumName(SPECIFIC_HEAT, SPECIFIC_HEAT_ENUM_ID , SPECIFIC_HEAT_ENUM_NAME, true);
DeclareMaterialEnumName(ALPHA_EXPANSION, ALPHA_EXPANSION_ENUM_ID, ALPHA_EXPANSION_ENUM_NAME, false);
DeclareMaterialEnumName(THERMAL_CONDUCTIVITY, THERMAL_CONDUCTIVITY_ENUM_ID, THERMAL_CONDUCTIVITY_ENUM_NAME, false);
DeclareMaterialEnumName(YOUNGS_MODULUS, YOUNGS_MODULUS_ENUM_ID, YOUNGS_MODULUS_ENUM_NAME, false);
DeclareMaterialEnumName(POISSONS_RATIO, POISSONS_RATIO_ENUM_ID, POISSONS_RATIO_ENUM_NAME, false);
DeclareMaterialEnumName(SPRING_CONSTANT, SPRING_CONSTANT_ENUM_ID, SPRING_CONSTANT_ENUM_NAME, true);
DeclareMaterialEnumName(EMISSIVITY, EMISSIVITY_ENUM_ID, EMISSIVITY_ENUM_NAME, true);


// ************************ define the base class for material cards *******************
// define a enumeration class for material property cards
DeclareEnumClass(MaterialPropertyCardEnum);


#endif // __fesystem_material_property_name_h__
