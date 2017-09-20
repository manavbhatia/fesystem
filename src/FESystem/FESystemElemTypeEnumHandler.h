// $Id: FESystemElemTypeEnumHandler.h,v 1.2 2006-09-05 20:17:36 manav Exp $

#ifndef __fesystem_elem_type_enum_handler_h__
#define __fesystem_elem_type_enum_handler_h__

// C/C++ includes
#include <string>
#include <map>


// FESystem includes
#include "Utilities/NameEnumHandler.h"

// libmesh includes
#include "geom/elem.h"

namespace FESystemElem
{
  
  class FESystemElemNameEnumerationHandler: public NameEnumerationHandler
  {
public: 
    inline void registerEnumData(const unsigned int ID, 
                                 const std::string& name,
                                 const ElemType elem_type);
    
    
    inline ElemType getElemType(const unsigned int ID) const;
protected:
      
      /// map for storing extra information about the elements, about which kind of 
      /// geometric element they depend on
      typedef std::map<unsigned int, ElemType> ElemKindMapType;
    static ElemKindMapType elem_kind_map;
    
  };
  
  inline
    void
    FESystemElemNameEnumerationHandler::registerEnumData(const unsigned int ID, 
                                                         const std::string& name, 
                                                         const ElemType elem_type)
    {
      this->registerNameAndID(ID, name);
      
      // find the name in the map. If it does not exist, add it, else return the 
      // id for this name
      FESystemElemNameEnumerationHandler::ElemKindMapType::const_iterator it, end;
      it = this->elem_kind_map.find(ID);
      end = this->elem_kind_map.end();
      
      // make sure that the id and name are unique
      Assert(it == end, ExcFESystemNameAlreadyExists(ID, name));
      
      // now add the ID to the map
      std::pair<FESystemElem::FESystemElemNameEnumerationHandler::ElemKindMapType::iterator, bool> 
        return_pair = this->elem_kind_map.insert
        (FESystemElemNameEnumerationHandler::ElemKindMapType::value_type(ID, elem_type));
      
      assert (return_pair.second == true); 
    }
  
  
  inline
    ElemType
    FESystemElemNameEnumerationHandler::getElemType(const unsigned int ID) const
    {
      // get the iterator from the map and return the ID
      FESystemElemNameEnumerationHandler::ElemKindMapType::const_iterator it, end;
      it = this->elem_kind_map.find(ID);
      end = this->elem_kind_map.end();
      
      Assert(it != end, ExcFESystemIDDoesNotExist(ID));
      
      return it->second;
    }
  
  /// this class provides a way to define property names anywhere in the code, and still
  /// be used by the PropertyCard to store values in a map. This is an abstract class,
  /// hence the property names will be derived from this class. 
  class FESystemElemTypeEnum
    {
public: 
      /// constructor. 
      /// @param ID  the enumeration ID of this class
      /// @param name the enumeration name of this class
      /// @param type ElemType from libmesh defining which geometric element this is based on
      inline FESystemElemTypeEnum(unsigned int ID,
                                  const std::string& name,
                                  const ElemType type);
      
      /// this provides a method to obtain the string name of the enumeration
      /// with the enum id is provided
      /// @param enum_id integer ID of the enumeration
      static inline const std::string& enumName(const unsigned int enum_id);
      
      /// this provides a method to obtain the enum id of the enumeration
      /// when the string name is provided
      static inline unsigned int enumID(const std::string& name);
      
      /// returns whether the enumeration is only an isotropic quantity
      static inline ElemType elemType(const unsigned int ID);
      
      /// returns whether the enumeration is only an isotropic quantity
      static inline ElemType elemType(const std::string& name);
      
      
      /// @returns the enum id of this object
      inline unsigned int enumID() const;
      
      /// @returns the enum name of this object
      inline const std::string& enumName() const;
      
      /// @returns if the quantity is isotropic or not
      inline ElemType elemType() const;
      
protected:
        
        /// enumeration ID of this object
        const unsigned int enum_id;
      
      /// enumeration name of this object
      const std::string enum_name;
      
      /// whether this property is a pure isotropic quantity, in any dimensional space
      const ElemType elem_type;
      
      /// this stores the enumeration constants for all derived classes
      static FESystemElemNameEnumerationHandler enum_handler;
      
      
      /// this defines an object that will serve the purpose of registering 
      /// the object to the object handler. It is important that a static object
      /// of this class be added inside every derived enum class, since it will register the
      /// enum class with this base.
      template <typename EnumClass>
        struct EnumRegistration
        { EnumRegistration()
          {
          FESystemElem::FESystemElemTypeEnum::enum_handler.registerEnumData
          (EnumClass::num(),
           EnumClass::name(),
           EnumClass::elem_geom_type);
          }
        };
    };
}

inline 
FESystemElem::FESystemElemTypeEnum::FESystemElemTypeEnum(const unsigned int ID, 
                                                         const std::string& name,
                                                         const ElemType type):
enum_id(ID),
enum_name(name),
elem_type(type)
{}





inline 
const std::string& 
FESystemElem::FESystemElemTypeEnum::enumName(const unsigned int enum_id)
{return FESystemElem::FESystemElemTypeEnum::enum_handler.getEnumerationName(enum_id);}





inline
unsigned int 
FESystemElem::FESystemElemTypeEnum::enumID(const std::string& name)
{return FESystemElem::FESystemElemTypeEnum::enum_handler.getEnumerationID(name);}




inline
ElemType 
FESystemElem::FESystemElemTypeEnum::elemType(const unsigned int ID) 
{return FESystemElem::FESystemElemTypeEnum::enum_handler.getElemType(ID);}



inline
ElemType 
FESystemElem::FESystemElemTypeEnum::elemType(const std::string& name)
{
  unsigned int ID = FESystemElem::FESystemElemTypeEnum::enumID(name);
  return FESystemElem::FESystemElemTypeEnum::elemType(ID);
}





inline unsigned int 
FESystemElem::FESystemElemTypeEnum::enumID() const
{return this->enum_id;}




inline const std::string& 
FESystemElem::FESystemElemTypeEnum::enumName() const
{return this->enum_name;}





inline 
ElemType
FESystemElem::FESystemElemTypeEnum::elemType() const
{return this->elem_type;}

//***************************************************************************************//



//***************************************************************************************//

/// define a macro that will instantiate the enumeration name class and 
/// also register it
#define DeclareFESystemElemTypeEnumName(ClassName, id, string_name,    \
                                        type)              \
namespace FESystemElem{ \
  class ClassName: public FESystemElem::FESystemElemTypeEnum           \
  {                                                          \
public:                                                   \
    ClassName():FESystemElem::FESystemElemTypeEnum(id , string_name, type)\
    {}                                                         \
    \
    static std::string name()                                \
    {return string_name;}                                   \
    \
    static unsigned int num()                                \
    {return id;}      \
    \
    static const ElemType elem_geom_type = type;             \
      \
protected:                                                 \
      static FESystemElem::FESystemElemTypeEnum::EnumRegistration<ClassName> registration; \
  }; \
}


#endif // __fesystem_elem_type_enum_handler_h__
