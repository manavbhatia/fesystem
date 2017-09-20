// $Id: NameEnumHandler.h,v 1.6 2006-11-13 00:04:30 manav Exp $

#ifndef __fesystem_name_enum_handler_h__
#define __fesystem_name_enum_handler_h__


// C/C++ includes
#include <iostream>
#include <string>
#include <map>
#include <cassert>


// FESystem includes
#include "FESystem/FESystemExceptions.h"

/// This class provides an interface to define names that can be used as 
/// enumerations. It stores the names and assigns a unique ID to each.
/// Then, each ID can be used as an enumeration constant, with an
/// associated name.
class NameEnumerationHandler
{
public:
  /// this method accepts a class object, whose name is registered in the 
  /// enumeration
  inline void registerNameAndID(const unsigned int ID, 
                                const std::string& name);
  
  /// this method accepts a class object, whose enumeration constant has been 
  /// returned
  inline unsigned int getEnumerationID(const std::string& name) const;
  
  /// this method returns an enumeration name for the given constant
  inline const std::string& getEnumerationName(unsigned int ID) const;
  
protected:
    
    // the following typenames have been defined for local reference
    typedef std::map<std::string, unsigned int> NameToIDMap;
  typedef std::map<unsigned int, std::string> IDToNameMap;
  
  /// this stores the ID for each name
  NameToIDMap name_to_id_map;
  
  /// this stores the name for each ID
  IDToNameMap id_to_name_map;
  
  // decleration of some exceptions
  /// this exception is thrown if the an ID for an enumeration class is already 
  /// taken by another class. If this error occurs then a different ID should be chosen
  DeclException2(ExcFESystemIDAlreadyExists, unsigned int, std::string, 
                 << " Specified ID: " << arg1 << " for Enumeration: " 
                 << arg2 << " already taken by another class.");
  
  
  /// this exception is thrown if the an name for an enumeration class is already 
  /// taken by another class. If this error occurs then a different name should be chosen
  DeclException2(ExcFESystemNameAlreadyExists, unsigned int, std::string, 
                 << " Specified Name: " << arg2 << " for Enumeration ID: " 
                 << arg1 << " already taken by another class.");
  
  /// this expception is thrown if the requested name does not exist in the map
  DeclException1(ExcFESystemNameDoesNotExist, std::string, 
                 << "Requested EnumName: " << arg1 << " does not exist in map.");
  
  /// this expception is thrown if the requested name does not exist in the map
  DeclException1(ExcFESystemIDDoesNotExist, unsigned int, 
                 << "Requested EnumID: " << arg1 << " does not exist in map.");
};


inline void
NameEnumerationHandler::registerNameAndID(const unsigned int ID, 
                                          const std::string& name)
{
  // find the name in the map. If it does not exist, add it, else return the 
  // id for this name
  NameEnumerationHandler::NameToIDMap::const_iterator name_it, name_end;
  name_it = this->name_to_id_map.find(name);
  name_end = this->name_to_id_map.end();
  
  NameEnumerationHandler::IDToNameMap::const_iterator id_it, id_end;
  id_it = this->id_to_name_map.find(ID);
  id_end = this->id_to_name_map.end();
  
  // make sure that the id and name are unique
  Assert(name_it == name_end, ExcFESystemNameAlreadyExists(ID, name));
  Assert(id_it == id_end, ExcFESystemIDAlreadyExists(ID, name));
  
  
  // now add the name and ID to the map
  std::pair<NameEnumerationHandler::NameToIDMap::iterator, bool> 
    return_pair = this->name_to_id_map.insert
    (NameEnumerationHandler::NameToIDMap::value_type(name, ID));
  
  assert (return_pair.second == true); 
  
  std::pair<NameEnumerationHandler::IDToNameMap::iterator, bool>
    reverse_map_return_pair = this->id_to_name_map.insert
    (NameEnumerationHandler::IDToNameMap::value_type(ID, name));
  
  assert (return_pair.second == true);    
}


inline 
unsigned int 
NameEnumerationHandler::getEnumerationID(const std::string& name) const
{
  // get the iterator from the map and return the ID
  NameEnumerationHandler::NameToIDMap::const_iterator it, end;
  it = this->name_to_id_map.find(name);
  end = this->name_to_id_map.end();
  
  AssertThrow(it != end, ExcFESystemNameDoesNotExist(name));
  
  return it->second;
}



inline
const std::string& NameEnumerationHandler::getEnumerationName(unsigned int ID) const
{
  // also add in the reverse id to name map. Get the iterators for this ID
  NameEnumerationHandler::IDToNameMap::const_iterator reverse_map_it, reverse_map_end;
  reverse_map_it = this->id_to_name_map.find(ID);
  reverse_map_end = this->id_to_name_map.end();
  
  // make sure that the value exists
  AssertThrow(reverse_map_it != reverse_map_end, ExcFESystemIDDoesNotExist(ID));
  
  return reverse_map_it->second;
}



/// a macro is defined to create an enumeratoin class. The purpose of this class is to 
/// create a base class for a category of enums, so that all of them can he referred to 
/// in a single category. Every class that is created using this method must also 
/// initialize the static object EnumClass::enum_handler 
#define DeclareEnumClass(EnumClass) \
class EnumClass \
{ \
public:  \
  static const std::string& enumName(const unsigned int enum_id) \
  {return EnumClass::enum_handler.getEnumerationName(enum_id);} \
  static unsigned int enumID(const std::string& name) \
  {return EnumClass::enum_handler.getEnumerationID(name);} \
protected: \
  template <typename EnumNameClass> \
  struct EnumRegistration \
  {EnumRegistration() \
    {EnumClass::enum_handler.registerNameAndID(EnumNameClass::num(), \
                                               EnumNameClass::name());} };\
   static NameEnumerationHandler enum_handler; \
}



/// define a macro that will instantiate an enum name
/// class and also register it
#define DeclareEnumName(ClassName, BaseClass ,id, string_name)  \
class ClassName: public BaseClass                         \
{                                                          \
public:                                                   \
  static const std::string name()                                \
  {return string_name;}                                   \
  static const unsigned int num()                                \
  {return id;}                                              \
  static const bool only_isotropic;                         \
protected:                                                 \
    static BaseClass::EnumRegistration<ClassName> registration; \
}                                                          
//extern BaseClass::EnumRegistration<ClassName> ClassName::registration


#endif // __fesystem_name_enum_handler_h__
