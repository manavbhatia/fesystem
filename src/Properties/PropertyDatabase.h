// $Id: PropertyDatabase.h,v 1.4 2006-09-05 20:41:51 manav Exp $

#ifndef __fesystem_property_database_h__
#define __fesystem_property_database_h__

// C++ includes
#include <iostream>
#include <map>


// FESystem includes
#include "FESystem/FESystemExceptions.h"


/// PropertyList class performs the function of reading the property cards
/// from a specified input stream, storing them, and returning the card when
/// asked for by the user

namespace FESystem
{
  class FESystemController;
}

class ElemDataCard;

namespace Property
{
  // Forward decleratoins
  class PropertyCardBase;
  template<typename T> class PropertyCard; 
    
  class PropertyDatabase
  {
public:
    /// default constructor
    PropertyDatabase(FESystem::FESystemController& controller);
    
    /// destructor
    ~PropertyDatabase();
    
    /// returns the material property card for this ID
    template<typename T>
    Property::PropertyCard<T>& getMaterialPropertyCardFromID(const unsigned int card_id);  
    
    /// returns the element property card for this ID
    ElemDataCard& getElemDataCardFromID(const unsigned int card_id);
    
    
    /// this method initializes all element data and property cards at the local parameters 
    /// provided in the input. 
    void initAllCardsForLocalParameters(std::map<unsigned int, double>* value_map = NULL);
    
    /// this method initializes all element data cards at the global parameters provided in the
    /// input. 
    void initAllCardsForGlobalParameters(std::map<unsigned int, double>& value_map);

    
    /// clears the local parameter initialization of all property card and elem data cards
    void clearLocalParameterInitializationForAllCards();

    /// clears the global parameter initialization of all property card and elem data cards
    void clearGlobalParameterInitializationForAllCards();
    
    /// method reads the data from an input stream
    std::istream& readFromInputStream( std::istream& input);
    
    /// overloaded input stream operator
    friend inline std::istream& operator>> 
      (std::istream& input, Property::PropertyDatabase& database);
    
protected:	
      
    /// adds a material property card to the map  
    void readAndStoreMaterialCard(std::istream& input,
                                  const unsigned int card_type_enum_ID);
    
    /// adds an element property card to the map
    void readAndStoreElemDataCard(std::istream& input, 
                                  const unsigned int card_type_enum_ID);  
        
    /// local type definitions
    typedef std::map<unsigned int, Property::PropertyCardBase*> PropertyCardMap;
    typedef std::map<unsigned int, ElemDataCard*> ElemDataCardMap;
    
    /// FESystem controller
    FESystem::FESystemController& fesystem_controller;
    
    /// this map stores the material property cards
    PropertyCardMap ID_to_material_property_card_map;	
    
    /// this map stores the elem property cards
    ElemDataCardMap ID_to_elem_data_card_map; 
  };
  
  
  inline std::istream& operator>> 
    (std::istream& input, Property::PropertyDatabase& database)
    {
      database.readFromInputStream(input);
      return input;
    }
  
  
}




template<typename T>
Property::PropertyCard<T>& 
Property::PropertyDatabase::getMaterialPropertyCardFromID(const unsigned int ID)
{
	// get the iterator for this ID, and then make sure that the card exists
	Property::PropertyDatabase::PropertyCardMap::const_iterator 
  card_it, card_end;
  card_it =	this->ID_to_material_property_card_map.find(ID);
	card_end =	this->ID_to_material_property_card_map.end();
  
	Assert(card_it != card_end,
         FESystemExceptions::ExcIDDoesNotExist("Material Property Card",ID));
  
  // next, dynamic cast the map, and return it
  Property::PropertyCard<T> *card = 
    dynamic_cast<Property::PropertyCard<T>*> (card_it->second);
  
  return *card;
}


#endif // __fesystem_property_database_h__

