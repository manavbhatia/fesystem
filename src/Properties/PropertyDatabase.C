// $Id: PropertyDatabase.C,v 1.6 2007-01-05 03:08:50 manav Exp $

// C++ inlcudes
#include <string>


// FESystem inlcudes
#include "Properties/PropertyDatabase.h"
#include "FESystem/FESystemController.h"
#include "Utilities/InputOutputUtility.h"
#include "Properties/MaterialPropertyNameEnums.h"
#include "Properties/IsotropicMaterialPropertyCard.h"
#include "Properties/Isotropic1DElemDataCard.h"
#include "Properties/Isotropic2DElemDataCard.h"
#include "Properties/Isotropic3DElemDataCard.h"
#include "Properties/ForcedConvection1DElemDataCard.h"

Property::PropertyDatabase::PropertyDatabase(FESystem::FESystemController& controller):
fesystem_controller(controller)
{
	
}







Property::PropertyDatabase::~PropertyDatabase()
{
  {
    Property::PropertyDatabase::PropertyCardMap::iterator it, end;
    it = this->ID_to_material_property_card_map.begin();
    end = this->ID_to_material_property_card_map.end();
    
    for (; it!= end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
  }
  
  {
    Property::PropertyDatabase::ElemDataCardMap::iterator it, end;
    it = this->ID_to_elem_data_card_map.begin();
    end = this->ID_to_elem_data_card_map.end();
    
    for (; it!= end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
  }
}










std::istream& 
Property::PropertyDatabase::readFromInputStream( std::istream& input)
{
	std::string tag;
	unsigned int n_cards=0;
	
  FESystemIOUtility::readFromInput(input, "PROPERTY_DATABASE");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  {
    FESystemIOUtility::readFromInput(input, "MATERIAL_PROPERTY_CARDS");
    FESystemIOUtility::readFromInput(input, "BEGIN");
    
    FESystemIOUtility::readFromInput(input, "N_MATERIAL_PROPERTY_CARDS", n_cards);
    
    unsigned int card_enum_ID = 0;
    
    for (unsigned int i=0; i < n_cards; i++)
      {
      FESystemIOUtility::peekFromInput(input, tag);
      card_enum_ID = MaterialPropertyCardEnum::enumID(tag);
      
      this->readAndStoreMaterialCard(input, card_enum_ID);
      }
    
    FESystemIOUtility::readFromInput(input, "MATERIAL_PROPERTY_CARDS");
    FESystemIOUtility::readFromInput(input, "END");
  }
  
  {
    FESystemIOUtility::readFromInput(input, "ELEM_DATA_CARDS");
    FESystemIOUtility::readFromInput(input, "BEGIN");
    
    FESystemIOUtility::readFromInput(input, "N_ELEM_DATA_CARDS", n_cards);
    
    unsigned int card_enum_ID = 0;
    
    for (unsigned int i=0; i < n_cards; i++)
      {
      FESystemIOUtility::peekFromInput(input, tag);
      card_enum_ID = ElemDataCardEnum::enumID(tag);
      
      this->readAndStoreElemDataCard(input, card_enum_ID);
      }
    
    FESystemIOUtility::readFromInput(input, "ELEM_DATA_CARDS");
    FESystemIOUtility::readFromInput(input, "END");
  }
	
	
  FESystemIOUtility::readFromInput(input, "PROPERTY_DATABASE");
  FESystemIOUtility::readFromInput(input, "END");
	
	return input;
}





ElemDataCard& 
Property::PropertyDatabase::getElemDataCardFromID
(const unsigned int ID)
{
	// get the iterator for this ID, and then make sure that the card exists
	Property::PropertyDatabase::ElemDataCardMap::const_iterator 
  card_it, card_end;
  card_it =	this->ID_to_elem_data_card_map.find(ID);
	card_end =	this->ID_to_elem_data_card_map.end();
  
	Assert(card_it != card_end,
         FESystemExceptions::ExcIDDoesNotExist("ElemDataCard", ID));
	
	return *(card_it->second);
}



void
Property::PropertyDatabase::readAndStoreMaterialCard
(std::istream& input,
 const unsigned int card_type_enum_ID)
{
  Property::PropertyCardBase *card = NULL;
  
  // first read in the card
  switch (card_type_enum_ID)
    {
    case ISOTROPIC_MATERIAL_PROPERTY_CARD_ENUM_ID:
      {
        IsotropicMaterialPropertyCard *mat_card = new IsotropicMaterialPropertyCard;
        input >> (*mat_card);
        card = mat_card;
      }
      break;
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (MaterialPropertyCardEnum::enumName(card_type_enum_ID)));
    }
  
  // attach the function database to the card
  card->attachFunctionDatabase(*(this->fesystem_controller.function_database.get()));
  
  // now add it
  unsigned int ID = card->getID();
  
  // first make sure that the property ID does not already exist in the map
  bool insert_success = 
    this->ID_to_material_property_card_map.insert
    (PropertyDatabase::PropertyCardMap::value_type(ID, card)).second;
  
	if (!insert_success)
    delete card;
    
  Assert(insert_success,
         FESystemExceptions::ExcDuplicateID("MaterialPropertyCard", ID));	
}




void
Property::PropertyDatabase::readAndStoreElemDataCard
(std::istream& input,
 const unsigned int card_type_enum_ID)
{
  ElemDataCard *card = NULL;
  
  // first read in the card
  switch (card_type_enum_ID)
    {
    case ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_ID:
      {
        Isotropic1D_ElemDataCard *data_card = new Isotropic1D_ElemDataCard;
        input >> (*data_card);
        card = data_card;
      }
      break;

    case ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_ID:
      {
        Isotropic2D_ElemDataCard *data_card = new Isotropic2D_ElemDataCard;
        input >> (*data_card);
        card = data_card;
      }
      break;

    case ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_ID:
      {
        Isotropic3D_ElemDataCard *data_card = new Isotropic3D_ElemDataCard;
        input >> (*data_card);
        card = data_card;
      }
      break;

    case FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_ID:
      {
        ForcedConvection_1D_ElemDataCard *data_card = new ForcedConvection_1D_ElemDataCard;
        input >> (*data_card);
        card = data_card;
      }
      break;

    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardEnum::enumName(card_type_enum_ID)));
    }
  
  
  // attach the function database to the card
  card->attachFunctionDatabase(*(this->fesystem_controller.function_database.get()));
  
  // attach the property database to this card
  card->attachPropertyDatabase(*this);
  
  // now add it
  unsigned int ID = card->getID();
  
  // first make sure that the property ID does not already exist in the map
  bool insert_success = 
    this->ID_to_elem_data_card_map.insert
    (PropertyDatabase::ElemDataCardMap::value_type(ID, card)).second;
  
	if (!insert_success)
    delete card;
  
  Assert(insert_success,
         FESystemExceptions::ExcDuplicateID("ElemDataCard" ,ID));	
}



void
Property::PropertyDatabase::initAllCardsForLocalParameters
(std::map<unsigned int, double>* value_map)
{
  // get iterators for all elem data cards, and initialize them
  {
    Property::PropertyDatabase::ElemDataCardMap::iterator it, end;
    it = this->ID_to_elem_data_card_map.begin();
    end = this->ID_to_elem_data_card_map.end();
    
    for (; it != end; it++)
      it->second->reinitLocalParameters(value_map);
  }
  
  // get iterators for all material property cards, and initialize them
  {
    Property::PropertyDatabase::PropertyCardMap::iterator it, end;
    it = this->ID_to_material_property_card_map.begin();
    end = this->ID_to_material_property_card_map.end();
    
    for (; it != end; it++)
      it->second->reinitLocalParameters(value_map);
  }
  
}





void
Property::PropertyDatabase::initAllCardsForGlobalParameters
(std::map<unsigned int, double>& value_map)
{
  // get iterators for all elem data cards, and initialize them
  {
    Property::PropertyDatabase::ElemDataCardMap::iterator it, end;
    it = this->ID_to_elem_data_card_map.begin();
    end = this->ID_to_elem_data_card_map.end();
    
    for (; it != end; it++)
      it->second->reinitGlobalParameters(value_map);
  }
  
  // get iterators for all material property cards, and initialize them
  {
    Property::PropertyDatabase::PropertyCardMap::iterator it, end;
    it = this->ID_to_material_property_card_map.begin();
    end = this->ID_to_material_property_card_map.end();
    
    for (; it != end; it++)
      it->second->reinitGlobalParameters(value_map);
  }
  
}



void
Property::PropertyDatabase::clearLocalParameterInitializationForAllCards()
{
  // get iterators for all elem data cards, and initialize them
  {
    Property::PropertyDatabase::ElemDataCardMap::iterator it, end;
    it = this->ID_to_elem_data_card_map.begin();
    end = this->ID_to_elem_data_card_map.end();
    
    for (; it != end; it++)
      it->second->clearLocalParameterInitialization();
  }
  
  // get iterators for all material property cards, and initialize them
  {
    Property::PropertyDatabase::PropertyCardMap::iterator it, end;
    it = this->ID_to_material_property_card_map.begin();
    end = this->ID_to_material_property_card_map.end();
    
    for (; it != end; it++)
      it->second->clearLocalParameterInitialization();
  }
  
}




void
Property::PropertyDatabase::clearGlobalParameterInitializationForAllCards()
{
  // get iterators for all elem data cards, and initialize them
  {
    Property::PropertyDatabase::ElemDataCardMap::iterator it, end;
    it = this->ID_to_elem_data_card_map.begin();
    end = this->ID_to_elem_data_card_map.end();
    
    for (; it != end; it++)
      it->second->clearGlobalParameterInitialization();
  }
  
  // get iterators for all material property cards, and initialize them
  {
    Property::PropertyDatabase::PropertyCardMap::iterator it, end;
    it = this->ID_to_material_property_card_map.begin();
    end = this->ID_to_material_property_card_map.end();
    
    for (; it != end; it++)
      it->second->clearGlobalParameterInitialization();
  }
  
}





