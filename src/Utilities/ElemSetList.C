// $Id: ElemSetList.C,v 1.4 2006-10-29 05:09:18 manav Exp $

// C++ inlcudes
#include <string>
#include <cassert>


// FESystem inlcudes
#include "Utilities/ElemSetList.h"
#include "Utilities/InputOutputUtility.h"


ElemSetList::ElemSetList()
{
	
}







ElemSetList::~ElemSetList()
{
  this->clear();
}







void ElemSetList::clear()
{
  std::map<unsigned int, ElemSet*>::iterator it, end;
  it = this->ID_to_elem_set_map.begin();
  end = this->ID_to_elem_set_map.end();

  for (; it != end; it++)
    {
      delete it->second;
      it->second = NULL;
    }

  this->ID_to_elem_set_map.clear();
}






std::istream& ElemSetList::readFromInputStream( std::istream& input)
{
  unsigned int n_cards=0;
	
  FESystemIOUtility::readFromInput(input, "ELEM_SET_LIST");
  FESystemIOUtility::readFromInput(input, "BEGIN");
		
  FESystemIOUtility::readFromInput(input, "N_ELEM_SETS", n_cards);
	
  for (unsigned int card_incr=0; card_incr < n_cards; card_incr++)
    {
      ElemSet* elem_set = new ElemSet;;
      elem_set->readFromInputStream(input);
		
      this->addElemSetToMap(elem_set->ID(),elem_set);
    }
	
	
  FESystemIOUtility::readFromInput(input, "ELEM_SET_LIST");
  FESystemIOUtility::readFromInput(input, "END");
	
  return input;
}






const ElemSet& ElemSetList::getElemSetFromID( const unsigned int ID) const
{
	// get the iterator for this ID, and then make sure that the card exists
	std::map<unsigned int, ElemSet*>::const_iterator card_it = 
	this->ID_to_elem_set_map.find(ID);
	
	assert(card_it != this->ID_to_elem_set_map.end());
	
	return *(card_it->second);
}





void ElemSetList::addElemSetToMap(const unsigned int ID,  ElemSet* card)
{
	// first make sure that the property ID does not already exist in the map
	assert (this->ID_to_elem_set_map.find(ID) ==
			this->ID_to_elem_set_map.end());	
	
	// then insert the ID, ElemSet pair
	std::pair<std::map<unsigned int, ElemSet*>::iterator, bool> 
	insert_return_pair = 
	this->ID_to_elem_set_map.insert(std::map<unsigned int, ElemSet*>::value_type(ID,card));
	
	// and check the return value for a successful insertion
	assert (insert_return_pair.second == true);
}
