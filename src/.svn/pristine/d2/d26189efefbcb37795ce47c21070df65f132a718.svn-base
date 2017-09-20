// $Id: ElemSetList.h,v 1.2 2006-09-05 20:41:44 manav Exp $

#ifndef __elem_set_list_h__
#define __elem_set_list_h__

// C++ includes
#include <map>

// FESystem inlcudes
#include "ElemSet.h"

/// Elem_SetList class performs the function of reading the elem_set cards
/// from a specified input stream, storing them, and returning the card when
/// asked for by the user

class ElemSetList
{
public:

	ElemSetList();
	
	~ElemSetList();
	
	void clear();
	
	std::istream& readFromInputStream( std::istream& );
	
	/// returns the elem_set card for this ID
	const ElemSet& getElemSetFromID( const unsigned int ) const;
	
protected:	
	
	void addElemSetToMap(const unsigned int , ElemSet* );
	
	std::map<unsigned int, ElemSet*> ID_to_elem_set_map;		
};


#endif // __elem_set_list_h__
