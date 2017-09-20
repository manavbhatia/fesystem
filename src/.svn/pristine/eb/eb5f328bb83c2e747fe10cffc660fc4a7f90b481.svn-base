// $Id: IsotropicMaterialPropertyCard.C,v 1.5.6.1 2007-03-14 22:05:03 manav Exp $

// C++ includes 
#include <iostream>

// FESystem includes
#include "Properties/PropertyCard.h"
#include "Properties/IsotropicMaterialPropertyCard.h"
#include "Utilities/InputOutputUtility.h"


IsotropicMaterialPropertyCard::IsotropicMaterialPropertyCard():
Property::PropertyCard<double>()
{
}



IsotropicMaterialPropertyCard::~IsotropicMaterialPropertyCard()
{
  
}



std::ostream& 
operator<< (std::ostream& os, const IsotropicMaterialPropertyCard& card)
{
  // this needs to be implemented, but for now, this statement will prevent
  // the compiler to give a warning about unused parameter
  (void) card;
  //   TensorBase* value;
  //   unsigned int enum_id;
  
  //   // output the ID of the card
  //   os << "\nMATERIAL_PROPERTY_CARD_DATA_BEGIN" << std::endl
  //      << "\tID " << card.cardID << std::endl
  //      << "\tTYPE " << card.getPropertyCardKindEnumName() << std::endl;
  
  //   // output the parameters for the material card
  //   {
  //     os << "PARAMETERS" << std::endl;
  
  //     MaterialPropertyCard<>::ParameterValueMapType::const_iterator it, end;
  //     it = card.parameter_values.begin();
  //     end = card.parameter_values.end();
  //     for (; it != end; it++)
  //       { os << it->first << "\t" << it->second.ref_value << std::endl;}
  //   }
	
  //   // output the property values at the reference parameters
  //   {
  //     os << "PROPERTY_VALUES_AT_PARAMETER_REFERENCE" << std::endl;
  
  //     MaterialPropertyCard::PropertyMapType::const_iterator p_it, p_end;
  //     p_it = card.property_reference_values.begin();
  //     p_end = card.property_reference_values.end();
	
  //     for (; p_it != p_end; p_it++)
  //       {
  // 	enum_id = p_it->first;
  // 	value = p_it->second;
  
  // 	const std::string& tag = MaterialPropertyNameBase::enumName(enum_id);
	
  // 	os << tag << std::endl 
  // 	   << value << std::endl;
  //       }
  //   }
  
  
  // output the table of function IDs for each property-parameter combination, which is 
  // stored in the table
  //   {
  //     os << "PARAMETER_FUNCTION_ID_TABLE" << std::endl;
  //     unsigned int param_it,prop_it, n_param, n_prop;
  //     n_param = this->parameter_function_ID_table.size(1);
  //     n_prop = this->parameter_function_ID_table.size(2);
  //     for (prop_it = 1; prop_it <= n_prop; prop_it++)
  //       {
  // 	os << 
  // 	for (param_it = 1; param_it <= n_param; param_it++)
  // 	  {
  // 	    os << ;
  // 	      }
  //       }
  //   }
  
  //  os << "MATERIAL_PROPERTY_CARD_DATA_END" << std::endl;
	
  return os;
}



std::istream& operator>> (std::istream& input , 
                          IsotropicMaterialPropertyCard& card)
{
  card.readFromInputStream(input);
  return input; 
}



