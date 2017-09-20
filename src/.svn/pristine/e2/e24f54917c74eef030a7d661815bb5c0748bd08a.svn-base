// $Id: PostProcessQty.C,v 1.5 2006-09-05 20:41:50 manav Exp $

// C++ includes
#include <cassert>


// FESystem includes
#include "PostProcess/PostProcessQty.h"


// libMesh includes


PostProcessQty::PostProcessQty()
{

}
  
	

PostProcessQty::~PostProcessQty()
{

}
	

void 
PostProcessQty::addTensor(const std::string& name, const TensorValue<double>& value )
{
  // check if the tensor already exists in the map or not
  std::map<std::string, TensorValue<double> >::const_iterator qty_it = this->name_tensor_map.find(name);
  assert (qty_it == this->name_tensor_map.end());

  bool insert = this->name_tensor_map.insert(std::map<std::string, TensorValue<double> >::value_type(name, value)).second;

  assert (insert == true);
}



TensorValue<double>& 
PostProcessQty::getTensor(const std::string& name)
{
  // check if the tensor already exists in the map or not
  std::map<std::string, TensorValue<double> >::iterator qty_it = this->name_tensor_map.find(name);
  assert (qty_it != this->name_tensor_map.end());

  return (qty_it->second);  
}


