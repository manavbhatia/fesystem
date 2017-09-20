// $Id: PostProcessQty.h,v 1.5 2006-09-05 20:41:50 manav Exp $

#ifndef __post_process_qty_h__
#define __post_process_qtyh_h__

// C++ include
#include <map>
#include <string>


// FESystem include

// libMesh include
#include "numerics/tensor_value.h"


class PostProcessQty
{
public:
  
  /// constructor
  /// @param order of the tensor
  PostProcessQty();
  
	
  /// destructor
  virtual ~PostProcessQty();
	
  /// add the value of tensor, stored by name
  /// @param name of the tensor
  /// @param value of the tensor
  void addTensor(const std::string& name, const TensorValue<double>& value );

  /// gets the value of the tensor
  TensorValue<double>& getTensor(const std::string& name);


protected:
  
  /// map in which all the tensors are stored
  std::map<std::string, TensorValue<double> >  name_tensor_map;
};


#endif // __post_process_qty_h__
