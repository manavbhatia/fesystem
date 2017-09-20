// $Id: FunctionDatabase.h,v 1.2 2006-09-05 20:41:52 manav Exp $

#ifndef __fesystem_function_database_h__
#define __fesystem_function_database_h__

// C++ includes
#include <iostream>
#include <vector>
#include <map>


// FESystem inlcudes
#include "Numerics/FunctionBase.h"

///
/// this class provides the functionality of an interface to the different function
/// classes, in deal.II and FESystem. This is needed since deal.II has multiple classes
/// for functions, i.e. Polynomial and Function, which provide the same basic functionality.
/// One way would have been to create a base class named FunctionBase and derive all these
/// classes out of that so that all of them could be stored as a single pointer. However, 
/// they are implemented with different method names and hence, inheritence would have 
/// been confusing. Hence, this helper class was conceived to hangle the hetrogeniety 
/// in these functions, but still provide a uniform interface. This class is more like 
/// a database of functions, and will store all functions with an ID associated with them. 
/// Then, each function can be accessed from this database, or the interface of this class
/// can be used to access the methods of each individual class
///

class FunctionDatabase
{
 public:
  
  /// constructor
  FunctionDatabase();

  /// Destructor
  ~FunctionDatabase();

  
  /// returns the value of the function whose ID has been specified. The coordinates 
  /// where the function needs to be evaluated is also specified in the input. This method 
  /// returns the calculated data in created objects, and hence, has to create it each time
  /// the function is called. Hence, to save on computational cost, the more specialized functions 
  /// in the API should be called.
  inline void getScalarFunctionValue(const unsigned int func_ID, const double input,
                                     double& func_value) const;

//  /// this method calculates the value of a function at multiple points and returns it in the 
//  /// same vector that contained the input values. This works only for uni-component functions, 
//  /// in 1-D spaces, i.e. f = f(x). This method is advantageous when the database will be called
//  /// multiple times to obtain funciton values, since the same vector can then be used to store 
//  /// abcissa and ordinate values.
//  void getScalarFunctionValues(const unsigned int func_ID, std::vector<double>& values) const;
//  
//  /// this method calculates the value of the multiple functions at the given value and returns
//  /// the calculated values in vector, that is passed as input. This works only for uni-components
//  /// in 1-D spaces, i.e. f = f(x)
//  void getScalarFunctionValues(const std::vector<unsigned int>& func_IDs, 
//                               const double abcissa,
//                               std::vector<double>& values)  const;
  
  /// returns the derivative value of the function whose ID has been specified. The coordinates 
  /// where the function derivative needs to be evaluated is also specified in the input. This method 
  /// returns the calculated data in created objects, and hence, has to create it each time
  /// the function is called. Hence, to save on computational cost, the more specialized functions 
  /// in the API should be called.
  inline void getScalarFunctionDerivative(const unsigned int func_ID, 
                                     const double input, double& func_value) const;
  

//  /// this method calculates the value of a function derivative at multiple points and returns it in the 
//  /// same vector that contained the input values. This works only for uni-component functions, 
//  /// in 1-D spaces, i.e. f = f(x). This method is advantageous when the database will be called
//  /// multiple times to obtain funciton values, since the same vector can then be used to store 
//  /// abcissa and ordinate values.
//  void getScalarFunctionDerivatives(const unsigned int func_ID, std::vector<double>& values) const;
//  
//  /// this method calculates the derivative value of the multiple functions at the given value and returns
//  /// the calculated values in vector, that is passed as input. This works only for uni-components
//  /// in 1-D spaces, i.e. f = f(x)
//  void getScalarFunctionDerivatives(const std::vector<unsigned int>& func_IDs, 
//			    const double abcissa,
//			    std::vector<double>& values)  const;


  /// this method returns a constant reference to the function base object defined by the ID
  inline const FunctionBase& getFunction(const unsigned int ID) const;

  /// reads from the input stream 
  std::istream& readFromInputStream(std::istream& input);

  /// this operator reads the card from an input stream
  friend std::istream&  operator>>(std::istream& input, FunctionDatabase& obj);

 protected:
  
  /// this exception is thrown if an ID requested by the user does not exist in the database
  DeclException1(ExcFunctionIDDoesNotExist, unsigned int, 
                 << " Requested ID: " << arg1 << " does not exist in FunctionDatabase");
  
  /// local type definitions
  typedef std::map<unsigned int, FunctionBase*>  FunctionMapType;
  
  /// map that stores the functions against their ID
  FunctionMapType id_to_function_map;
};



inline
const FunctionBase& 
FunctionDatabase::getFunction(const unsigned int ID) const
{
  // make sure that this ID is valid
  Assert(ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(ID));
  
  // make sure that this function ID exists in the map, and if it does, return it
  FunctionDatabase::FunctionMapType::const_iterator it, end;
  it = this->id_to_function_map.find(ID);
  end = this->id_to_function_map.end();
  
  Assert(it != end,
         FunctionDatabase::ExcFunctionIDDoesNotExist(ID));
  
  return *(it->second);
}




inline
void
FunctionDatabase::getScalarFunctionValue
(const unsigned int func_ID, const double input,
 double& value) const
{
  const FunctionBase & func_base = this->getFunction(func_ID);
  func_base.getFunctionValue(input, value);
}



inline
void
FunctionDatabase::getScalarFunctionDerivative
(const unsigned int func_ID, const double input,
 double& value) const
{
  const FunctionBase & func_base = this->getFunction(func_ID);
  
  func_base.getFunctionDerivative(input, value);
}


#endif // __fesystem_function_database_h__
