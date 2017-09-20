// $Id: FunctionBase.h,v 1.2 2006-09-05 20:41:45 manav Exp $

#ifndef __fesystem_function_base_h__
#define __fesystem_function_base_h__

// C++ includes
#include <iostream>
#include <string>


// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"
#include "Utilities/NameEnumHandler.h"



/// an enumeration class is declared to handle the function type ennumerations
DeclareEnumClass(FunctionTypeEnum);

/// this class defines a base class for all functions in both deal.II and FESystem, 
/// so that the can be
/// referred to by a pointer to this base class
class FunctionBase
{
public:
  
  /// default constructor
  inline FunctionBase();
  
  /// destructor
  virtual inline ~FunctionBase();
  
  /// @returns the function ID
  inline const unsigned int getID() const;
  
  /// @returns the type of this function as the enumeration ID of the enumeration type of this
  /// object
  virtual const unsigned int getFunctionTypeID() const = 0;
  
  /// @returns the type of this function as a string
  virtual const std::string getFunctionTypeName() const = 0;
  
  /// @returns value of the function at the given abcissa 
  /// @param abcissa value at which the function should be evaluated
  /// @param ordinate evaluated function value
  virtual void getFunctionValue(const double& abcissa,
                                double& ordinate) const = 0;

  
  /// @returns derivative of the function at the given abcissa 
  /// @param abcissa value at which the function should be evaluated 
  /// @param ordinate evaluated function value 
  virtual void getFunctionDerivative(const double& abcissa,
                                     double& ordinate) const = 0;
  
  
protected:
    
  /// method to read from the function from the input stream
  virtual std::istream& readFromInputStream(std::istream& input) = 0;
    
    
  /// ID of this function object
  unsigned int ID;
  
};





inline 
FunctionBase::FunctionBase():
ID(FESystemNumbers::InvalidID)
{}





inline 
FunctionBase::~FunctionBase()
{}




inline const unsigned int 
FunctionBase::getID() const
{
  Assert (this->ID > 0, FESystemExceptions::ExcInvalidID(this->ID));
  return this->ID;
}


#endif // __fesystem_function_base_h__
