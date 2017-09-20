// $Id: MultilinearFunction.h,v 1.3.6.1 2008-04-06 04:02:41 manav Exp $

#ifndef __fesystem_multilinear_function_h__
#define __fesystem_multilinear_function_h__

// C++ includes
#include <vector>
#include <map>
#include <cmath> 

// FESystem includes
#include "FESystem/FESystemExceptions.h"
#include "Numerics/FunctionBase.h"
#include "Utilities/InputOutputUtility.h"


/// define the enumeration of this function type

#ifndef MULTILINEAR_FUNCTION_ENUM_ID
#define MULTILINEAR_FUNCTION_ENUM_ID 1
#else
#error 
#endif

#ifndef MULTILINEAR_FUNCTION_ENUM_NAME
#define MULTILINEAR_FUNCTION_ENUM_NAME "MULTILINEAR_FUNCTION"
#else
#error
#endif

DeclareEnumName(MULTILINEAR_FUNCTION, FunctionTypeEnum, 
                MULTILINEAR_FUNCTION_ENUM_ID, 
                MULTILINEAR_FUNCTION_ENUM_NAME);

/// this class defines a 1-D function that is an assortment of linear functions 
/// over adjacent intervals. This is used to define, for example, the dependence
/// of a material property on a parameter, which is available as a table of ordinate
/// values at discrete abssicas. This method asserts that only one function value 
/// be specified at a given abcissa
class MultilinearFunction: public FunctionBase
{
public:

  
  /// constructor. This is the default constructor, and if this is used, then
  /// the function should be reinitialized with the data table before calling 
  /// any of the methods to obtain the values from the function.
  MultilinearFunction();

  /// constructor. This takes the table of values as input and initializes the
  /// data structures, so that methods to obtain the values of the function 
  /// can be called directly after this
  MultilinearFunction(const std::map<double, double>& values_table);

  /// destructor
  virtual ~MultilinearFunction();

  /// returns the type of this function as the enumeration ID of the enumeration type of this
  /// object
  inline const unsigned int getFunctionTypeID() const;

  /// returns the type of this function as a string
  inline const std::string getFunctionTypeName() const;
  
  /// this method re-initializes this object to the new data set. The data should 
  /// have atleast two data points
  inline void reinit(const std::map<double, double>& values_table);

  /// method returns the value of the function at the specified point. If the value
  /// of abcissa in the input parameter is outside the bounds of where the function is 
  /// defined, then linear extrapolation is used.
  inline void getFunctionValue(const double& x,
                               double& ordinate) const;


  /// @returns derivative of the function at the given abcissa 
  /// @param abcissa value at which the function should be evaluated 
  /// @param ordinate evaluated function value 
  inline void getFunctionDerivative(const double& x,
                                    double& ordinate) const;
  
  
  
  /// overloaded input stream operator
  friend std::istream& operator>> (std::istream& input, MultilinearFunction& function);
  
protected:

  /// method to read from the function from the input stream
  virtual std::istream& readFromInputStream(std::istream& input);

  /// this exception will be thrown if two same abcissa values are specified 
  /// in the input
  DeclException1(ExcDuplicateAbcissa, double, 
                 << "Duplicate abcissa value: " << arg1 
                 << " specified for the function");
    
  /// if the method has been initialized or not
  bool if_initialized;
  
  /// local type definitions
  typedef std::map<double, double> ValueMap;
  
  /// table of abcissa-ordinate values for this function
  ValueMap value_table;
};








inline const unsigned int 
MultilinearFunction::getFunctionTypeID() const
{
  return MULTILINEAR_FUNCTION::num();
}



inline const std::string 
MultilinearFunction::getFunctionTypeName() const
{
  return MULTILINEAR_FUNCTION::name();
}



inline void 
MultilinearFunction::reinit(const std::map<double, double>& values_table)
{
  unsigned int n_points = values_table.size();
  // perform a sanity check on the data. 
  Assert(n_points > 1, ExcLowerRange(n_points, 2));

  this->value_table.clear();
  this->value_table = values_table;
  this->if_initialized = true;
}



inline 
void
MultilinearFunction::getFunctionValue(const double& x,
                                      double& func_value) const
{
  Assert (this->if_initialized ,ExcNotInitialized());

  // find the iterator whose abcissa value is less than the given abcissa value
  MultilinearFunction::ValueMap::const_iterator it, end;
  MultilinearFunction::ValueMap::const_reverse_iterator rit;
  it = this->value_table.begin();
  end = this->value_table.end();
  rit = this->value_table.rbegin();

  double lower_x = 0.0 , upper_x = 0.0, lower_y = 0.0, upper_y = 0.0, delta_x = 0.0;

  if ( x <= it->first)
    {
      lower_x = it->first; lower_y = it->second;
      it++;
      upper_x = it->first; upper_y = it->second;      
    }
  else if (x > rit->first)
    {
      upper_x = rit->first; upper_y = rit->second;
      rit++;
      lower_x = rit->first; lower_y = rit->second;
    }
  else
    {
      // now find the iterator whose value is greater than the given abcissa
      it++;
      for (; it != end; it++)
	{
	  if (x <= it->first)
	    // x lies between the previous and this iterator
	    {
	      upper_x = it->first; upper_y = it->second;
	      it--;
	      lower_x = it->first; lower_y = it->second;
        break;
      }
	}
    }

  // now calculate the ordinate value and return
  delta_x = upper_x-lower_x;
  if (fabs(delta_x) <= FESystemNumbers::Epsilon)
    func_value = lower_y;
  else
    func_value = (lower_y + (upper_y-lower_y)/delta_x * (x-lower_x));
}





inline 
void
MultilinearFunction::getFunctionDerivative(const double& x,
                                      double& func_derivative) const
{
  Assert(this->if_initialized ,ExcNotInitialized());
  
  // find the iterator whose abcissa value is less than the given abcissa value
  MultilinearFunction::ValueMap::const_iterator it, end;
  MultilinearFunction::ValueMap::const_reverse_iterator rit;
  it = this->value_table.begin();
  end = this->value_table.end();
  rit = this->value_table.rbegin();
  
  double lower_x = 0.0 , upper_x = 0.0, lower_y = 0.0, upper_y = 0.0, delta_x = 0.0;
  
  if ( x <= it->first)
    {
    lower_x = it->first; lower_y = it->second;
    it++;
    upper_x = it->first; upper_y = it->second;      
    }
  else if (x > rit->first)
    {
    upper_x = rit->first; upper_y = rit->second;
    rit++;
    lower_x = rit->first; lower_y = rit->second;
    }
  else
    {
    // now find the iterator whose value is greater than the given abcissa
    it++;
    for (; it != end; it++)
      {
      if (x <= it->first)
        // x lies between the previous and this iterator
        {
	      upper_x = it->first; upper_y = it->second;
	      it--;
	      lower_x = it->first; lower_y = it->second;
          break;
        }
      }
    }
  
  // now calculate the ordinate value and return
  delta_x = upper_x-lower_x;
  Assert(fabs(delta_x) > FESystemNumbers::Epsilon,
         ExcDivideByZero());
    
  func_derivative = (upper_y-lower_y)/delta_x;
}


#endif // __fesystem_multilinear_function_h__
