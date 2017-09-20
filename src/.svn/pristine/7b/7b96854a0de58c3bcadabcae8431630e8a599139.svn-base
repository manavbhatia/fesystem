// $Id: MultilinearFunction.C,v 1.2 2006-09-05 20:41:45 manav Exp $

// FESystem includes
#include "Numerics/MultilinearFunction.h"


MultilinearFunction::MultilinearFunction():
if_initialized(false)
{
  
}




MultilinearFunction::MultilinearFunction(const std::map<double, double>& values_table):
if_initialized(false)
{
  this->reinit(values_table);
}



MultilinearFunction::~MultilinearFunction()
{
  
}


std::istream& 
MultilinearFunction::readFromInputStream(std::istream& input)
{
  unsigned int n_values = 0;
  double value1 = 0.0, value2 = 0.0;
  
  FESystemIOUtility::readFromInput(input, MULTILINEAR_FUNCTION::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  FESystemIOUtility::readFromInput(input, "N_PAIRS", n_values);
  
  MultilinearFunction::ValueMap value_map;
  
  bool insert_success = 0.0;
  
  for (unsigned int i=0; i < n_values; i++)
    {
    value1 = 0.0; value2 = 0.0;
    input >> value1; input >> value2;
    
    insert_success = value_map.insert(MultilinearFunction::ValueMap::value_type
                                      (value1, value2)).second;
    
    Assert(insert_success,
           MultilinearFunction::ExcDuplicateAbcissa(value1));
    }
  
  // now initialize this function using these value
  this->reinit(value_map);
  
  FESystemIOUtility::readFromInput(input, MULTILINEAR_FUNCTION::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}




std::istream& operator>> (std::istream& input, MultilinearFunction& function)
{
  function.readFromInputStream(input);
  return input;
}



