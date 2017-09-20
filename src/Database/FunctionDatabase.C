// $Id: FunctionDatabase.C,v 1.2 2006-09-05 20:41:52 manav Exp $

// C++ includes
#include <string>


// FESystem includes
#include "FESystem/FESystemExceptions.h"
#include "Database/FunctionDatabase.h"
#include "Utilities/InputOutputUtility.h"

#include "Numerics/MultilinearFunction.h"


FunctionDatabase::FunctionDatabase()
{
}


FunctionDatabase::~FunctionDatabase()
{
  // iterate over all the function objects and delete them
  FunctionDatabase::FunctionMapType::iterator it, end;
  it = this->id_to_function_map.begin();
  end = this->id_to_function_map.end();
  
  for (; it != end; it++)
    {
      delete it->second;
      it->second = NULL;
    }
}






std::istream& 
FunctionDatabase::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0, enumID =0;
  
  // beginning of function cards
  FESystemIOUtility::readFromInput(input, "FUNCTION_DATABASE");
  FESystemIOUtility::readFromInput(input, "BEGIN");

  
  // number of function cards
  FESystemIOUtility::readFromInput(input, "N_FUNCTIONS", num);
  
  FunctionBase *func_base = NULL;
  bool insert_success = false;
  
  for (unsigned int i=0; i < num; i++)
    {
    tag.clear();

    // input the type of function card
    FESystemIOUtility::peekFromInput(input, tag);
    
    enumID = FunctionTypeEnum::enumID(tag);
    
    switch (enumID)
      {
      case  MULTILINEAR_FUNCTION_ENUM_ID:
        {
          MultilinearFunction *func = new MultilinearFunction;
          input >> (*func);
          
          func_base = func;
        }
        break;
        
      default:
        {
          Assert(false, FESystemExceptions::ExcEnumCaseNotImplemented(tag));
        }
        break;
      }
    
    // now add the function to the map
    insert_success = this->id_to_function_map.insert(FunctionDatabase::FunctionMapType::value_type
                                                     (func_base->getID(), func_base)).second;
    }
  
  // make sure it ends right  
  FESystemIOUtility::readFromInput(input, "FUNCTION_DATABASE");
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}






std::istream& operator>>(std::istream& input, FunctionDatabase& obj)
{
  obj.readFromInputStream(input);
  return input;
}






//void 
//FunctionDatabase::getScalarFunctionValues
//(const unsigned int func_ID, std::vector<double>& values) const
//{
//
//}
//  




//void 
//FunctionDatabase::getScalarValues(const std::vector<unsigned int>& func_IDs, 
//                                  const double abcissa,
//                                  std::vector<double>& values)  const
//{
//  
//}
//



//void 
//FunctionDatabase::getScalarDerivatives(const unsigned int func_ID, std::vector<double>& values) const
//{
//  
//}
//  
//


//
//void 
//FunctionDatabase::getScalarDerivatives(const std::vector<unsigned int>& func_IDs, 
//				       const double abcissa,
//				       std::vector<double>& values)  const
//{
//
//}

