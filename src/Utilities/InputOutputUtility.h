// $Id: InputOutputUtility.h,v 1.4 2006-09-10 05:48:06 manav Exp $

#ifndef __fesystem_input_output_utility_h__
#define __fesystem_input_output_utility_h__

// C++ includes
#include <iostream>
#include <string>

// FESystem includes
#include "FESystem/FESystemExceptions.h"

namespace FESystemIOUtility
{
  /// makes sure that the input has the specified tag
  inline std::istream& readFromInput(std::istream& input, 
                                     const std::string& tag);
  
  /// reads a string from input and stores it in the input parameter.
  inline std::istream& readFromInput(std::istream& input, 
                                     const std::string& tag,
                                     std::string& str);
  
  /// reads an integer from input and stores it in the input parameter
  inline std::istream& readFromInput(std::istream& input, 
                                     const std::string& tag,
                                     unsigned int& num);

  /// reads a real from input and stores it in the input parameter
  inline std::istream& readFromInput(std::istream& input, 
                                     const std::string& tag,
                                     double& num);
  
  /// reads a bool from input and stores it in the input parameter
  inline std::istream& readFromInput(std::istream& input, 
                                     const std::string& tag,
                                     bool& val);

  /// peeks the next string from the input. The function returns the input 
  /// stream whose seek pointer is unchanged 
  inline std::istream& peekFromInput(std::istream& input, 
                                     std::string& tag);

}



inline 
std::istream& 
FESystemIOUtility::readFromInput(std::istream& input, 
                                 const std::string& tag)
{
  static std::string in_tag;
  in_tag.clear();
  
  // read the input and make sure that it is same as the 
  // tag specified in the input
  input >> in_tag;
  Assert(in_tag == tag,
         FESystemExceptions::ExcIOBadTag(in_tag, tag));
  
  return input;
}



inline 
std::istream& 
FESystemIOUtility::readFromInput(std::istream& input, 
                                 const std::string& tag,
                                 std::string& str)
{
  static std::string in_tag;
  in_tag.clear();

  // read the input and make sure that it is same as the 
  // tag specified in the input
  input >> in_tag;
  Assert(in_tag == tag,
         FESystemExceptions::ExcIOBadTag(in_tag, tag));

  // now read the next value
  input >> str;

  return input;
}



inline 
std::istream& 
FESystemIOUtility::readFromInput(std::istream& input, 
                                 const std::string& tag,
                                 unsigned int& num)
{
  static std::string in_tag;
  in_tag.clear();
  
  // read the input and make sure that it is same as the 
  // tag specified in the input
  input >> in_tag;
  Assert(in_tag == tag,
         FESystemExceptions::ExcIOBadTag(in_tag, tag));
  
  // now read the next value
  input >> num;
  
  return input;
}



inline 
std::istream& 
FESystemIOUtility::readFromInput(std::istream& input, 
                                 const std::string& tag,
                                 double& num)
{
  static std::string in_tag;
  in_tag.clear();
  
  // read the input and make sure that it is same as the 
  // tag specified in the input
  input >> in_tag;
  Assert(in_tag == tag,
         FESystemExceptions::ExcIOBadTag(in_tag, tag));
  
  // now read the next value
  input >> num;
  
  return input;
}


inline 
std::istream& 
FESystemIOUtility::readFromInput(std::istream& input, 
                                 const std::string& tag,
                                 bool& val)
{
  static std::string in_tag;
  in_tag.clear();
  
  // read the input and make sure that it is same as the 
  // tag specified in the input
  input >> in_tag;
  Assert(in_tag == tag,
         FESystemExceptions::ExcIOBadTag(in_tag, tag));
  
  // now read the next value
  in_tag.clear();
  input >> in_tag;
  if ((in_tag == "TRUE") || (in_tag == "true") || (in_tag == "True"))
    val = true;
  else if ((in_tag == "FALSE") || (in_tag == "false") || (in_tag == "False"))
    val = false;
  else 
    Assert(false, 
           FESystemExceptions::ExcIOBadTag(in_tag, "TRUE/FALSE"));
  
  return input;
}



inline 
std::istream& 
FESystemIOUtility::peekFromInput(std::istream& input, 
                                 std::string& tag)
{
  tag.clear();
  
  std::streampos pos = input.tellg();
  
  input >> tag;
  
  input.seekg(pos);
  
  return input;
}


#endif // __fesystem_input_output_utility_h__

