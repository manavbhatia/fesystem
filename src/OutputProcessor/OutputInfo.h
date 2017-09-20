// $Id: OutputInfo.h,v 1.3 2006-09-05 20:41:42 manav Exp $

#ifndef __fesystem_output_info_h__
#define __fesystem_output_info_h__

// C++ includes
#include <iostream>
#include <string>


// FESystem includes
#include "Utilities/NameEnumHandler.h"
#include "OutputProcessor/OutputProcessor.h"


class OutputInfo
{
public:
  OutputInfo();
  ~OutputInfo();

  inline unsigned int getOutputFormatEnumID() const;
  inline std::string getOutputFormatEnumName() const;
  inline const std::string& getOutputFileName() const;
  
  
  std::istream& readFromInputStream(std::istream& input);
  
  friend std::istream& operator >> (std::istream& input, OutputInfo& info);
  
protected:
  
  unsigned int output_format_enum_ID;
  std::string output_file_name;
};


inline
unsigned int
OutputInfo::getOutputFormatEnumID() const
{
  return this->output_format_enum_ID;
}



inline
std::string
OutputInfo::getOutputFormatEnumName() const
{
  return OutputFormatEnum::enumName(this->output_format_enum_ID);
}



inline
const std::string& 
OutputInfo::getOutputFileName() const
{
  return this->output_file_name;
}




#endif //__fesystem_output_info_h__

