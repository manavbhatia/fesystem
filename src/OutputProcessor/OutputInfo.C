// $Id: OutputInfo.C,v 1.2 2006-09-05 20:41:42 manav Exp $

// FESystem includes
#include "OutputProcessor/OutputInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "FESystem/FESystemNumbers.h"


OutputInfo::OutputInfo():
output_format_enum_ID(FESystemNumbers::InvalidID),
output_file_name()
{
  
}



OutputInfo::~OutputInfo()
{
  
}



std::istream& 
OutputInfo::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  FESystemIOUtility::readFromInput(input, "OUTPUT_INFO");
  FESystemIOUtility::readFromInput(input, "BEGIN");

  FESystemIOUtility::readFromInput(input, "OUTPUT_FORMAT", tag);
  this->output_format_enum_ID = OutputFormatEnum::enumID(tag);
  
  FESystemIOUtility::readFromInput(input, "OUTPUT_FILE_NAME", this->output_file_name);
  
  FESystemIOUtility::readFromInput(input, "OUTPUT_INFO");
  FESystemIOUtility::readFromInput(input, "END");

  return input;
}




std::istream& 
operator >> (std::istream& input, OutputInfo& info)
{
  info.readFromInputStream(input);
  return input;
}

