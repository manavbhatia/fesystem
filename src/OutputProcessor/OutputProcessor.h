// $Id: OutputProcessor.h,v 1.4 2006-09-05 20:41:42 manav Exp $

#ifndef __fesystem_output_processor_h__
#define __fesystem_output_processor_h__

// C++ includes
#include <string>
#include <memory>

// FESystem includes
#include "Utilities/NameEnumHandler.h"

// Forward decleration
namespace FESystem
{
  class FESystemController;
}

class OutputInfo;

DeclareEnumClass(OutputFormatEnum);



class OutputProcessor
{
public:	
  
	OutputProcessor(FESystem::FESystemController& ,
                  const OutputInfo& info,
                  const unsigned int output_enum_ID);
	virtual ~OutputProcessor();
	
	virtual void writeData()=0;
  
	inline unsigned int getOutputFormatEnumID() const;

	inline std::string getOutputFormatEnumName() const;

  static std::auto_ptr<OutputProcessor> 
    createOutputProcessor(FESystem::FESystemController& controller, 
                          const OutputInfo& info);
  
protected:	
		
    FESystem::FESystemController& fesystem_controller;
	const OutputInfo &output_info;
	unsigned int output_format_enum_ID;
};



inline 
unsigned int
OutputProcessor::getOutputFormatEnumID() const
{
  return this->output_format_enum_ID;
}



inline 
std::string
OutputProcessor::getOutputFormatEnumName() const
{
  return OutputFormatEnum::enumName(this->output_format_enum_ID);
}


#endif // __fesystem_output_processor_h__
