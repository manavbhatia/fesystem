// $Id: OutputProcessor.C,v 1.7 2006-09-25 08:19:54 manav Exp $



//FESystem includes
#include "OutputProcessor/OutputProcessor.h"
#include "FESystem/FESystemController.h"
#include "OutputProcessor/GMVOutputProcessor.h"
#include "OutputProcessor/TecPlotOutputProcessor.h"
#include "OutputProcessor/GmshOutputProcessor.h"
#include "OutputProcessor/OutputInfo.h"

OutputProcessor::OutputProcessor(FESystem::FESystemController& fesys_controller,
                                 const OutputInfo& info,
                                 const unsigned int output_enum_ID):
fesystem_controller(fesys_controller),
output_info(info),
output_format_enum_ID(output_enum_ID)
{
	
}



OutputProcessor::~OutputProcessor()
{
	
}




std::auto_ptr<OutputProcessor> 
OutputProcessor::createOutputProcessor(FESystem::FESystemController& controller, 
                                       const OutputInfo& info)
{
  std::auto_ptr<OutputProcessor> processor;
  
  switch (info.getOutputFormatEnumID())
    {
    case GMV_OUTPUT_PROCESSOR_ENUM_ID:
      processor.reset(new GMVOutputProcessor(controller, info));
      break;
      
    case TECPLOT_OUTPUT_PROCESSOR_ENUM_ID:
      processor.reset(new TecPlotOutputProcessor(controller, info));
      break;
      
    case GMSH_OUTPUT_PROCESSOR_ENUM_ID:
      processor.reset(new GmshOutputProcessor(controller, info));
      break;

    default:
      abort();
    }

  return processor;
}
