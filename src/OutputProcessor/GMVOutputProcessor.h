// $Id: GMVOutputProcessor.h,v 1.6 2006-09-05 20:41:42 manav Exp $

#ifndef __fesystem_gmv_output_processor_h__
#define __fesystem_gmv_output_processor_h__

// C++ includes
#include <string>
#include <map>
#include <iostream>
#include <set>

// FESystem includes
#include "OutputProcessor/OutputProcessor.h"

// libmesh inclues
#include "mesh.h"

namespace MeshDS
{
  class FEMeshData;
  class FEMesh;
}


#ifndef GMV_OUTPUT_PROCESSOR_ENUM_ID
#define GMV_OUTPUT_PROCESSOR_ENUM_ID 1
#else
#error
#endif

#ifndef GMV_OUTPUT_PROCESSOR_ENUM_NAME
#define GMV_OUTPUT_PROCESSOR_ENUM_NAME "GMV_OUTPUT_PROCESSOR"
#else
#error
#endif


DeclareEnumName(GMV_OUTPUT_PROCESSOR, OutputFormatEnum,
                GMV_OUTPUT_PROCESSOR_ENUM_ID, GMV_OUTPUT_PROCESSOR_ENUM_NAME);


class GMVOutputProcessor: public OutputProcessor
{
 public:	
  GMVOutputProcessor(FESystem::FESystemController& ,
                     const OutputInfo& info);
  virtual ~GMVOutputProcessor();
	
  virtual void writeData();
	
	
 protected:	
	
  void writeHeader(std::ostream& );
  void writeNodes(std::ostream& );
  void writeElems(std::ostream& );
  void writeDisplayNodeIDs(std::ostream& );
  void writeDisplayElemIDs(std::ostream& );
  void writeVariables(std::ostream& , const unsigned int discipline_enum_ID);
  void writeFooter(std::ostream& );
	
  MeshDS::FEMesh* mesh;
  MeshDS::FEMeshData* mesh_data;
  std::set<Elem*>* dummy_elem_set;
	
};


#endif // __fesystem_gmv_output_processor_h__
