// $Id: TecPlotOutputProcessor.h,v 1.5.6.1 2007-03-14 22:05:03 manav Exp $

#ifndef __fesystem_tecplot_output_processor_h__
#define __fesystem_tecplot_output_processor_h__

// C++ includes
#include <string>
#include <map>
#include <iostream>
#include <set>

// FESystem includes
#include "OutputProcessor/OutputProcessor.h"

// libmesh includes
#include "mesh/mesh.h"

namespace MeshDS
{
  class FEMeshData;
  class FEMesh;
}
namespace FESystem
{
  class FESystemController;
}



#ifndef TECPLOT_OUTPUT_PROCESSOR_ENUM_ID
#define TECPLOT_OUTPUT_PROCESSOR_ENUM_ID 2
#else
#error
#endif

#ifndef TECPLOT_OUTPUT_PROCESSOR_ENUM_NAME
#define TECPLOT_OUTPUT_PROCESSOR_ENUM_NAME "TECPLOT_OUTPUT_PROCESSOR"
#else
#error
#endif


DeclareEnumName(TECPLOT_OUTPUT_PROCESSOR, OutputFormatEnum,
                TECPLOT_OUTPUT_PROCESSOR_ENUM_ID, TECPLOT_OUTPUT_PROCESSOR_ENUM_NAME);


class TecPlotOutputProcessor: public OutputProcessor
{
 public:	
  
  /// constructor
  TecPlotOutputProcessor(FESystem::FESystemController& ,
                         const OutputInfo& info);

  
  /// destructor
  virtual ~TecPlotOutputProcessor();
	
  /// writes data to output file(s)
  virtual void writeData();
	
	
 protected:	
  /// prepares the elem type map
  void prepareElemTypeMap();

  /// writes the thermal data
  void writeNodalData(std::ostream& output, const unsigned int discipline, 
		      const std::vector<std::string>& dof_names);
	
  // writes the structural data 
  void writeStructuralData(std::ostream& output);

  MeshDS::FEMesh* mesh;
  MeshDS::FEMeshData* mesh_data;
  std::set<Elem*>* dummy_elem_set;
  std::map<ElemType, std::vector<Elem*> > elem_vec_map;
	
};


#endif // __fesystem_tecplot_output_processor_h__
