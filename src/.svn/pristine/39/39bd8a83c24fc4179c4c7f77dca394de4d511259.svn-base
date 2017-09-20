// $Id: GmshOutputProcessor.h,v 1.1.6.1 2007-03-14 22:05:03 manav Exp $

#ifndef __fesystem_gmsh_output_processor_h__
#define __fesystem_gmsh_output_processor_h__

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

namespace FESystemDatabase
{
  class DataInfoBase;
}

template<class T> class NumericVector;

#ifndef GMSH_OUTPUT_PROCESSOR_ENUM_ID
#define GMSH_OUTPUT_PROCESSOR_ENUM_ID 3
#else
#error
#endif

#ifndef GMSH_OUTPUT_PROCESSOR_ENUM_NAME
#define GMSH_OUTPUT_PROCESSOR_ENUM_NAME "GMSH_OUTPUT_PROCESSOR"
#else
#error
#endif


DeclareEnumName(GMSH_OUTPUT_PROCESSOR, OutputFormatEnum,
                GMSH_OUTPUT_PROCESSOR_ENUM_ID, GMSH_OUTPUT_PROCESSOR_ENUM_NAME);


class GmshOutputProcessor: public OutputProcessor
{
public:	
  
  /// constructor
  GmshOutputProcessor(FESystem::FESystemController& ,
                         const OutputInfo& info);
  
  
  /// destructor
  virtual ~GmshOutputProcessor();
	
  /// writes data to output file(s)
  virtual void writeData();
	
	
protected:	
    
    /// @returns the name of the element
    std::string getElemTypeName(ElemType type, 
                                const unsigned int tensor_order);
    
  void prepareElemTypeMap();
  
  void writeViewHeader(std::ostream& output, 
                       const unsigned int n_scalars,
                       const unsigned int n_vectors,
                       const unsigned int n_tensors);
    
  /// writes the nodal data
  void writeNodalData(std::ostream& output,
                      const unsigned int discipline,
                      const unsigned int tensor_order);
  
  void writeSolutionVector(std::ostream& output,
                           const std::vector<NumericVector<double>*>* sol_vec,
			   const unsigned int n_dofs);
  
  MeshDS::FEMesh* mesh;
  MeshDS::FEMeshData* mesh_data;
  std::set<Elem*>* dummy_elem_set;
  std::map<ElemType, std::vector<Elem*> > elem_vec_map;
};


#endif // __fesystem_tecplot_output_processor_h__
