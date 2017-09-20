/*
 *  VTKOutputProcessor.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 11/16/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

#ifndef __vtk_output_processor_h__
#define __vtk_output_processor_h__


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

#ifndef VTK_OUTPUT_PROCESSOR_ENUM_ID
#define VTK_OUTPUT_PROCESSOR_ENUM_ID 4
#else
#error
#endif

#ifndef VTK_OUTPUT_PROCESSOR_ENUM_NAME
#define VTK_OUTPUT_PROCESSOR_ENUM_NAME "VTK_OUTPUT_PROCESSOR"
#else
#error
#endif


DeclareEnumName(VTK_OUTPUT_PROCESSOR, OutputFormatEnum,
                VTK_OUTPUT_PROCESSOR_ENUM_ID, VTK_OUTPUT_PROCESSOR_ENUM_NAME);


class VTKOutputProcessor: public OutputProcessor
  {
  public:	
    
    /// constructor
    VTKOutputProcessor(FESystem::FESystemController& ,
                        const OutputInfo& info);
    
    
    /// destructor
    virtual ~VTKOutputProcessor();
    
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


#endif // __vtk_output_processor_h__

