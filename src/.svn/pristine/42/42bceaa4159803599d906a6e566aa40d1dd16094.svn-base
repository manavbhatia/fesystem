// $Id: FEMeshData.C,v 1.6 2006-09-05 20:41:34 manav Exp $

// C++ includes


// FESystem inlucdes
#include "FEMeshData.h"


MeshDS::FEMeshData::FEMeshData(const MeshDS::FEMesh& m):
mesh(m)
{
	
}



MeshDS::FEMeshData::~FEMeshData()
{
  
}








void 
MeshDS::FEMeshData::addElemAndIDData(const Elem* elem, 
                                  const unsigned int foreign_elem_ID,
                                  const unsigned int property_ID,
                                  const unsigned int elem_kind_enum_ID)
{
  Assert(elem != NULL,
         ExcEmptyObject());
  Assert(foreign_elem_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(foreign_elem_ID));
  

  // create the entries in the differet maps and vectors.
  // first get the internal ID for this elem.
  unsigned int internal_ID = this->elem_vector.size();
  
  bool insert_success = 
    this->foreign_ID_to_elem_map.insert(MeshDS::FEMeshData::IDToElemMap::value_type
                                        (foreign_elem_ID, elem)).second;
  // if the insert failed, then the ID is not unique
  Assert(insert_success,
         MeshDS::FEMeshData::ExcDuplicateID("Elem", foreign_elem_ID));
	
  // next create the ID	data structure
  MeshDS::FEMeshData::ElemPropIDs elem_ids;
  elem_ids.foreign_ID = foreign_elem_ID;
  elem_ids.internal_ID = internal_ID;
  elem_ids.property_ID = property_ID;
  elem_ids.kind_enum_ID = elem_kind_enum_ID;
  
  // insert this in the map
  insert_success = 
    this->elem_ID_map.insert(MeshDS::FEMeshData::ElemToIDMap::value_type
                                        (elem, elem_ids)).second;

  // if the insert failed, then the elem is not unique
  Assert(insert_success,
         MeshDS::FEMeshData::ExcDuplicatePointer("Elem"));

  // finally, insert the elem in the vector
  this->elem_vector.push_back(elem);
}





void 
MeshDS::FEMeshData::addNodeAndIDData(const Node* node, 
                                  const unsigned int foreign_node_ID)
{
  Assert(node != NULL,
         ExcEmptyObject());
  Assert(foreign_node_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(foreign_node_ID));
  
  
  // create the entries in the differet maps and vectors.
  // first get the internal ID for this elem.
  unsigned int internal_ID = this->node_vector.size();
  
  bool insert_success = 
    this->foreign_ID_to_node_map.insert(MeshDS::FEMeshData::IDToNodeMap::value_type
                                        (foreign_node_ID, node)).second;
  // if the insert failed, then the ID is not unique
  Assert(insert_success,
         MeshDS::FEMeshData::ExcDuplicateID("Node", foreign_node_ID));
	
  // next create the ID	data structure
  NodeIDs node_ids;
  node_ids.foreign_ID = foreign_node_ID;
  node_ids.internal_ID = internal_ID;
  
  // insert this in the map
  insert_success = 
    this->node_ID_map.insert(MeshDS::FEMeshData::NodeToIDMap::value_type
                                        (node, node_ids)).second;
  
  // if the insert failed, then the elem is not unique
  Assert(insert_success,
         MeshDS::FEMeshData::ExcDuplicatePointer("Node"));
  
  // finally, insert the elem in the vector
  this->node_vector.push_back(node);
}




void
MeshDS::FEMeshData::writeMeshDataToStream(std::ostream& output)
{
  // iterate over each external node ID, get the node from the mesh data for that
  // and write its data
  IDToNodeMap::const_iterator it, end;
  it = this->foreign_ID_to_node_map.begin();
  end = this->foreign_ID_to_node_map.end();
  
  for ( ; it != end; it++)
    {
      output << it->first << "  " 
      << (*(it->second))(0) << " "
      << (*(it->second))(1) << " "
      << (*(it->second))(2) << " ";
      // now the dof indices
      // system is assumed to be zero
      for (unsigned int i=0; i<it->second->n_vars(0); i++)
        output << it->second->dof_number(0,i,0) << "  ";
      output << std::endl;
    }
  
}




