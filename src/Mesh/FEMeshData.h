// $Id: FEMeshData.h,v 1.8 2006-09-05 20:41:34 manav Exp $

#ifndef __fesystem_fe_mesh_data_h__
#define __fesystem_fe_mesh_data_h__

// C++ includes
#include <iostream>
#include <map>
#include <vector>
#include <string>


// FESystem inlucdes
#include "FESystem/FESystemExceptions.h"
#include "FESystem/FESystemElem.h"


// Forward declerations
class Elem;
class Node;

/// this class provides the necessary data structure to store the elements and nodes, 
/// and the property data for each element
namespace MeshDS
{
  class FEMesh;
  class MeshList;

  class FEMeshData
  {
public:
    /// constructor
    FEMeshData(const MeshDS::FEMesh& );
    
    
    /// destructor
    ~FEMeshData();
	  
    /// @returns the node whose foreign ID is given in the input parameter
    inline const Node* getNodeFromForeignID (const unsigned int fid) const;
    
    /// @returns the elem whose foreign ID is given in the input parameter
    inline const Elem* getElemFromForeignID (const unsigned int fid) const;
    
    /// @returns the foreign ID of the element
    inline unsigned int getForeignIDFromElem (const Elem* elem) const;
    
    /// @returns the foreign ID of the node
    inline unsigned int getForeignIDFromNode (const Node* elem) const;
    
    /// @returns the property ID of the element
    inline unsigned int getElemPropertyID(const Elem* elem) const;
    
    /// @returns the elem kind enum ID
    inline unsigned int getElemKindEnumID(const Elem* elem) const;
    
    /// @returns the elem kind enum name
    inline std::string getElemKindEnumName(const Elem* elem) const;
    
    /// @returns the node whose internal ID is given in the input parameter
    inline const Node* getNodeFromInternalID (const unsigned int ) const;
    
    /// @returns the internal ID of the given node
    inline unsigned int getInternalIDFromNode (const Node* ) const;

    /// @returns the foreign ID of a node whose internal ID is given
    inline unsigned int getNodeForeignIDFromInternalID (const unsigned int ) const;
    
    /// writes the mesh data to the provided output stream. The node ID, location, and the 
    /// global dof IDs are written for each node
    void writeMeshDataToStream(std::ostream& output);

protected:
      
      
      /// this method adds the node and its foreign ID to the mesh data
      void addNodeAndIDData (const Node* node, 
                             const unsigned int foreign_node_id);
    
    /// this method adds the elem and its foreign ID to the mesh data
    void addElemAndIDData (const Elem* elem, 
                           const unsigned int foreign_elem_ID,
                           const unsigned int property_ID, 
                           const unsigned int elem_kind_enum_ID);
    
    /// @returns the element whose internal ID is given in the input parameter
    inline const Elem* getElemFromInternalID (const unsigned int ) const;
    
    /// @returns the internal ID of the given element 
    inline unsigned int getInternalIDFromElem (const Elem* ) const;
    
    
    /// @returns the internal ID of a node whose foreign ID is given 
    inline unsigned int getNodeInternalIDFromForeignID (const unsigned int ) const;
    
    /// @returns the internal ID of an elem whose foreign ID is given 
    inline unsigned int getElemInternalIDFromForeignID (const unsigned int ) const;
    
    
    /// @returns the foreign ID of an elem whose internal ID is given
    inline unsigned int getElemForeignIDFromInternalID (const unsigned int ) const;
    
    /// this exception is thrown if a duplicate ID is specified
    DeclException2(ExcDuplicateID, std::string, unsigned int,
                   << "Duplicate " << arg1 << " ID: " << arg2 << " specified." );
    
    /// this exception is thrown if a duplicate ID is specified
    DeclException1(ExcDuplicatePointer, std::string, 
                   << "Duplicate " << arg1 << " pointer specified." );
    
    /// this exception is thrown if the ID does not exist in the data
    DeclException2(ExcIDDoesNotExist, std::string, unsigned int,
                   << "Specified " << arg1 << " ID: " << arg2 << " does not exist." );
    
    /// this exception is thrown if the pointer does not exist in the data
    DeclException1(ExcPointerDoesNotExist, std::string,
                   << "Specified " << arg1 << " pointer does not exist." );
    
    /// this exception is thrown if the internal ID does not exist in the data
    DeclException2(ExcInternalIDDoesNotExist, std::string, unsigned int,
                   << "Specified " << arg1 << "  internal ID: " << arg2 << " does not exist." );
    
    /// the mesh whose data is stored in this mesh data
    const MeshDS::FEMesh& mesh;
    
    
    /// a structure for storing the element property and kind enum ID
    struct ElemPropIDs
      {
        unsigned int foreign_ID;
        unsigned int internal_ID;
        unsigned int property_ID;
        unsigned int kind_enum_ID;
      };
    
    /// a structure for storing the node IDs
    struct NodeIDs
      {
        unsigned int foreign_ID;
        unsigned int internal_ID;
      };
    
    
    // local type definitions
    typedef std::map<unsigned int, const Elem*> IDToElemMap;
    typedef std::map<unsigned int, const Node*> IDToNodeMap;
    typedef std::map< const Elem*, ElemPropIDs> ElemToIDMap;
    typedef std::map< const Node*, NodeIDs> NodeToIDMap;
    typedef std::vector<const Elem*> ElemVector;
    typedef std::vector<const Node*> NodeVector;
    
    /// map for storing elems vs their external IDs
    IDToElemMap foreign_ID_to_elem_map;
    
    /// map for storing nodes vs their external IDs
    IDToNodeMap foreign_ID_to_node_map;
    
    /// map for storing elem IDs vs their pointers
    ElemToIDMap elem_ID_map;
    
    /// map for storing node IDs vs their pointers
    NodeToIDMap node_ID_map;
    
    /// vector of elem, where the location of the elem in the vector 
    /// is its internal ID
    ElemVector elem_vector;
    
    /// vector of nodes, where the location of the node in the vector 
    /// is its internal ID
    NodeVector node_vector;  
    
    
    /// the MeshList class is a friend since it performs the input
    friend class MeshList;
    
  };
}



inline 
const Node* 
MeshDS::FEMeshData::getNodeFromForeignID (const unsigned int fid) const
{
  MeshDS::FEMeshData::IDToNodeMap::const_iterator it, end;
  it = this->foreign_ID_to_node_map.find(fid);
  end = this->foreign_ID_to_node_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcIDDoesNotExist("Node", fid));
  
  return it->second;
}



inline
unsigned int
MeshDS::FEMeshData::getForeignIDFromNode (const Node* node) const
{
  MeshDS::FEMeshData::NodeToIDMap::const_iterator it, end;
  it = this->node_ID_map.find(node);
  end = this->node_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Node"));
  
  return it->second.foreign_ID;
}




inline
const Elem* 
MeshDS::FEMeshData::getElemFromForeignID (const unsigned int fid) const
{
  MeshDS::FEMeshData::IDToElemMap::const_iterator it, end;
  it = this->foreign_ID_to_elem_map.find(fid);
  end = this->foreign_ID_to_elem_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcIDDoesNotExist("Elem", fid));
  
  return it->second;
}




inline 
unsigned int 
MeshDS::FEMeshData::getForeignIDFromElem (const Elem* elem) const
{
  MeshDS::FEMeshData::ElemToIDMap::const_iterator it, end;
  it = this->elem_ID_map.find(elem);
  end = this->elem_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Elem"));
  
  return it->second.foreign_ID;
}





inline
const Node* 
MeshDS::FEMeshData::getNodeFromInternalID (const unsigned int internal_ID) const
{
  Assert(internal_ID < this->node_vector.size(),
         MeshDS::FEMeshData::ExcInternalIDDoesNotExist("Node", internal_ID));
  return this->node_vector[internal_ID];
}





inline 
const Elem*
MeshDS::FEMeshData::getElemFromInternalID (const unsigned int internal_ID) const
{
  Assert(internal_ID < this->elem_vector.size(),
         MeshDS::FEMeshData::ExcInternalIDDoesNotExist("Elem", internal_ID));
  return this->elem_vector[internal_ID];
}




inline 
unsigned int
MeshDS::FEMeshData::getInternalIDFromElem (const Elem* elem) const
{
  MeshDS::FEMeshData::ElemToIDMap::const_iterator it, end;
  it = this->elem_ID_map.find(elem);
  end = this->elem_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Elem"));
  
  return it->second.internal_ID;
}




inline 
unsigned int
MeshDS::FEMeshData::getInternalIDFromNode (const Node* node) const
{
  MeshDS::FEMeshData::NodeToIDMap::const_iterator it, end;
  it = this->node_ID_map.find(node);
  end = this->node_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Node"));
  
  return it->second.internal_ID;
}





inline 
unsigned int
MeshDS::FEMeshData::getNodeInternalIDFromForeignID (const unsigned int fid) const
{  
	const Node* node = this->getNodeFromForeignID(fid);
  
  return this->getInternalIDFromNode(node);
}




inline
unsigned int
MeshDS::FEMeshData::getElemInternalIDFromForeignID (const unsigned int fid) const
{
	const Elem* elem = this->getElemFromForeignID(fid);
  
  return this->getInternalIDFromElem(elem);
}






inline
unsigned int
MeshDS::FEMeshData::getNodeForeignIDFromInternalID (const unsigned int internal_ID) const
{
  const Node* node = this->getNodeFromInternalID(internal_ID);
  return this->getForeignIDFromNode(node);
}






inline
unsigned int
MeshDS::FEMeshData::getElemForeignIDFromInternalID (const unsigned int internal_ID) const
{
  const Elem* elem = this->getElemFromInternalID(internal_ID);
  return this->getForeignIDFromElem(elem);
}







inline 
unsigned int
MeshDS::FEMeshData::getElemPropertyID(const Elem* elem) const
{
  MeshDS::FEMeshData::ElemToIDMap::const_iterator it, end;
  it = this->elem_ID_map.find(elem);
  end = this->elem_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Elem"));
  
  return it->second.property_ID;
}



inline
unsigned int 
MeshDS::FEMeshData::getElemKindEnumID(const Elem* elem) const
{
  MeshDS::FEMeshData::ElemToIDMap::const_iterator it, end;
  it = this->elem_ID_map.find(elem);
  end = this->elem_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Elem"));
  
  return it->second.kind_enum_ID;
}



inline
std::string
MeshDS::FEMeshData::getElemKindEnumName(const Elem* elem) const
{
  MeshDS::FEMeshData::ElemToIDMap::const_iterator it, end;
  it = this->elem_ID_map.find(elem);
  end = this->elem_ID_map.end();
	
  Assert(it != end,
         MeshDS::FEMeshData::ExcPointerDoesNotExist("Elem"));
  
  return FESystemElem::FESystemElemTypeEnum::enumName(it->second.property_ID);
}



#endif // __fesystem_fe_mesh_data_h__
