// $Id: MeshList.h,v 1.9.6.1 2008-02-25 04:27:01 manav Exp $

#ifndef __fesystem_mesh_list_h__
#define __fesystem_mesh_list_h__

// C++ includes
#include <memory>
#include <iostream>
#include <map>
#include <set>

// FESystem includes
#include "FESystem/FESystemConfig.h"
#include "FESystem/FESystemExceptions.h"

// libMesh includes
#include "mesh/mesh.h"
#include "base/dof_map.h"

// Forward declerations
class RadiationCavity;
class Mesh;
class DofMap;
class Elem;


namespace MeshDS
{
  
  class FEMesh: public Mesh
  {
public: 
    FEMesh(unsigned int dimension): Mesh(dimension){}
    ~FEMesh(){}
  };

  class FEDofMap: public DofMap
    {
public:
      FEDofMap(): DofMap(0){}
      ~FEDofMap(){}
    };

  class FEMeshData;
  
class MeshList
{
 public:
  /// constructor
  MeshList();
	
  /// destructor
  ~MeshList();
	  
  /// returns mesh with the given ID
  MeshDS::FEMesh* getMeshFromID(const unsigned int );

  /// returns mesh data for the given ID
  MeshDS::FEMeshData* getMeshDataFromID(const unsigned int );

  /// returns the DofMap associated with this ID
  MeshDS::FEDofMap* getDofMapFromID(const unsigned int );

  /// returns the radiation cavity associated with the ID
  RadiationCavity* getRadiationCavity(const unsigned int );

  /// returns the dummy elem set associated with ID
  std::set<Elem*>* getDummyElemSetFromID(const unsigned int mesh_ID);

  /// method to read the mesh data from input stream
  std::istream& readFromInputStream(std::istream& );
  
 protected:
	
  /// this method clears the data structures
  void clear();
  
  /// adds new mesh and related data structures to this class
  void addNewMeshDataStructures(const unsigned int mesh_ID,
                                const unsigned int dimension);

  /// adds a radiation cavity to this class
  void addRadiationCavity(RadiationCavity* );
	
  // local type definitions
  typedef std::pair<MeshDS::FEMesh*, MeshDS::FEMeshData*> MeshAndMeshDataPair;
  typedef std::map<unsigned int, MeshAndMeshDataPair> IDToMeshPairMap;
  typedef std::map<unsigned int, MeshDS::FEDofMap*> IDToDofMapMap;
  typedef std::map<unsigned int, RadiationCavity*> IDToRadiationCavityMap;
  typedef std::map<unsigned int, std::set<Elem*>*> IDToElemSetMap;
  
  /// Mesh and MeshData maps for the MeshList object
  IDToMeshPairMap meshpair_ID_map;
	
  /// map of mesh ID to DofMap  
  IDToDofMapMap mesh_ID_to_dof_map;

  /// map of ID to  Radiation Cavity
  IDToRadiationCavityMap radiation_cavity_ID_map;

  IDToElemSetMap dummy_elem_set_map;
};
}


inline 
MeshDS::FEMesh* 
MeshDS::MeshList::getMeshFromID(const unsigned int ID) 
{
  // get the pair, and make sure the ID existed
  MeshDS::MeshList::IDToMeshPairMap::const_iterator map_it = this->meshpair_ID_map.find(ID);
	
  Assert(map_it != this->meshpair_ID_map.end(),
         FESystemExceptions::ExcIDDoesNotExist("Mesh",ID));
	
  return map_it->second.first;
}




inline
MeshDS::FEMeshData* 
MeshDS::MeshList::getMeshDataFromID(const unsigned int ID) 
{
  // get the pair, and make sure the ID existed
  MeshDS::MeshList::IDToMeshPairMap::const_iterator map_it = this->meshpair_ID_map.find(ID);
	
  Assert(map_it != this->meshpair_ID_map.end(),
         FESystemExceptions::ExcIDDoesNotExist("Mesh",ID));
	
  return map_it->second.second;
}



inline
RadiationCavity* 
MeshDS::MeshList::getRadiationCavity(const unsigned int ID) 
{
  // get the pair, and make sure the ID existed
  MeshDS::MeshList::IDToRadiationCavityMap::const_iterator map_it = this->radiation_cavity_ID_map.find(ID);
	
  Assert(map_it != this->radiation_cavity_ID_map.end(),
         FESystemExceptions::ExcIDDoesNotExist("RadiationCavity",ID));
	
  return map_it->second;
}



inline
MeshDS::FEDofMap* 
MeshDS::MeshList::getDofMapFromID(const unsigned int mesh_id)
{
  MeshDS::MeshList::IDToDofMapMap::const_iterator map_it = 
  this->mesh_ID_to_dof_map.find(mesh_id);
  
  Assert(map_it != this->mesh_ID_to_dof_map.end(),
         FESystemExceptions::ExcIDDoesNotExist("DofMap", mesh_id));
  
  return map_it->second;
}




inline
std::set<Elem*>* 
MeshDS::MeshList::getDummyElemSetFromID(const unsigned int mesh_ID)
{
  // if the mesh ID does not exist in the map, return a NULL
  MeshDS::MeshList::IDToElemSetMap::const_iterator 
  it = this->dummy_elem_set_map.find(mesh_ID);
  
  Assert(it != this->dummy_elem_set_map.end(),
         FESystemExceptions::ExcIDDoesNotExist("Dummy Elem Set", mesh_ID));

  return it->second;
}



#endif // __fesystem_mesh_list_h__
