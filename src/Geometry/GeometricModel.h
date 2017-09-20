// $Id: GeometricModel.h,v 1.3 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_model_h__
#define __geometric_model_h__

// C++ includes
#include <map>
#include <iostream>

// FESystem includes
#include "Geometry/GeometricPoint.h"
#include "Geometry/GeometricStraightLine.h"
#include "Geometry/GeometricCircularArc.h"
#include "Geometry/GeometricEllipticArc.h"
#include "Geometry/GeometricSurface.h"
#include "Geometry/GeometricSet.h"
#include "Mesh/FEMeshData.h"
#include "FESystem/FESystemConfig.h"
#include "Loads/load.h"

// Libmesh includes
#include "mesh/mesh.h"
#include "geom/elem.h"



class GeometricModel
{
 public:
  
  /// constructor 
  GeometricModel();

  /// destructor 
  ~GeometricModel();

  /// returns a point corresponding to the unique ID
  inline GeometricPoint& getPointForID(unsigned int ) const;

  /// returns a line corresponding to the unique ID
  inline GeometricLine& getLineForID(unsigned int ) const;
  
  /// returns a surface corresponding to the unique ID
  inline GeometricSurface& getSurfaceForID(unsigned int ) const;

  /// add geometry set to this model
  GeometricSet& getGeometricSet(const std::string& );

  /// reads the input configuration file
  void readModelConfigurationInputFile(std::istream&);

  /// reads from an input stream
  void readFromInput(std::istream& );

  /// reads the Gmsh mesh output
  void readGmshMeshOutput(std::istream& );

  /// writes FESystem input file
  void writeFESystemInputFile(std::ostream& );

  /// writes nastran input file
  void writeNastranInputFile(std::ostream& );

 protected:
  
  /// method adds geometric entity to a geometric set 
  template<typename EntityType >
  void addGeometricEntityToSet(std::string& , EntityType * );


  /// method to create sets for tags and add the entity to the set
  /// @param input stream from which to read
  /// @param geometric entity to be added to the set
  template<typename EntityType>
    void readTagsAndCreateSets(std::istream& , EntityType* );


  /// prepares a list of writable elements for the model
  void createWritableElemList();

  /// marks all the duplicates in the model
  void markDuplicateEntities();

  /// method to add geometric point to model
  void addGeometricEntityToModel(GeometricPoint* );

  /// method to add geometric Line to model
  void addGeometricEntityToModel(GeometricLine* );

  /// method to add geometric surface to model
  void addGeometricEntityToModel(GeometricSurface* );

  /// method to add geometric point to tag
  void addGeometricEntityToTag(GeometricPoint* );

  /// method to add geometric Line to tag
  void addGeometricEntityToTag(GeometricLine* );

  /// method to add geometric surface to tag
  void addGeometricEntityToTag(GeometricSurface* );

  /// adds a tag to the model
  void addGeometryTag(std::string&);

  /// map of points. Here, the points are stored against their unique IDs
  std::map<unsigned int, GeometricPoint*> point_map;

  /// map of lines against unique IDs
  std::map<unsigned int, GeometricLine*> line_map;

  /// map of surfaces, against unique IDs
  std::map<unsigned int, GeometricSurface*> surface_map;

  /// map of geometry set
  std::map<std::string, GeometricSet*> set_map;

  /// mesh and mesh data for the finite element mesh
  MeshDS::FEMesh* mesh;
  MeshDS::FEMeshData* mesh_data;

  /// name of the geometry input file
  std::string geometry_file;

  /// name of the mesh input file
  std::string mesh_file;

  /// map of line entity ID and the associated 1-D elements
  std::multimap<unsigned int, Elem*> line_entity_ID_elem_map;

  /// map of surface entity ID and the associated 2-D elements
  std::multimap<unsigned int, Elem*> surface_entity_ID_elem_map;

  /// vector of sets that will be written to the FE model
  std::vector<std::string> writable_set_vector;

  /// map of sets and their material IDs
  std::map<std::string , unsigned int> set_material_ID_map;

  /// map of sets and boundary conditions on them
  std::map<std::string, Load> set_boundary_condition_map;

  /// map of ID and elem which will be written in the file
  std::map<unsigned int, Elem*> writable_elem_map;

  /// map of ID and nodes which will be written
  std::map<unsigned int, Node*> writable_node_map;
};


inline
GeometricPoint& GeometricModel::getPointForID(unsigned int ID) const
{
  std::map<unsigned int, GeometricPoint*>::const_iterator it;
  it = this->point_map.find(ID);
  assert (it != this->point_map.end());

  return *(it->second);
}


inline
GeometricLine& GeometricModel::getLineForID(unsigned int ID) const
{
  std::map<unsigned int, GeometricLine*>::const_iterator it;
  it = this->line_map.find(ID);
  assert (it != this->line_map.end());
  
  return *(it->second);
}


inline
GeometricSurface& GeometricModel::getSurfaceForID(unsigned int ID) const
{
  std::map<unsigned int, GeometricSurface*>::const_iterator it;
  it = this->surface_map.find(ID);
  assert (it != this->surface_map.end());

  return *(it->second);  
}



template<typename EntityType >
void GeometricModel::addGeometricEntityToSet(std::string& name, EntityType * entity)
{
  // get the set by this name
  GeometricSet& set = this->getGeometricSet(name);

  // insert this point in it
  set.addEntity(entity);
}




template<typename EntityType>
void GeometricModel::readTagsAndCreateSets(std::istream& input, EntityType* entity)
{
  std::string tag;
  unsigned int n_tags=0;

  input >> tag;
  assert (tag == "TAGS");
  
  input >> n_tags;
  
  for (unsigned int i=0; i< n_tags; i++)
    {
      tag.clear();
      input >> tag;
      
      this->addGeometricEntityToSet(tag, entity);
    }
  
}



#endif // __geometric_model_h__
