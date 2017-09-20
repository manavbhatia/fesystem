// $Id: GeometricModel.C,v 1.3 2006-09-05 20:41:56 manav Exp $

// C++ includes


// FESystem includes
#include "Geometry/GeometricModel.h"
#include "Geometry/GeometricSet.h"
#include "Mesh/FEMeshData.h"
#include "Mesh/MeshList.h"


// libMesh includes
#include "geom/edge_edge2.h"
#include "geom/face_quad4.h"
#include "geom/face_tri3.h"
#include "geom/cell_hex8.h"


GeometricModel::GeometricModel():
  mesh(NULL),
  mesh_data(NULL)
{
    this->mesh = new MeshDS::FEMesh();
    this->mesh_data = new MeshDS::FEMeshData(*(this->mesh));

}


GeometricModel::~GeometricModel()
{
  // iterate over all the sets, and delete them
  // iterate over all the surfaces and delte them
  // then the lines, 
  // then the points

  // delete the mesh and mesh data

  
}


void GeometricModel::readFromInput(std::istream& input)
{
  std::string tag, set_name;
  // Look for the first begin_geometry tag. The entire data should be between
  // the first begin_geometry and last end_geometry tag

  while (tag != "BEGIN_GEOMETRY_DATA")
    input >> tag; 
  tag.clear();

  while (tag != "END_GEOMETRY_DATA")
    {
      input.clear(); input >> tag;
      assert (tag == "BEGIN_GEOMETRY");
      
      // read the geometry section till end_geometry is not encountered
      while (tag != "END_GEOMETRY")
	{
	  set_name.clear();
	  input >> set_name; // this is the name of the section. All other entities will be a part of this
	  assert (tag != "END_GEOMETRY"); // it should contain the name of the entity
	  
	  tag.clear();
	  input >> tag;
	  if (tag == "POINT")
	    {
	      // create a point and read it in
	      GeometricPoint* point = new GeometricPoint(*this);
	      input >> (*point);

	      // add it to the model
	      this->addGeometricEntityToModel(point);

	      // also, add it to the high level
	      this->addGeometricEntityToSet(set_name, point);

	      // read the tags
	      this->readTagsAndCreateSets(input, point);
	    }
	  else if (tag == "STRAIGHT_LINE")
	    {
	      // create a point and read it in
	      GeometricStraightLine* line = new GeometricStraightLine(*this);
	      input >> (*line);

	      // add it to the model
	      this->addGeometricEntityToModel(line);

	      // also, add it to the high level
	      this->addGeometricEntityToSet(set_name, line);

	      // read the tags
	      this->readTagsAndCreateSets(input, line);
	    }
	  else if (tag == "CIRCULAR_ARC")
	    {
	      // create a point and read it in
	      GeometricCircularArc* line = new GeometricCircularArc(*this);
	      input >> (*line);

	      // add it to the model
	      this->addGeometricEntityToModel(line);

	      // also, add it to the high level
	      this->addGeometricEntityToSet(set_name, line);

	      // read the tags
	      this->readTagsAndCreateSets(input, line);
	    }
	  else if (tag == "ELLIPTIC_ARC")
	    {
	      // create a point and read it in
	      GeometricEllipticArc* line = new GeometricEllipticArc(*this);
	      input >> (*line);

	      // add it to the model
	      this->addGeometricEntityToModel(line);

	      // also, add it to the high level
	      this->addGeometricEntityToSet(set_name, line);

	      // read the tags
	      this->readTagsAndCreateSets(input, line);
	    }
	  else if (tag == "SURFACE")
	    {
	      // create a point and read it in
	      GeometricSurface* surface = new GeometricSurface(*this);
	      input >> (*surface);

	      // add it to the model
	      this->addGeometricEntityToModel(surface);

	      // also, add it to the high level
	      this->addGeometricEntityToSet(set_name, surface);

	      // read the tags
	      this->readTagsAndCreateSets(input, surface);
	    }
	  
	}
      
    }

  // after having red in the geometric model, remove all the duplicate entities
  this->markDuplicateEntities();
}



void GeometricModel::markDuplicateEntities()
{
  // *************** first work on all the points *************

  // keep two iterators to the points, one for keeping track of unique IDs, and 
  // the other for duplicates
  std::map<unsigned int , GeometricPoint*>::iterator unique_point_it, duplicate_point_it;
  std::map<unsigned int, GeometricPoint*>::const_iterator point_end;
  
  unique_point_it = this->point_map.begin();
  point_end = this->point_map.end();

  for (; unique_point_it != point_end;  unique_point_it++)
    {
      // if the unique point ID has not been marked as a duplicate, 
      // then iterate over all other points and mark the duplicates

      if (unique_point_it->second->isDuplicate() == false)
	{
	  // set the duplicate point iterator as one after the unique point it
	  duplicate_point_it = unique_point_it;
	  duplicate_point_it++;

	  // next, iterate over all the points in the dupicate iterator, and mark the
	  // dupliates
	  for (; duplicate_point_it != point_end; duplicate_point_it++)
	    {
	      // if the point has not been marked as duplicate, compare, else do nothing
	      if (duplicate_point_it->second->isDuplicate() != false)
		{
		  duplicate_point_it->second->checkIfDuplicate( (*unique_point_it->second));
		}
	    }
	}
      
    }

  // ************ now the lines ***************

  // keep two iterators to the lines, one for keeping track of unique IDs, and 
  // the other for duplicates
  std::map<unsigned int , GeometricLine*>::iterator unique_line_it, duplicate_line_it;
  std::map<unsigned int, GeometricLine*>::const_iterator line_end;
  
  unique_line_it = this->line_map.begin();
  line_end = this->line_map.end();

  for (; unique_line_it != line_end;  unique_line_it++)
    {
      // if the unique line ID has not been marked as a duplicate, 
      // then iterate over all other lines and mark the duplicates

      if (unique_line_it->second->isDuplicate() == false)
	{
	  // set the duplicate line iterator as one after the unique line it
	  duplicate_line_it = unique_line_it;
	  duplicate_line_it++;

	  // next, iterate over all the lines in the dupicate iterator, and mark the
	  // dupliates
	  for (; duplicate_line_it != line_end; duplicate_line_it++)
	    {
	      // if the line has not been marked as duplicate, compare, else do nothing
	      if (duplicate_line_it->second->isDuplicate() != false)
		{
		  duplicate_line_it->second->checkIfDuplicate( (*unique_line_it->second));
		}
	    }
	}
      
    }


  // ************* now the surfaces  ************************

  // keep two iterators to the surfaces, one for keeping track of unique IDs, and 
  // the other for duplicates
  std::map<unsigned int , GeometricSurface*>::iterator unique_surface_it, duplicate_surface_it;
  std::map<unsigned int, GeometricSurface*>::const_iterator surface_end;
  
  unique_surface_it = this->surface_map.begin();
  surface_end = this->surface_map.end();

  for (; unique_surface_it != surface_end;  unique_surface_it++)
    {
      // if the unique surface ID has not been marked as a duplicate, 
      // then iterate over all other surfaces and mark the duplicates

      if (unique_surface_it->second->isDuplicate() == false)
	{
	  // set the duplicate surface iterator as one after the unique surface it
	  duplicate_surface_it = unique_surface_it;
	  duplicate_surface_it++;

	  // next, iterate over all the surfaces in the dupicate iterator, and mark the
	  // dupliates
	  for (; duplicate_surface_it != surface_end; duplicate_surface_it++)
	    {
	      // if the surface has not been marked as duplicate, compare, else do nothing
	      if (duplicate_surface_it->second->isDuplicate() != false)
		{
		  duplicate_surface_it->second->checkIfDuplicate( (*unique_surface_it->second));
		}
	    }
	}
      
    }


}



void GeometricModel::addGeometricEntityToModel(GeometricPoint* point)
{
  // get the ID of this entity
  unsigned int ID = point->ID();

  // make sure that this ID does not already exist in the map
  assert (this->point_map.find(ID) == this->point_map.end());

  // add the entity to the map
  bool insert = this->point_map.insert(std::map<unsigned int, GeometricPoint*>::value_type(ID, point)).second;
  assert (insert == true);
}



void GeometricModel::addGeometricEntityToModel(GeometricLine* line)
{
  // get the ID of this entity
  unsigned int ID = line->ID();

  // make sure that this ID does not already exist in the map
  assert (this->line_map.find(ID) == this->line_map.end());

  // add the entity to the map
  bool insert = this->line_map.insert(std::map<unsigned int, GeometricLine*>::value_type(ID, line)).second;
  assert (insert == true);

}




void GeometricModel::addGeometricEntityToModel(GeometricSurface* surface)
{
  // get the ID of this entity
  unsigned int ID = surface->ID();

  // make sure that this ID does not already exist in the map
  assert (this->surface_map.find(ID) == this->surface_map.end());

  // add the entity to the map
  bool insert = this->surface_map.insert(std::map<unsigned int, GeometricSurface*>::value_type(ID, surface)).second;
  assert (insert == true);
}




GeometricSet& GeometricModel::getGeometricSet(const std::string& name)
{
  // get the iterator for this name
  std::map<std::string, GeometricSet*>::iterator set_it = this->set_map.find(name);  

  // if the set does not exist, create it
  if (set_it == this->set_map.end())
    {
      // create a new set 
      GeometricSet* set = new GeometricSet(name);

      // insert the set
      std::pair<std::map<std::string, GeometricSet*>::iterator, bool > insert_pair = 
	this->set_map.insert(std::map<std::string, GeometricSet*>::value_type(name, set));
      assert (insert_pair.second == true);
      
      set_it = insert_pair.first;
    }

  return *(set_it->second);
}





void GeometricModel::readGmshMeshOutput(std::istream& input)
{
  std::string tag;
  unsigned int n_nodes=0, n_elem=0;

  // read the nodes
  input >> tag;
  assert(tag == "$NOD");

  input >> n_nodes;

  //loop over all the nodes and read the data
  {
    unsigned int node_id = 0; double x_loc = 0.0, y_loc = 0.0, z_loc = 0.0;
    for (unsigned int i = 0; i < n_nodes; i++)
      {
	//read the nodal data
	input >> node_id;	// node ID
	input >> x_loc;	// x location of the node
	input >> y_loc;	// y location
	input >> z_loc;	// z location
				
	//create the point and add it to the mesh for this FEsystem
	Point this_point(x_loc, y_loc, z_loc);
//	Node* this_node = this->mesh->add_point(this_point);
				
	//use the ID of the node and create an entry in mesh_data
//	this->mesh_data->addNodeAndIDData(this_node, node_id);
      }
  }
  
  tag.clear(); 
  input >> tag;
  assert (tag == "$ENDNOD");


  //******************	read the element data******************
  //
  tag.clear();
  input >> tag;
  assert (tag == "$ELM");
			
  input >> n_elem;
		
  //loop over all the elements and read the data
  {
//    const unsigned int DUMMY_NUM = 1;
    unsigned int elem_id, elem_type, reg_phys, reg_elem, n_nodes, node_id;
    Elem* this_elem;
    for (unsigned int i = 0; i < n_elem; i++)
      {
	//read the element data
	input >> elem_id;	// element ID
	input >> elem_type;
				//use the element elem_type to appropriately read the file
	if (elem_type == 1)
	  {
	    //	2 noded Edge element
	    this_elem = mesh->add_elem(new Edge2);
	  }
	else if (elem_type == 2)
	  {
	    //	3 noded Tri element
	    this_elem = mesh->add_elem(new Tri3);
	  }
	else if (elem_type == 3)
	  {
	    //	4 noded Quad element
	    this_elem = mesh->add_elem(new Quad4);
	  }
	else if (elem_type == 5)
	  {
	    //	8 noded hex element
	    this_elem = mesh->add_elem(new Hex8);
	  }
	else 
	  {
	    std::cout << "Error!! Invalid Element type specified." << std::endl;
	    abort();
	  }
				
	
	//now, repeat and read the ID, property ID and connectivity
	input >> reg_phys;	// registered physical element
	input >> reg_elem;	// registered geometric entity
	input >> n_nodes;       // number of nodes for this elem

	// make sure that the elems match
	assert (n_nodes == this_elem->n_nodes());

	// the current status of the code cannot handle physical elements. Hence, the reg_elem and reg_phys
	// should be the same numbers
	assert (reg_phys == reg_elem);
	
	//read the node IDs and set the element connectivity
	for (unsigned int j=0; j < n_nodes; j++)
	  {
	    input >> node_id;
	    this_elem->set_node(j) = 
	      const_cast<Node*>(this->mesh_data->getNodeFromForeignID(node_id));
	  }
	
	//next, insert the element and its ID into mesh data
//	this->mesh_data->addElemAndIDData(this_elem, elem_id);
//	this->mesh_data->addElemPropertyIDAndTag(this_elem, DUMMY_NUM,"DUMMY_TAG");
	
	// also insert this elem into the entity ID to elem map
	switch (this_elem->type())
	  {
	  case EDGE2:
	    {
	      this->line_entity_ID_elem_map.insert(std::multimap<unsigned int, Elem*>::value_type(reg_elem, this_elem));
	    }
	    break;

	  case TRI3:
	  case QUAD4:
	    {
	      this->surface_entity_ID_elem_map.insert(std::multimap<unsigned int, Elem*>::value_type(reg_elem, this_elem));
	    }
	    break;

	  default:
	    abort();
	    break;
	  }

      }
    //now clear the tag
    tag.clear();
    input >> tag;
    assert (tag == "$ENDELM");

  }
}


void GeometricModel::createWritableElemList()
{
  // get the list of sets for which the elements will have to be written
  std::vector<std::string>::const_iterator set_it, set_end;
  set_it = this->writable_set_vector.begin();
  set_end = this->writable_set_vector.end();
 
  std::map<unsigned int, GeometricLine*> unique_lines;
  std::map<unsigned int, GeometricSurface*> unique_surfaces;

  for ( ; set_it != set_end; set_it++)
    {
      GeometricSet& set = this->getGeometricSet(*set_it);

      // get the geometric entities from these sets
      const std::vector<GeometricLine*>& lines = set.getLines();
      const std::vector<GeometricSurface*>& surfaces = set.getSurfaces();

      // get the iterator of the lines
      std::vector<GeometricLine*>::const_iterator line_it, line_end;
      std::vector<GeometricSurface*>::const_iterator surface_it, surface_end;

      line_it = lines.begin(); line_end = lines.end();
      surface_it = surfaces.begin(); surface_end = surfaces.end();
      unsigned int line_id = 0, surface_id = 0;
      
      for ( ; line_it != line_end; line_it++ )
	{
	  // if this line ID does not exist in the unique line map, insert it in
	  line_id = (*line_it)->baseEntityID();
	  
	  std::map<unsigned int,GeometricLine*>::const_iterator it = unique_lines.find(line_id);
	  if (it == unique_lines.end())
	    {
	      // insert the ID
	      bool insert = unique_lines.insert(std::map<unsigned int, GeometricLine*>::value_type(line_id,*line_it)).second;
	      assert (insert == true);
	    }
	} 
      
      // repeat for the surfaces
      for ( ; surface_it != surface_end; surface_it++ )
	{
	  // if this surface ID does not exist in the unique surface map, insert it in
	  surface_id = (*surface_it)->baseEntityID();
	  
	  std::map<unsigned int,GeometricSurface*>::const_iterator it = unique_surfaces.find(surface_id);
	  if (it == unique_surfaces.end())
	    {
	      // insert the ID
	      bool insert = unique_surfaces.insert(std::map<unsigned int, GeometricSurface*>::value_type(surface_id,*surface_it)).second;
	      assert (insert == true);
	    }
	} 
    }


  // look for the elements that belong to these unique IDs
  // it is assumed that 1-D elements will belong to 1-D entities, and so on for 2-D and 3-D
  // for this the map of 1-D entities and nodes will be used. 
  // 
  // First process the line elements
  // iterate over the line entities, for each line entity get the ID, and in the map of 
  // line-entity ID vs elems, get the vector of elems that belong to this entity.
  // add those elements in the vector of writable elements

  unsigned int elem_id = 0;

  std::map<unsigned int, GeometricLine*>::const_iterator line_it, line_end;
  line_it = unique_lines.begin();
  line_end = unique_lines.end();
  
  std::multimap<unsigned int, Elem*>::const_iterator elem_it, elem_end;
  std::pair<std::multimap<unsigned int, Elem*>::const_iterator, 
    std::multimap<unsigned int, Elem*>::const_iterator> it_pair;

  for (; line_it != line_end; line_it++)
    {
      it_pair = this->line_entity_ID_elem_map.equal_range(line_it->first);
      elem_it = it_pair.first;
      elem_end = it_pair.second;
      
      // make sure that some elements exist for this entity
      assert (elem_it != elem_end);
      assert (elem_end != this->line_entity_ID_elem_map.end());

      for (; elem_it != elem_end; elem_end++)
	{
	  // get the ID of this element
	  elem_id = this->mesh_data->getForeignIDFromElem(elem_it->second);

	  // insert this element to the writable element map
	  bool insert = 
	    this->writable_elem_map.insert(std::map<unsigned int, Elem*>::value_type(elem_id, elem_it->second)).second;
	  assert (insert == true);
	}      
    }
  
  // Following this, process the surface elements
  std::map<unsigned int, GeometricSurface*>::const_iterator surf_it, surf_end;
  surf_it = unique_surfaces.begin();
  surf_end = unique_surfaces.end();

  for (; surf_it != surf_end; surf_it++)
    {
      it_pair = this->surface_entity_ID_elem_map.equal_range(surf_it->first);
      elem_it = it_pair.first;
      elem_end = it_pair.second;

      // make sure that some elements exist for this entity
      assert (elem_it != elem_end);
      assert (elem_end != this->line_entity_ID_elem_map.end());

      for (; elem_it != elem_end; elem_end++)
	{
	  // get the ID of this element
	  elem_id = this->mesh_data->getForeignIDFromElem(elem_it->second);

	  // insert this element to the writable element map
	  bool insert = 
	    this->writable_elem_map.insert(std::map<unsigned int, Elem*>::value_type(elem_id, elem_it->second)).second;
	  assert (insert == true);
	}
    }

  // create the list of nodes that need to be written
  // for iterate over all elements in the writable element map, find their node, insert them in the 
  // writable elem map it does not already exist
  std::map<unsigned int, Elem*>::const_iterator model_elem_it, model_elem_end;
  model_elem_it = this->writable_elem_map.begin();
  model_elem_end = this->writable_elem_map.end();

  unsigned int n_nodes = 0, node_id = 0;
  Node* node= NULL;

  for (; model_elem_it != model_elem_end; model_elem_it++)
    {
      n_nodes = model_elem_it->second->n_nodes();
      for (unsigned int i=0; i < n_nodes; i++)
	{
	  node = model_elem_it->second->get_node(i);
	  node_id = this->mesh_data->getForeignIDFromNode(node);

	  bool insert = 
	    this->writable_node_map.insert(std::map<unsigned int, Node*>::value_type(node_id, node)).second;
	  assert (insert == true);
	}
    }
}


void GeometricModel::writeFESystemInputFile(std::ostream& output)
{
  // 

  // write the nodes

  // for elements, only those elements will be written for which a tag has been
  // specified in the input file

  // next, for loads, check the tags for which the loads need to be written
  // and write them one by one

}



void GeometricModel::writeNastranInputFile(std::ostream& output)
{

}



void GeometricModel::readModelConfigurationInputFile(std::istream& input)
{
  // the information needed from this file includes:
  // -- name of the input file which includes the model data
  // -- name of the file which includes the mesh data
  // -- sets to be included
  // -- sets and their corresponding materials
  // -- sets on which boundary conditions will be applied, and the kind of BC
  // -- sets on which loads will be included, and the kind of loads ????

  std::string tag; 
  unsigned int num = 0;

  // read the gmsh geometry data file name first
  tag.clear();
  input >> tag;
  assert (tag == "GEOMETRY_INPUT_FILE");

  input >> this->geometry_file;



  // next read the mesh data file name
  tag.clear();
  input >> tag;
  assert (tag == "MESH_INPUT_FILE");
  
  input >> this->mesh_file;
  
  // next read the sets to be included in the geometry
  tag.clear();
  input >> tag;
  assert (tag == "BEGIN_PHYSICAL_GEOMETRY_SETS");

  input >> num;
  for (unsigned int i=0; i < num; i++)
    {
      tag.clear();
      input >> tag;
      this->writable_set_vector.push_back(tag);	
    }


  tag.clear();
  input >> tag;
  assert (tag == "END_PHYSICAL_GEOMETRY_SETS");


  // next read the sets and the material properties
  tag.clear();
  input >> tag;
  assert (tag == "BEGIN_MATERIAL_PROPERTY_SETS");

  input >> num;
  unsigned int mat_id = 0;
  for (unsigned int i=0; i < num; i++)
    {
      tag.clear();
      input >> tag;
      input >> mat_id;
      bool insert = 
	this->set_material_ID_map.insert(std::map<std::string, unsigned int>::value_type(tag, mat_id)).second;
      assert (insert == true);
    }

  tag.clear();
  input >> tag;
  assert (tag == "END_MATERIAL_PROPERTY_SETS");

  // next read the sets and the boundary conditions on them
  tag.clear();
  input >> tag;
  assert (tag == "BEGIN_BOUNDARY_CONDITION_SETS");
  
  unsigned int n_dofs;

  input >> num;
  for (unsigned int i=0; i< num; i++)
    {
      tag.clear();
      input >> tag;
      input >> n_dofs;

      Load load(Load::BOUNDARY_CONDITION, n_dofs);
      input >> load;

      bool insert = 
	this->set_boundary_condition_map.insert(std::map<std::string, Load>::value_type(tag, load)).second;
      assert (insert == true);
    }


  tag.clear();
  input >> tag;
  assert (tag == "END_BOUNDARY_CONDITION_SETS");
    

  // next read the sets and the kind of loads on them ???? 
}
