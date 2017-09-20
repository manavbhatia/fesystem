// $Id: RadiationCavity.C,v 1.10.6.2 2008-08-21 00:56:19 manav Exp $


// C++ includes
#include <iostream>
#include <string>



// FESystem includes
#include "Radiation/RadiationCavity.h"
#include "Utilities/InputOutputUtility.h"


// libMesh include
#include "geom/node.h"
#include "geom/elem.h"

RadiationCavity::RadiationCavity():
  radiation_cavity_ID(0),
if_approximate_AB_inv(false),
n_AB_inv_taylor_approx_terms(FESystemNumbers::InvalidID),
n_iters_before_AB_inv_recalculation(FESystemNumbers::InvalidID),
initialized(false)
{
  
}




RadiationCavity::~RadiationCavity()
{
  // iterate over all the elements, and delete all of them
  std::map<unsigned int, RadiationCavityFiniteElem*>::iterator 
    elem_it , elem_end;
  elem_it = this->ID_to_elem_map.begin();
  elem_end = this->ID_to_elem_map.end();

	
  for (; elem_it != elem_end; elem_it++)
    {
      delete elem_it->second;
      elem_it->second = NULL;
    }
  
  this->ID_to_elem_map.clear();
  
  // node delete the rad elems
  std::map<unsigned int, std::vector<RadiationElement*> >::iterator 
    rad_map_it, rad_map_end;

  rad_map_it = this->rad_elem_vector_map.begin();
  rad_map_end = this->rad_elem_vector_map.end();
  
  for ( ; rad_map_it != rad_map_end; rad_map_it++)
    {
      std::vector<RadiationElement*>::iterator rad_elem_it, rad_elem_end;
      rad_elem_it = rad_map_it->second.begin();
      rad_elem_end = rad_map_it->second.end();
      
      for (; rad_elem_it != rad_elem_end; rad_elem_it++)
	{
	  delete (*rad_elem_it);
	}
    
      rad_map_it->second.clear();
    }
  this->rad_elem_vector_map.clear();
}





std::istream& RadiationCavity::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_cards=0;
	
  FESystemIOUtility::readFromInput(input, "RADIATION_CAVITY");
  FESystemIOUtility::readFromInput(input, "BEGIN");

  FESystemIOUtility::readFromInput(input, "RADIATION_CAVITY_ID", this->radiation_cavity_ID);
	
  unsigned int id = 0;
  FESystemIOUtility::readFromInput(input, "FINITE_ELEMENT_MESH_ID", id);
  this->setFiniteElementMeshID(id);

  FESystemIOUtility::readFromInput(input, "N_ELEMENTS", n_cards);

  unsigned int internal_ID = 0;
	
  for (unsigned int i=0; i<n_cards; i++)
    {
      // read in all the element
      RadiationCavityFiniteElem* elem = new RadiationCavityFiniteElem();
      assert (elem != NULL);
		
      input >> elem->elem_ID;
      input >> elem->side_num;
      input >> elem->n_x_divs;
      input >> elem->n_y_divs;
      elem->internal_ID = internal_ID;
      internal_ID ++;
      
      
      std::pair<std::map<unsigned int, RadiationCavityFiniteElem*>::const_iterator, bool> 
	insert_return_pair =
	this->ID_to_elem_map.insert(std::map<unsigned int, RadiationCavityFiniteElem*>::value_type
				    (elem->elem_ID, elem));
      
      assert (insert_return_pair.second == true);
      
    }
  
  FESystemIOUtility::readFromInput(input, "IF_AB_INVERSE_APPROXIMATE",
                                   this->if_approximate_AB_inv);
  
  FESystemIOUtility::readFromInput(input, "N_AB_INV_TAYLOR_TERMS",
                                   this->n_AB_inv_taylor_approx_terms);

  FESystemIOUtility::readFromInput(input, "IF_AB_INV_APPROX_BEFORE_RECALCULATE",
                                   this->n_iters_before_AB_inv_recalculation);

  FESystemIOUtility::readFromInput(input, "RADIATION_CAVITY");
  FESystemIOUtility::readFromInput(input, "END");
  
  
  return input;
}	





void RadiationCavity::createFENodeVectorAndRadiationElems(unsigned int DV_ID)
{
  // iterate over all finite elems in the cavity
  // get the nodes for the elem. Insert them in the node set
  // once the set has been prepared, put it in the node vector
  // and record its location with each
 

  std::vector<Node*>* local_ID_to_node_vec = NULL;
  std::vector<RadiationElement*>* rad_elem_vec = NULL;
  MeshDS::FEMeshData* mesh_data = this->getMeshData(DV_ID);
  {
    std::map<unsigned int, std::vector<Node*> >::const_iterator node_map_it =
      this->local_ID_to_node_map.find(DV_ID);
    assert (node_map_it == this->local_ID_to_node_map.end());

    assert(this->local_ID_to_node_map.insert
	   (std::map<unsigned int, std::vector<Node*> >::value_type
	    (DV_ID, std::vector<Node*>())).second == true);

    local_ID_to_node_vec = &(this->local_ID_to_node_map[DV_ID]);

    std::map<unsigned int, std::vector<RadiationElement*> >::const_iterator rad_elem_map_it =
      this->rad_elem_vector_map.find(DV_ID);
    assert (rad_elem_map_it == this->rad_elem_vector_map.end());

    assert(this->rad_elem_vector_map.insert
	   (std::map<unsigned int, std::vector<RadiationElement*> >::value_type
	    (DV_ID, std::vector<RadiationElement*>())).second == true);

    rad_elem_vec = &(this->rad_elem_vector_map[DV_ID]);
  }

  Node* node = NULL;
  unsigned int foreign_node_ID = 0, internal_node_ID = 0;
  
  // only for the base design point, create the vector of local to global and
  // global to local ID maps. Once this is done, the map can be used to 
  // create the vector of nodes
  if (DV_ID == 0)
    {
      std::map<unsigned int, RadiationCavityFiniteElem*>::const_iterator it, end;
      it = this->ID_to_elem_map.begin();
      end = this->ID_to_elem_map.end();

      std::map<unsigned int, unsigned int>::iterator node_it;

      for (; it != end; it++)
	{
	  const Elem* elem = 
	    mesh_data->getElemFromForeignID(it->second->elem_ID);
      
	  // the elem can only be a 2-D planar elem
	  assert (elem->dim() < 3);
      
	  // iterate over the nodes for this elem
	  for (unsigned int i=0; i < elem->n_nodes(); i++)
	    {
	      node = elem->get_node(i);
	      foreign_node_ID = mesh_data->getForeignIDFromNode(node);
	  
	      // if this node does not already exist in the map, insert it
	      node_it = this->foreign_to_local_node_ID.find(foreign_node_ID);
	      if (node_it == this->foreign_to_local_node_ID.end())
		{
		  // push the node in the vectors
		  this->local_to_foreign_node_ID.push_back(foreign_node_ID);

		  bool insert = 
		    this->foreign_to_local_node_ID.insert
		    (std::map<unsigned int, unsigned int>::value_type
		     (foreign_node_ID, internal_node_ID)).second;
		  assert (insert == true);
	      
		  internal_node_ID++;
		}  
	    }
	}
    }

  // now iterate over all the nodes and create the node vector
  {
    std::vector<unsigned int>::const_iterator node_ID_it, node_ID_end;
    node_ID_it = this->local_to_foreign_node_ID.begin();
    node_ID_end = this->local_to_foreign_node_ID.end();
    
    for ( ; node_ID_it != node_ID_end; node_ID_it++)
      {
	node =  const_cast<Node*>(mesh_data->getNodeFromForeignID(*node_ID_it));
	local_ID_to_node_vec->push_back(node);
      }
  }

  // create the radiation elements
  {
    std::map<unsigned int, RadiationCavityFiniteElem*>::const_iterator it, end;
    it = this->ID_to_elem_map.begin();
    end = this->ID_to_elem_map.end();

    RadiationElement* rad_elem = NULL;
    unsigned int local_rad_elem_ID = 0;

    for (; it != end; it++)
      {
	// for each elem, get the number of divisions in x and y direction
	// iterate over each row of divisions and create the radiation elems
	// calculate the local coordinates of the elem in the FE, and 
	// also the global coordinates of the elem vertices

	const Elem* elem = 
	  mesh_data->getElemFromForeignID(it->second->elem_ID);

	for (unsigned int i=0; i < it->second->n_x_divs; i++) // increment along x
	  for (unsigned int j = 0; j < it->second->n_y_divs; j++) // increment along y
	    {
	      // create the local radiation elem
	      rad_elem = new RadiationElement(*(it->second), elem);
	      rad_elem->setData(local_rad_elem_ID,
				i, j);
	    
	    
	      // insert the radiation elem in the elem vector
	      rad_elem_vec->push_back(rad_elem);

	      local_rad_elem_ID++;
	    }
      }
  }  
  this->initialized = true;  
}



void RadiationCavity::getInternalNodeIDForRadElem(unsigned int rad_ID,
						  std::vector<unsigned int>& node_vec)
{
  // get the FE elem ID for this rad elem
  unsigned int elem_ID = 
    (this->rad_elem_vector_map[0])[rad_ID]->FeElementID(),
    node_external_ID = 0;

  const Elem* elem = 
    (this->mesh_data_map[0])->getElemFromForeignID(elem_ID);

  std::map<unsigned int, unsigned int>::const_iterator it, end;
  end = this->foreign_to_local_node_ID.end();
  
  // iterate over the nodes in this elem, and get their node IDs
  // from the ID map
  node_vec.clear();
  for (unsigned int i=0; i<elem->n_nodes(); i++)
    {
      Node* node = elem->get_node(i);
      node_external_ID = 
	(this->mesh_data_map[0])->getForeignIDFromNode(node);

      it = this->foreign_to_local_node_ID.find(node_external_ID);
      assert (it != end);

      node_vec.push_back(it->second);
    }
}
