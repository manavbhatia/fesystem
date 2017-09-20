// $Id: MeshList.C,v 1.17.6.1 2008-02-25 04:27:01 manav Exp $

// C++ includes


// FESystem inlcudes
#include "Utilities/InputOutputUtility.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Radiation/RadiationCavity.h"


//libMesh includes
#include "mesh.h"
#include "dof_map.h"


//different element kinds
#include "edge_edge2.h"
#include "face_quad4.h"
#include "face_tri3.h"
#include "cell_tet4.h"
#include "cell_hex8.h"
#include "cell_prism6.h"



MeshDS::MeshList::MeshList()
{
  
}


MeshDS::MeshList::~MeshList()
{
  this->clear();
}


void 
MeshDS::MeshList::clear()
{
  // get the iterator for the meshID map
  MeshDS::MeshList::IDToMeshPairMap::iterator map_it = this->meshpair_ID_map.begin();
  MeshDS::MeshList::IDToMeshPairMap::iterator map_end = this->meshpair_ID_map.end();
	
  // iterate over all the mesh and meshdata pair
  for (; map_it != map_end; map_it++)
    {
    // delete the pointers
    delete map_it->second.first; map_it->second.first = NULL;
    delete map_it->second.second; map_it->second.second = NULL;
    }
	
  // finally, clear the map
  this->meshpair_ID_map.clear();
  
  // clear the DofMap map
  MeshDS::MeshList::IDToDofMapMap::iterator dof_map_it, dof_map_end;
  dof_map_it = this->mesh_ID_to_dof_map.begin();
  dof_map_end = this->mesh_ID_to_dof_map.end();
  
  
  for (; dof_map_it != dof_map_end; dof_map_it++)
    {
    delete dof_map_it->second; dof_map_it->second = NULL;
    }
  this->mesh_ID_to_dof_map.clear();
	
  // also clear the radiation cavities
	
  MeshDS::MeshList::IDToRadiationCavityMap::iterator rad_map_it, rad_map_end;
  rad_map_it = this->radiation_cavity_ID_map.begin();
  rad_map_end = this->radiation_cavity_ID_map.end();
	
  for (; rad_map_it != rad_map_end; rad_map_it++)
    {
    delete rad_map_it->second; rad_map_it->second = NULL;
    }
	
  this->radiation_cavity_ID_map.clear();
	
  // iterate over all the dummy elem sets and delete them
  MeshDS::MeshList::IDToElemSetMap::iterator elem_set_it, elem_set_end;
  elem_set_it = this->dummy_elem_set_map.begin();
  elem_set_end = this->dummy_elem_set_map.end();
  
  for (; elem_set_it != elem_set_end; elem_set_it++)
    {
    delete elem_set_it->second;
    elem_set_it->second = NULL;
    }
  
  this->dummy_elem_set_map.clear();
}



std::istream& 
MeshDS::MeshList::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num_cards = 0;
  unsigned int num_mesh = 0; // number of mesh in the file
  unsigned int mesh_ID = 0, dimension = 0;
  
	
  // the first tag should be BEGIN_MESH, followed by the number of mesh
  // in the file
  FESystemIOUtility::readFromInput(input, "MESH_LIST");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "N_MESH", num_mesh);
  
	
  // iterate over all the mesh and read in the data for each
  for (unsigned int mesh_it=0; mesh_it < num_mesh; mesh_it++)
    {
		
    // next, start reading the data for each mesh
    // the first tag should be BEGIN_MESH, next, read the ID for this mesh
    FESystemIOUtility::readFromInput(input, "MESH");
    FESystemIOUtility::readFromInput(input, "BEGIN");
    FESystemIOUtility::readFromInput(input, "MESH_ID", mesh_ID);
    FESystemIOUtility::readFromInput(input, "DIMENSION", dimension);
    AssertThrow(dimension <= 3, ExcInternalError());
		
    // now add a mesh and meshdata pair to this mesh list
    this->addNewMeshDataStructures(mesh_ID, dimension);
		
    MeshDS::FEMesh* mesh = this->getMeshFromID(mesh_ID);
    MeshDS::FEMeshData* mesh_data = this->getMeshDataFromID(mesh_ID);
				
    // next, start reading the nodes
    // it should start with the BEGIN_NODES tag	
    FESystemIOUtility::readFromInput(input, "NODES");
    FESystemIOUtility::readFromInput(input, "BEGIN");
    FESystemIOUtility::readFromInput(input, "N_NODES", num_cards);
			
		
    //loop over all the nodes and read the data
    {
      unsigned int node_id = 0; double x_loc = 0.0, y_loc = 0.0, z_loc = 0.0;
      for (unsigned int i = 0; i < num_cards; i++)
        {
        //read the nodal data
        input>> node_id;	// node ID
        input>> x_loc;	// x location of the node
        input>> y_loc;	// y location
        input>> z_loc;	// z location
				
        //create the point and add it to the mesh for this FEsystem
        Point this_point(x_loc, y_loc, z_loc);
        Node* this_node = mesh->add_point(this_point);
				
        //use the ID of the node and create an entry in mesh_data
        mesh_data->addNodeAndIDData(this_node, node_id);
        }
    }
		
		
    //make sure that the nodal data ends with the END_NODES tag
    FESystemIOUtility::readFromInput(input, "NODES");
    FESystemIOUtility::readFromInput(input, "END");
		
		
    //******************	read the element data******************
    //
    //it should start with the BEGIN_ELEMENTS tag 
    FESystemIOUtility::readFromInput(input, "ELEMENTS");
    FESystemIOUtility::readFromInput(input, "BEGIN");
    FESystemIOUtility::readFromInput(input, "N_ELEMENTS", num_cards);
		
    
    //loop over all the elements and read the data
    {
      unsigned int elem_id, prop_id, node_id, elem_type_enum_ID;
      Elem* this_elem;
      ElemType elem_type;
      for (unsigned int i = 0; i < num_cards; i++)
        {
        //read the element data
        tag.clear();
        input >> tag;	// element tag
        elem_type_enum_ID = FESystemElem::FESystemElemTypeEnum::enumID(tag);
        elem_type = FESystemElem::FESystemElemTypeEnum::elemType(tag);
        
        
        switch (elem_type)
          {
          case EDGE2:
            this_elem = mesh->add_elem(new Edge2);
            break;
            
          case QUAD4:
            this_elem = mesh->add_elem(new Quad4);
            break;
            
          case TRI3:
            this_elem = mesh->add_elem(new Tri3);
            break;
            
          case PRISM6:
            this_elem = mesh->add_elem(new Prism6);
            break;
            
          case HEX8:
            this_elem = mesh->add_elem(new Hex8);
            break;
            
          default:
            abort();
          }
        
        				
        //now, repeat and read the ID, property ID and connectivity
        input >> elem_id;	// element ID
        input >> prop_id;	// property ID
				
        //read the node IDs and set the element connectivity
        for (unsigned int j=0; j < this_elem->n_nodes(); j++)
          {
          input >> node_id;
          this_elem->set_node(j) = 
            const_cast<Node*>(mesh_data->getNodeFromForeignID(node_id));
          }
				
        //next, insert the element and its ID into mesh data
        mesh_data->addElemAndIDData(this_elem, elem_id, prop_id, elem_type_enum_ID);
        }
    }
		
		
		
    //it should end with the END_ELEMENTS tag
    FESystemIOUtility::readFromInput(input, "ELEMENTS");
    FESystemIOUtility::readFromInput(input, "END");
		
		
    FESystemIOUtility::readFromInput(input, "MESH");
    FESystemIOUtility::readFromInput(input, "END");
		
    // now that the elements have been read,  close the foreign id maps for elems
    mesh->prepare_for_use();
    }
	
	
  FESystemIOUtility::readFromInput(input, "MESH_LIST");
  FESystemIOUtility::readFromInput(input, "END");
	
	
  // next, read in the radiation cavities
  FESystemIOUtility::readFromInput(input, "RADIATION_CAVITY_LIST");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "N_RADIATION_CAVITIES", num_cards);
	
  for (unsigned int i=0; i < num_cards; i++)
    {
    RadiationCavity *cavity = new RadiationCavity();
    cavity->readFromInputStream(input);
    this->addRadiationCavity(cavity);	
    }
		
  FESystemIOUtility::readFromInput(input, "RADIATION_CAVITY_LIST");
  FESystemIOUtility::readFromInput(input, "END");
	
	
  
  // add the dof map for each mesh
  MeshDS::MeshList::IDToMeshPairMap::iterator mesh_it, mesh_end;
  mesh_it = this->meshpair_ID_map.begin();
  mesh_end = this->meshpair_ID_map.end();
  
  // this is the internal mesh ID, since the IDs to init dofmap need to be 
  // the internal ID.
//  const unsigned int system_ID = 0; // all mesh have the same system ID = 0
  bool insert = false;
  
  for (; mesh_it != mesh_end; mesh_it++)
    {
    mesh_ID = mesh_it->first;
    MeshDS::FEDofMap* dof_map = new MeshDS::FEDofMap();
    insert = 
      this->mesh_ID_to_dof_map.insert(MeshDS::MeshList::IDToDofMapMap::value_type
                                      (mesh_ID, dof_map)).second;
    assert (insert == true);
    }
  
  return input;
}




void 
MeshDS::MeshList::addNewMeshDataStructures(const unsigned int mesh_ID,
                                           const unsigned int dimension)
{
  // first creat a new mesh and meshdata pointer, and then add its pair to 
  // the map
  MeshDS::FEMesh* mesh = new MeshDS::FEMesh(dimension);
  MeshDS::FEMeshData* mesh_data = new MeshDS::FEMeshData(*mesh);
	
  // now add this to the map and make sure that the map was added
  MeshDS::MeshList::MeshAndMeshDataPair pair(mesh,mesh_data);
  bool insert_success = this->meshpair_ID_map.insert
    (MeshDS::MeshList::IDToMeshPairMap::value_type(mesh_ID, pair)).second;
  
  Assert(insert_success,
         FESystemExceptions::ExcDuplicateID("FEMesh", mesh_ID));
  

  
  //  // also add a dummy set for this mesh
  std::auto_ptr<std::set<Elem*> > elem_set(new std::set<Elem*>);
  insert_success = 
    this->dummy_elem_set_map.insert(MeshDS::MeshList::IDToElemSetMap::
                                    value_type(mesh_ID, elem_set.release())).second;

  // this should always be true
  Assert(insert_success,
         FESystemExceptions::ExcDuplicateID("Mesh" ,mesh_ID));
}





void 
MeshDS::MeshList::addRadiationCavity(RadiationCavity* rad_cavity)
{
  const unsigned int ID = rad_cavity->getCavityID();
		
  bool insert_success =
    this->radiation_cavity_ID_map.insert(MeshDS::MeshList::IDToRadiationCavityMap::value_type
                                         (ID, rad_cavity)).second;
	
  Assert(insert_success,
         FESystemExceptions::ExcDuplicateID("Radiation Cavity", ID));
}


