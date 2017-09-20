// $Id: RadiationCavity.h,v 1.10.6.1 2007-06-13 14:59:19 manav Exp $

#ifndef __fesystem_radiation_cavity_h__
#define __fesystem_radiation_cavity_h__

// C++ includes
#include <iostream>
#include <map>
#include <vector>
#include <memory>

// FESystem includes
#include "Radiation/RadiationElement.h"
#include "Mesh/FEMeshData.h"
#include "Radiation/RadiationCavityFiniteElement.h"

// Forward declearations
class Node;
class Elem;


typedef std::vector<RadiationElement*>::const_iterator 
RadiationElemConstIterator;


class RadiationCavity
{
public:
  /// constructor
  RadiationCavity();

  /// destructor
  ~RadiationCavity();


  /// method to read from input stream
  std::istream& readFromInputStream(std::istream& );
  
  /// returns number of finite elems in the cavity
  inline unsigned int nFiniteElems() const;
  
  /// returns the cavity ID
  inline unsigned int getCavityID() const;
  

  /// sets the finite element mesh ID for the specific DV. If no
  /// DV ID is specified, then the mesh ID is set for the base 
  /// design point
  inline void setFiniteElementMeshID(unsigned int mesh_ID,
				     unsigned int DV_ID = 0);


  /// returns the finite element mesh ID for the specified DV. If no
  /// DV ID is specified, then the mesh ID of the base design point is 
  /// returned
  inline unsigned int getFiniteElementMeshID(unsigned int DV_ID=0) const;
  
  /// returns the number of radiation elems in the cavity
  inline unsigned int nRadiationElems() const;
  
  /// returns the nodes for finite element for the specified design variable.
  /// If no design variable is specified, then the nodes for the base mesh 
  /// is returned
  inline const std::vector<Node*>& getFENodes(unsigned int DV_ID =0) const;
  
  /// returns the finite element from the elem ID for the specified 
  /// design variable. If no DV is specified, then it returns an element
  /// for the base cavity.
  /// @param elem_ID ID of the elem to be retuned
  /// @param DV_ID  (optional) ID of the DV for which this quantity is needed
  inline const RadiationElement& 
    getElemFromID(const unsigned int elem_ID,
		  const unsigned int DV_ID=0) const;
  
  /// returns the radiation cavity finite element for the given ID
  inline RadiationCavityFiniteElem& 
    getRadiationCavityFiniteElement(unsigned int elem_ID);
  
  /// returns the geometric elem given by the FE ID, and DV ID.
  /// if no DV is specified, then the elem for base mesh is returned
  /// @param elem_ID ID of the FE that is being requested
  /// @param DV_ID design variable ID
  inline const Elem* getGeometricFE(unsigned int elem_ID,
				    unsigned int DV_ID=0);

  /// returns the geometric elem given by the radiation elem ID for the
  /// DV_ID. If no design variable is specified, then the element
  /// for the base mesh is returned
  /// @param elem_ID of the radiation element that is being requested
  /// @param DV_ID ID of the design variable for which the elem is to be 
  /// returned
  inline Elem* getGeometricRadiationElem(unsigned int elem_ID,
					 unsigned int DV_ID=0);

  /// returns the beginning iterator to the radiation elems of the mesh 
  /// corresponding to the DV. If no DV is specified, then the one for 
  /// base mesh is returned
  inline RadiationElemConstIterator getElemBeginIterator(unsigned int DV_ID =0) const;

  /// returns the end iterator to the radiation elems of the mesh
  /// corresponding to the DV. If no DV is specified, then the one for 
  /// base mesh is returned
  inline RadiationElemConstIterator getElemEndIterator(unsigned int DV_ID =0) const;
	
  /// returns the vector of internal node IDs for the given 
  /// internal rad elem ID
  /// @param rad_ID
  /// @param vector the node IDs will be returned in this 
  void getInternalNodeIDForRadElem(unsigned int rad_ID,
				   std::vector<unsigned int>& vector);

  /// sets the finite element mesh data for the cavity
  /// @param mesh_data this is the mesh data corresponding to 
  /// the mesh that is participating in the radiation. The mesh is 
  /// specified for a DV, and if no DV is passed in the arguement, then
  /// it is assumed for the base mesh
  inline void setMeshData(MeshDS::FEMeshData* mesh_data,
			  unsigned int DV_ID =0);
	

  /// returns a pointer the mesh data object for this cavity, qualified
  /// by the design variable.
  inline MeshDS::FEMeshData* getMeshData(unsigned int DV_ID=0) const;

  
  /// boolean about whether to approximate or recalculate ABInv factor
  bool ifApproximateABInvFactor() const;
  
  /// @returns the number of iterations before the ABInv factor should be recalculated
  unsigned int getNItertionsToRecalculateABInv() const;
  
  /// @returns the number of Taylor series terms to be used for approximation of the
  /// ABInv factor
  unsigned int getNTaylorTermsForABInvApproximation() const;

 protected:

  /// method creates a vector of nodes for finite elements for the 
  /// specified DV.
  void createFENodeVectorAndRadiationElems(unsigned int DV_ID);
  
  /// ID of this radiation cavity
  unsigned int radiation_cavity_ID;

  /// if ABinv factor should be approximated? instead of recalculated
  bool if_approximate_AB_inv;
  
  /// number of terms to use in the approximation
  unsigned int n_AB_inv_taylor_approx_terms;
  
  /// number of iterations before fully recalculating
  unsigned int n_iters_before_AB_inv_recalculation;
  
  /// map of finite element mesh ID in which radiation is occurring
  /// the key used for this map is the design variable, where, the 
  /// mesh ID for each DV perturbation is stored against the DV ID
  std::map<unsigned int, unsigned int>  fe_mesh_ID_map;
  
  /// checks if the cavity data has been initialized or not
  bool initialized;
  
  /// map of mesh data mesh data stored against the shape DV ID
  std::map<unsigned int ,MeshDS::FEMeshData*> mesh_data_map;

  
  /// map of the finite element ID to elem. This is the same for 
  /// each DV
  std::map<unsigned int, RadiationCavityFiniteElem*> ID_to_elem_map;
  
  //  /// vector of radiation elems in this cavity
  //  std::map<unsigned int, RadiationElement> radiation_element;
  
  /// map of foreign node ID to local node ID
  std::map<unsigned int, unsigned int > foreign_to_local_node_ID;
  
  /// vector of local to foreign node ID. Local node ID is the location of 
  /// the ID in this vector
  std::vector<unsigned int> local_to_foreign_node_ID;
  
  /// vector of IDs. The internal ID is the same as the location of the node
  /// in this vector. One vector is stored for each design variable
  std::map<unsigned int, std::vector<Node*> > local_ID_to_node_map;


  /// vector of radiation elements the local position of the elem in this 
  /// vector is also its internal ID. This is stored as a map for each DV
  std::map<unsigned int, std::vector<RadiationElement*> > rad_elem_vector_map;
  
};




inline
unsigned int RadiationCavity::nFiniteElems() const
{
  return this->ID_to_elem_map.size();	
}




inline
unsigned int RadiationCavity::getCavityID() const
{
  return this->radiation_cavity_ID;
}



inline 
void RadiationCavity::setFiniteElementMeshID(unsigned int mesh_ID,
					     unsigned int DV_ID) 
{
  std::map<unsigned int, unsigned int>::const_iterator it = 
    this->fe_mesh_ID_map.find(DV_ID);
  assert (it == this->fe_mesh_ID_map.end());

  bool insert_bool = 
    this->fe_mesh_ID_map.insert(std::map<unsigned int, unsigned int>::value_type
				(DV_ID, mesh_ID)).second;

  assert ( insert_bool == true);
}


inline
unsigned int RadiationCavity::getFiniteElementMeshID(unsigned int DV_ID) const
{
  std::map<unsigned int, unsigned int>::const_iterator it =
    this->fe_mesh_ID_map.find(DV_ID);
  assert (it != this->fe_mesh_ID_map.end());
  
  return it->second;
}





inline const RadiationElement& 
RadiationCavity::getElemFromID(const unsigned int elem_ID,
			       unsigned int DV_ID) const
{
  // make sure that this ID exists in the map
  std::map<unsigned int, std::vector<RadiationElement*> >::const_iterator 
    it = this->rad_elem_vector_map.find(DV_ID);
  assert (it != this->rad_elem_vector_map.end());

  assert (elem_ID <= it->second.size());

  return *(it->second[elem_ID]);
}




inline
RadiationElemConstIterator 
RadiationCavity::getElemBeginIterator(unsigned int DV_ID) const
{
  std::map<unsigned int, std::vector<RadiationElement*> >::const_iterator it =
    this->rad_elem_vector_map.find(DV_ID);
  
  assert ( it != this->rad_elem_vector_map.end());
  
  return it->second.begin();
}





inline
RadiationElemConstIterator 
RadiationCavity::getElemEndIterator(unsigned int DV_ID) const
{
  std::map<unsigned int, std::vector<RadiationElement*> >::const_iterator it =
    this->rad_elem_vector_map.find(DV_ID);
  
  assert ( it != this->rad_elem_vector_map.end());
  
  return  it->second.end();
}



inline 
void RadiationCavity::setMeshData(MeshDS::FEMeshData* meshdata,
				  unsigned int DV_ID)
{
  // the mesh data can be set only once
  assert (meshdata != NULL);
  
  std::map<unsigned int, MeshDS::FEMeshData*>::const_iterator it = 
    this->mesh_data_map.find(DV_ID);
  assert (it == this->mesh_data_map.end());

  bool insert = 
    this->mesh_data_map.insert(std::map<unsigned int, MeshDS::FEMeshData*>::value_type
			       (DV_ID, meshdata)).second;

  assert ( insert == true);
  
  this->createFENodeVectorAndRadiationElems(DV_ID);
}



inline 
MeshDS::FEMeshData* RadiationCavity::getMeshData(unsigned int DV_ID) const
{
  // make sure that the mesh data is still not null
  std::map<unsigned int, MeshDS::FEMeshData*>::const_iterator it = 
    this->mesh_data_map.find(DV_ID);
  assert (it != this->mesh_data_map.end());

  assert (it->second != NULL);

  return it->second;
}




inline const Elem* 
RadiationCavity::getGeometricFE(unsigned int elem_ID,
				unsigned int DV_ID)
{
  MeshDS::FEMeshData* mesh_data = this->getMeshData(DV_ID);
  
  return mesh_data->getElemFromForeignID(elem_ID);
}


inline RadiationCavityFiniteElem& 
RadiationCavity::getRadiationCavityFiniteElement(unsigned int elem_ID)
{
  std::map<unsigned int, RadiationCavityFiniteElem*>::const_iterator it;
  it = this->ID_to_elem_map.find(elem_ID);
  assert (it != this->ID_to_elem_map.end());

  return *(it->second);
}



inline Elem* 
RadiationCavity::getGeometricRadiationElem(unsigned int elem_ID,
					   unsigned int DV_ID)
{
  assert (this->initialized);
  
  std::map<unsigned int, std::vector<RadiationElement*> >::const_iterator 
    it = this->rad_elem_vector_map.find(DV_ID);

  assert ( it != this->rad_elem_vector_map.end());

  return (it->second[elem_ID])->getRadiationGeometricElem();
}


inline 
unsigned int RadiationCavity::nRadiationElems() const
{
  // return the size of the base element vector
  assert (this->initialized);
  std::map<unsigned int, std::vector<RadiationElement*> >::const_iterator 
    it = this->rad_elem_vector_map.find(0);
  return it->second.size();
}


inline
const std::vector<Node*>& RadiationCavity::getFENodes(unsigned int DV_ID) const
{
  assert (this->initialized);
  std::map<unsigned int, std::vector<Node*> >::const_iterator it = 
    this->local_ID_to_node_map.find(DV_ID);
  assert (it != this->local_ID_to_node_map.end());
  
  return it->second;
}


inline
bool RadiationCavity::ifApproximateABInvFactor() const
{
  assert (this->initialized);
  return this->if_approximate_AB_inv;
}



inline
unsigned int RadiationCavity::getNItertionsToRecalculateABInv() const
{
  assert (this->initialized);
  return this->n_iters_before_AB_inv_recalculation;
}


inline
unsigned int RadiationCavity::getNTaylorTermsForABInvApproximation() const
{
  assert (this->initialized);
  return this->n_AB_inv_taylor_approx_terms;
}



#endif // __radiation_cavity_h__

