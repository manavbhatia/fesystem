// $Id: FESystemElem.C,v 1.14.4.5 2008-06-03 05:19:41 manav Exp $

// C++ includes

// FESystem includes
#include "FESystem/FESystemElem.h"
#include "FESystem/FESystemController.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/PistonTheory.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Database/ElementDataStorage.h"
#include "Database/GlobalDataStorage.h"
#include "Properties/ElemDataCard.h"
#include "PostProcess/ElemPostProcessQty.h"

#include "StructuralElems/plate_MITC4.h"
#include "StructuralElems/membrane_tri3.h"
#include "StructuralElems/membrane_quad4.h"
#include "StructuralElems/bar.h"
#include "StructuralElems/beam2.h"
#include "StructuralElems/brick_hex8.h"
#include "StructuralElems/linear_spring.h"
#include "StructuralElems/PlateDKT.h"
#include "StructuralElems/Tri3VonKarman.h"

#include "ThermalElems/conduction_1d.h"
#include "ThermalElems/conduction_hex8.h"
#include "ThermalElems/conduction_prism6.h"
#include "ThermalElems/conduction_quad4.h"
#include "ThermalElems/conduction_tri3.h"
#include "ThermalElems/ForcedConvection1D.h"

#include "PanelMethods/PistonTheoryTri3.h"



// libMesh includes
#include "numerics/numeric_vector.h"
#include "geom/edge_edge2.h"
#include "geom/face_quad4.h"
#include "geom/face_tri3.h"
#include "utils/string_to_enum.h"


FESystemElem::FESystemElemBase::FESystemElemBase
(const unsigned int dim,
 const unsigned int elem_enum_ID,
 Discipline::AnalysisDisciplineBase& discipline):
analysis_discipline(discipline),
dimension(dim),
elem_type_enum_ID(elem_enum_ID),
discipline_enumID(discipline.getDisciplineEnumID()),
elem_ID(FESystemNumbers::InvalidID),
sensitivity_parameter_ID(FESystemNumbers::InvalidID),
sensitivity_parameter(FESystemNumbers::InvalidID),
perturbation(0.0),
fe_quadrature_data_initialized(false),
elem_is_initialized_for_property_sensitivity(false),
property_card_initialized(false),
//elem_is_initialized_for_post_processing(false),
elem_property_card(NULL),
elem_data_storage(NULL),
volume_loads(NULL),
surface_loads(NULL),
nodal_loads(NULL),
volume_load_sens(NULL),
surface_load_sens(NULL),
nodal_load_sens(NULL)
{  
  // set the element qty database
  this->elem_data_storage = 
    this->analysis_discipline.getFESystemController().element_data_storage.get();
  
  // now initialize the maps
  // the global axis geometric elements
  this->geometric_elems_for_DV_map[FESystemElem::BASE_ELEM::num()] = 
    NULL;
  this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    NULL;

  // the local axis geometric elements
  this->local_elem_for_DV_map[FESystemElem::BASE_ELEM::num()] = 
    NULL;
  this->local_elem_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    NULL;
  
  /// nodes created for the local element 
  this->local_elem_nodes_for_DV_map[FESystemElem::BASE_ELEM::num()] = 
    NULL;
  this->local_elem_nodes_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    NULL;
  
  /// boolean for keeping track of whether this element is initialized or not
  this->elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()] = 
    false;
  this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    false;
  
  /// boolean for keeping track of whether this element is initialized or not
  this->local_elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()] = 
    false;
  this->local_elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    false;
  
  /// keeps track of the initialized side that was initialized. This is stored for each design point.
  this->initialized_side_for_DV_map[FESystemElem::BASE_ELEM::num()] = -1;
  this->initialized_side_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = -1;
  
  /// finite element object for the element domain. This is stored for each design point.
  this->fe_base_map_for_DV[FESystemElem::BASE_ELEM::num()] = 
    std::map<FEFamily, FEBase*>();
  this->fe_base_map_for_DV[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    std::map<FEFamily, FEBase*>();
  
  /// finite element object for the side. This is stored for each design point.
  this->fe_side_map_for_DV[FESystemElem::BASE_ELEM::num()] = 
    std::map<FEFamily, FEBase*>();
  this->fe_side_map_for_DV[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    std::map<FEFamily, FEBase*>();

  /// quadrature rule for the element volume integration. This is stored for each design point.
  this->qbase_map_for_DV[FESystemElem::BASE_ELEM::num()] = 
    std::map<FEFamily, QBase*>();
  this->qbase_map_for_DV[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    std::map<FEFamily, QBase*>();
  
  /// quadrature rule for the side integration. This is stored for each design point.
  this->qbase_side_map_for_DV[FESystemElem::BASE_ELEM::num()] = 
    std::map<FEFamily, QBase*>();
  this->qbase_side_map_for_DV[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
    std::map<FEFamily, QBase*>();

//  /// vector of JxW for all quadrature points on the volume of this element. 
//  /// This is stored for each design point.
//  this->JxW_map_for_DV[FESystemElem::BASE_ELEM::num()] = 
//    std::map<FEFamily, QBase*>();
//  this->JxW_map_for_DV[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = 
//    std::map<FEFamily, QBase*>();
//  std::map<unsigned int, std::map<FEFamily, std::vector<Real>*> > JxW_map_for_DV;
//  
//  /// vector of shape functions for this element at all quadrature points.
//  /// This is stored for each design point.
//  std::map<unsigned int, std::map<FEFamily, std::vector<std::vector<Real> >*> > phi_map_for_DV;
//  
//  /// vector of gradient of shape functions at all quadrature points. 
//  /// This is stored for each design point.
//  std::map<unsigned int, std::map<FEFamily, std::vector<std::vector<RealGradient> >* > > dphi_map_for_DV;
//  
//  /// vector of shape functions on the element side. This is stored for each design point.
//  std::map<unsigned int, std::map<FEFamily, std::vector<std::vector<Real> >* > > phi_face_map_for_DV;
//  
//  /// vector of JxW on the element side. This is stored for each design point.
//  std::map<unsigned int, std::map<FEFamily, std::vector<Real>*> > JxW_face_map_for_DV;
}




FESystemElem::FESystemElemBase::~FESystemElemBase()
{
  // delete the elem only if it is not 3-D, since, in that case 
  // the elem and local_elem pointers are same, and deleting the local 
  // elem will delete elem, which will mess with mesh.
  this->clearSensitivityInitialization();
  this->clearElemInitialization();
  
  this->deleteElemLocalData(FESystemElem::BASE_ELEM::num());
  this->deleteElemLocalData(FESystemElem::BASE_PLUS_DELTA_ELEM::num());
  
  // now iterate over all the fe and qbase maps and delete the objects
  {
    std::map<unsigned int, std::map<FEFamily, FEBase*> >::iterator it, end;
    it = this->fe_base_map_for_DV.begin();
    end = this->fe_base_map_for_DV.end(); 
    
    std::map<FEFamily, FEBase*>::iterator fe_it, fe_end;

    for ( ; it != end; it++)
      {
      fe_it = it->second.begin();
      fe_end = it->second.end();
      for (; fe_it != fe_end; fe_it++)
        {if (fe_it->second != NULL) {delete fe_it->second; fe_it->second = NULL;}}
      }
    
    
    // now repeat the same for element sides
    it = this->fe_side_map_for_DV.begin();
    end = this->fe_side_map_for_DV.end(); 

    for ( ; it != end; it++)
      {
      fe_it = it->second.begin();
      fe_end = it->second.end();
      for (; fe_it != fe_end; fe_it++)
        {if (fe_it->second != NULL) {delete fe_it->second; fe_it->second = NULL;}}
      }
  }
  

  {
    std::map<unsigned int, std::map<FEFamily, QBase*> >::iterator it, end;
    it = this->qbase_map_for_DV.begin();
    end = this->qbase_map_for_DV.end(); 
    
    std::map<FEFamily, QBase*>::iterator qit, qend;
    
    for ( ; it != end; it++)
      {
      qit = it->second.begin();
      qend = it->second.end();
      for (; qit != qend; qit++)
        {if (qit->second != NULL) {delete qit->second; qit->second = NULL;}}
      }
    
    
    // now repeat the same for element sides
    it = this->qbase_side_map_for_DV.begin();
    end = this->qbase_side_map_for_DV.end(); 
    
    for ( ; it != end; it++)
      {
      qit = it->second.begin();
      qend = it->second.end();
      for (; qit != qend; qit++)
        {if (qit->second != NULL) {delete qit->second; qit->second = NULL;}}
      }
  }
}




void 
FESystemElem::FESystemElemBase::initialize_element(const unsigned int design_point_enum_ID,
                                                   const std::vector<Point>* points)
{
  // create local element
  this->create_local_elem(design_point_enum_ID);
	
  // reinit fe with this local element
  std::map<FEFamily, FEBase*>::iterator it, end; 
  it = this->fe_base_map_for_DV[design_point_enum_ID].begin();
  end = this->fe_base_map_for_DV[design_point_enum_ID].end();
  for (; it != end; it++)
    it->second->reinit(this->local_elem_for_DV_map[design_point_enum_ID], points);
	
  this->local_elem_is_initialized_for_DV_map[design_point_enum_ID] = true;
}





void 
FESystemElem::FESystemElemBase::initialize_element_side(unsigned int side_num,
                                                        const unsigned int design_point_enum_ID)
{
  // make sure that the element is initialized, otherwise initialize it
  if(!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID]) 
      this->initialize_element(design_point_enum_ID);
  
  static Elem* local_elem;
  local_elem = this->local_elem_for_DV_map[design_point_enum_ID];
  
  switch (local_elem->dim())
    {
    case 1:
      {
        abort();
      }
      break;
            
    case 2:
    case 3:
      {
        std::map<FEFamily, FEBase*>::iterator it, end; 
        it = this->fe_side_map_for_DV[design_point_enum_ID].begin();
        end = this->fe_side_map_for_DV[design_point_enum_ID].end();
        for (; it != end; it++)
          it->second->reinit(this->local_elem_for_DV_map[design_point_enum_ID], side_num);
        
        this->initialized_side_for_DV_map[design_point_enum_ID] = side_num;
      }
      break;
    }
}





void 
FESystemElem::FESystemElemBase::create_local_elem(const unsigned int design_point_enum_ID)
{
  // check the dimension of the element
  // and call the appropriate function
    
  switch (this->geometric_elems_for_DV_map[design_point_enum_ID]->dim())
    {
    case 1:
      this->create_local_1D_elem(design_point_enum_ID);
      break;
			
    case 2:
      this->create_local_2D_elem(design_point_enum_ID);
      break;
			
    case 3:
      this->local_elem_for_DV_map[design_point_enum_ID] = 
      this->geometric_elems_for_DV_map[design_point_enum_ID];
      break;
			
    default:
      abort();
      break;
    }
}




void 
FESystemElem::FESystemElemBase::create_local_1D_elem(const unsigned int design_point_enum_ID)
{
  Elem* elem = this->geometric_elems_for_DV_map[design_point_enum_ID];
  
  // check if the element was previously created or now. If not, create the element, 
  // otherwise, simply change the nodal coordinates of the previous element
  Node **local_nodes = NULL;
  Elem *local_elem = NULL;
  
  switch (this->local_elem_for_DV_map[design_point_enum_ID] == NULL)
    {
    case true:
      {
        static const Point zero_point(0. , 0. , 0.);
        // create and set the nodes
        local_nodes = new Node*[elem->n_nodes()];
        this->local_elem_nodes_for_DV_map[design_point_enum_ID] = local_nodes;
        
        local_elem = Elem::build(elem->type()).release();
        this->local_elem_for_DV_map[design_point_enum_ID] = local_elem;
        
        for (unsigned int i=0; i < elem->n_nodes(); i++) 
          {
          local_nodes[i] = new Node(zero_point, i);
          local_elem->set_node(i) = local_nodes[i];
          }
        
      }
      break;
      
    case false:
    default:
      {
        // use the same nodes and element
        local_nodes = this->local_elem_nodes_for_DV_map[design_point_enum_ID];
        local_elem = this->local_elem_for_DV_map[design_point_enum_ID];
      }
      break;
    }
  
  switch (elem->type())
    {
    case EDGE2:
      {
        //get the two nodes of this elem
        const Point& point0 = elem->point(0);
        const Point& point1 = elem->point(1);
        
        //calculate the length
        static Point vector01;
        vector01.assign(point1 - point0);
        
        //create a new element and set its two nodes at the
        //origin and the elem length
        static Point new_point0, new_point1;
        new_point0.zero(); new_point1.zero();
        new_point1(0) = vector01.size();
        
        *(local_nodes[0]) = new_point0;
        *(local_nodes[1]) = new_point1;
        
        // now create the local element and set its nodes
        local_elem->set_node(0) = local_nodes[0];
        local_elem->set_node(1) = local_nodes[1];
      }
      break;
      
    default:
      // cases of EDGE3, EDGE4 are not handled here
      abort();
      break;
    }

}





void FESystemElem::FESystemElemBase::create_local_2D_elem(const unsigned int design_point_enum_ID)
{
  Elem* elem = this->geometric_elems_for_DV_map[design_point_enum_ID];

  // initialize the local nodes and elem pointers
  Node **local_nodes = NULL;
  Elem *local_elem = NULL;
  switch (this->local_elem_for_DV_map[design_point_enum_ID] == NULL)
    {
    case true:
      {
        static const Point zero_point(0. , 0. , 0.);
        // create and set the nodes
        local_nodes = new Node*[elem->n_nodes()];
        this->local_elem_nodes_for_DV_map[design_point_enum_ID] = local_nodes;
        
        local_elem = Elem::build(elem->type()).release();
        this->local_elem_for_DV_map[design_point_enum_ID] = local_elem;
        
        for (unsigned int i=0; i < elem->n_nodes(); i++) 
          {
          local_nodes[i] = new Node(zero_point, i);
          local_elem->set_node(i) = local_nodes[i];
          }

      }
      break;
      
    case false:
    default:
      {
        // use the same nodes and element
        local_nodes = this->local_elem_nodes_for_DV_map[design_point_enum_ID];
        local_elem = this->local_elem_for_DV_map[design_point_enum_ID];
      }
      break;
    }
  
  // now initialize the nodes
  switch (elem->type())
    {
    case TRI3:
      {
        //get the three nodes of this elem
        const Point& point0 = elem->point(0);
        const Point& point1 = elem->point(1);
        const Point& point2 = elem->point(2);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, tmp_vector1, tmp_vector2, 
          new_unit_vec_x, new_unit_vec_y;
        
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        
        new_unit_vec_x = vector01.unit();
        
        tmp_vector1.assign(new_unit_vec_x.cross(vector02));
        tmp_vector2 = tmp_vector1.unit();
        
        new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
        
        //create a new element and set its two nodes at the
        //origin and the elem length
        static Point new_point0, new_point1, new_point2;
        new_point0.zero(); new_point1.zero(); new_point2.zero();
        new_point1(0)  = vector01.size();
        new_point2(0) = vector02 * new_unit_vec_x;
        new_point2(1) = vector02 * new_unit_vec_y;
        
        *(local_nodes[0]) = new_point0;
        *(local_nodes[1]) = new_point1;
        *(local_nodes[2]) = new_point2;
      }			
      break;
			
    case QUAD4:
      {
        //get the four nodes of this elem
        const Point& point0 = elem->point(0);
        const Point& point1 = elem->point(1);
        const Point& point2 = elem->point(2);
        const Point& point3 = elem->point(3);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, vector03, tmp_vector1, tmp_vector2,
          new_unit_vec_x, new_unit_vec_y;
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        vector03.assign(point3 - point0);
        
        new_unit_vec_x = vector01.unit();
        
        tmp_vector1.assign(new_unit_vec_x.cross(vector03));
        tmp_vector2 = tmp_vector1.unit();
        
        new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
        
        //create a new element and set its two nodes at the
        //origin and the elem length
        static Point new_point0, new_point1, new_point2, new_point3;
        new_point0.zero(); new_point1.zero(); new_point2.zero(); new_point3.zero();
        
        new_point1(0) = vector01.size();
        
        new_point2(0) = vector02 * new_unit_vec_x; 
        new_point2(1) = vector02 * new_unit_vec_y; 
        
        new_point3(0) = vector03 * new_unit_vec_x; 
        new_point3(1) = vector03 * new_unit_vec_y;
        
        *(local_nodes[0]) = new_point0;
        *(local_nodes[1]) = new_point1;
        *(local_nodes[2]) = new_point2;
        *(local_nodes[3]) = new_point3;
      }
      break;
			
    default:
      // cases for TRI6, QUAD8, QUAD9 are not handled
      abort();
      break;
    }
}





void
FESystemElem::FESystemElemBase::calculate_T_matrix(DenseMatrix<double>* T_mat,
                                                   const unsigned int design_point_enum_ID)
{  
  assert (T_mat != NULL);
  
  // initialize the element if it is not already initialized
  if (!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID])
      this->initialize_element(design_point_enum_ID);
		
  static Elem* local_elem;
  
  local_elem = this->geometric_elems_for_DV_map[design_point_enum_ID];
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  switch (T_mat->m() + T_mat->n() - 2 * (n_nodes * 3))
    {
    case 0:
      {
        // keep going
      }
      break; 
      
    default:
      T_mat->resize(n_nodes * 3, n_nodes * 3);
      break;
    }
  
  T_mat->zero();
  
  switch (local_elem->type())
    {
    case EDGE2:
      {	
        const Point& point0 = local_elem->point(0);
        const Point& point1 = local_elem->point(1);
        
        //calculate the length
        static Point local_x, local_y, local_z, tmp_vector1;
        local_x.assign(point1 - point0); // this is the local x
        local_x = (local_x.unit());
        
        // choose a local y axis
        // if the local_x is parallel to the x-axis, then choose  
        if ( fabs(Point(1.,0.,0.)*local_x - 1.0) <= 1.0e-12 ) 
          {
          // if this is true, then local-x and x are parallel
          // and choose y as the temp vector
          tmp_vector1.assign(Point(0.,1.,0.));
          }
        else
          {
          tmp_vector1.assign(Point(1.,0.,0.));
          }
        
        local_y.assign(local_x.cross(tmp_vector1));
        local_z.assign(local_x.cross(local_y));
        
        // now set the transformation matrix T_mat
        // the dofs are arranged as 
        // (u1,u2,v1,v2,w1,w2,theta_x1,theta_x2,theta_y1,theta_y2,theta_z1,theta_z2)
        
        (*T_mat)(0,0) = local_x(0); (*T_mat)(0,2) = local_x(1); (*T_mat)(0,4) = local_x(2);
        (*T_mat)(1,1) = local_x(0); (*T_mat)(1,3) = local_x(1); (*T_mat)(1,5) = local_x(2);
        
        (*T_mat)(2,0) = local_y(0); (*T_mat)(2,2) = local_y(1); (*T_mat)(2,4) = local_y(2);
        (*T_mat)(3,1) = local_y(0); (*T_mat)(3,3) = local_y(1); (*T_mat)(3,5) = local_y(2);
        
        (*T_mat)(4,0) = local_z(0); (*T_mat)(4,2) = local_z(1); (*T_mat)(4,4) = local_z(2);
        (*T_mat)(5,1) = local_z(0); (*T_mat)(5,3) = local_z(1); (*T_mat)(5,5) = local_z(2);
      }
      break;
			
    case TRI3:
      {
        //get the three nodes of this elem
        const Point& point0 = local_elem->point(0);
        const Point& point1 = local_elem->point(1);
        const Point& point2 = local_elem->point(2);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, tmp_vector1, tmp_vector2,
          new_unit_vec_x, new_unit_vec_y, new_unit_vec_z;
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        
        new_unit_vec_x = vector01.unit();
        
        tmp_vector1.assign(new_unit_vec_x.cross(vector02));
        tmp_vector2 = tmp_vector1.unit();
        
        new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
        
        new_unit_vec_z = new_unit_vec_x.cross(new_unit_vec_y);
        
        // now set the transformation matrix T_mat
        (*T_mat)(0,0) = new_unit_vec_x(0); (*T_mat)(0,3) = new_unit_vec_x(1); (*T_mat)(0,6) = new_unit_vec_x(2);
        (*T_mat)(1,1) = new_unit_vec_x(0); (*T_mat)(1,4) = new_unit_vec_x(1); (*T_mat)(1,7) = new_unit_vec_x(2);
        (*T_mat)(2,2) = new_unit_vec_x(0); (*T_mat)(2,5) = new_unit_vec_x(1); (*T_mat)(2,8) = new_unit_vec_x(2);
        
        (*T_mat)(3,0) = new_unit_vec_y(0); (*T_mat)(3,3) = new_unit_vec_y(1); (*T_mat)(3,6) = new_unit_vec_y(2);
        (*T_mat)(4,1) = new_unit_vec_y(0); (*T_mat)(4,4) = new_unit_vec_y(1); (*T_mat)(4,7) = new_unit_vec_y(2);
        (*T_mat)(5,2) = new_unit_vec_y(0); (*T_mat)(5,5) = new_unit_vec_y(1); (*T_mat)(5,8) = new_unit_vec_y(2);
        
        (*T_mat)(6,0) = new_unit_vec_z(0); (*T_mat)(6,3) = new_unit_vec_z(1); (*T_mat)(6,6) = new_unit_vec_z(2);
        (*T_mat)(7,1) = new_unit_vec_z(0); (*T_mat)(7,4) = new_unit_vec_z(1); (*T_mat)(7,7) = new_unit_vec_z(2);
        (*T_mat)(8,2) = new_unit_vec_z(0); (*T_mat)(8,5) = new_unit_vec_z(1); (*T_mat)(8,8) = new_unit_vec_z(2);
        
      }
      break;
      
    case QUAD4:
      {
        //get the four nodes of this elem
        const Point& point0 = local_elem->point(0);
        const Point& point1 = local_elem->point(1);
        const Point& point2 = local_elem->point(2);
        const Point& point3 = local_elem->point(3);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, vector03, tmp_vector1, tmp_vector2,
          new_unit_vec_x, new_unit_vec_y, new_unit_vec_z;
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        vector03.assign(point3 - point0);
        
        new_unit_vec_x = vector01.unit();
        
        tmp_vector1.assign(new_unit_vec_x.cross(vector03));
        tmp_vector2 = tmp_vector1.unit();
        
        new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
        new_unit_vec_z = new_unit_vec_x.cross(new_unit_vec_y);
        
        (*T_mat)(0,0) = new_unit_vec_x(0); (*T_mat)(0,4) = new_unit_vec_x(1); (*T_mat)(0,8) = new_unit_vec_x(2);
        (*T_mat)(1,1) = new_unit_vec_x(0); (*T_mat)(1,5) = new_unit_vec_x(1); (*T_mat)(1,9) = new_unit_vec_x(2);
        (*T_mat)(2,2) = new_unit_vec_x(0); (*T_mat)(2,6) = new_unit_vec_x(1); (*T_mat)(2,10) = new_unit_vec_x(2);
        (*T_mat)(3,3) = new_unit_vec_x(0); (*T_mat)(3,7) = new_unit_vec_x(1); (*T_mat)(3,11) = new_unit_vec_x(2);
        
        (*T_mat)(4,0) = new_unit_vec_y(0); (*T_mat)(4,4) = new_unit_vec_y(1); (*T_mat)(4,8) = new_unit_vec_y(2);
        (*T_mat)(5,1) = new_unit_vec_y(0); (*T_mat)(5,5) = new_unit_vec_y(1); (*T_mat)(5,9) = new_unit_vec_y(2);
        (*T_mat)(6,2) = new_unit_vec_y(0); (*T_mat)(6,6) = new_unit_vec_y(1); (*T_mat)(6,10) = new_unit_vec_y(2);
        (*T_mat)(7,3) = new_unit_vec_y(0); (*T_mat)(7,7) = new_unit_vec_y(1); (*T_mat)(7,11) = new_unit_vec_y(2);
        
        (*T_mat)(8,0) = new_unit_vec_z(0); (*T_mat)(8,4) = new_unit_vec_z(1); (*T_mat)(8,8) = new_unit_vec_z(2);
        (*T_mat)(9,1) = new_unit_vec_z(0); (*T_mat)(9,5) = new_unit_vec_z(1); (*T_mat)(9,9) = new_unit_vec_z(2);
        (*T_mat)(10,2) = new_unit_vec_z(0); (*T_mat)(10,6) = new_unit_vec_z(1); (*T_mat)(10,10) = new_unit_vec_z(2);
        (*T_mat)(11,3) = new_unit_vec_z(0); (*T_mat)(11,7) = new_unit_vec_z(1); (*T_mat)(11,11) = new_unit_vec_z(2);
        
      }
      break;
			
    default:
      //nothing to be done if not a 1-D OR 2-D
      break;
			
    }
}








void 
FESystemElem::FESystemElemBase::calculateShapeFunctionFactors
(DenseMatrix<double>* matrix,
 const unsigned int qty_name,
 const unsigned int design_point_enum_ID,
 const unsigned int domain,
 const FEFamily& family)
{
  
  assert (matrix != NULL);
  // also, make sure that the family exists in the map. It is assumed that 
  // if this family exists in one map, then it should exist in all maps in this elem.
  Assert(this->fe_base_map_for_DV[design_point_enum_ID].find(family) 
         != this->fe_base_map_for_DV[design_point_enum_ID].end(),
         ExcInternalError());

  
  
  if (!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID])
    this->initialize_element(design_point_enum_ID);

  // this needs to be called after the function has been initialized
  static unsigned int n_shape_funcs; 
  n_shape_funcs = 
    this->fe_base_map_for_DV[design_point_enum_ID][family]->n_shape_functions();  
  
  if (matrix->m() + matrix->n() != 2 * n_shape_funcs)
    matrix->resize(n_shape_funcs, n_shape_funcs);
  
  matrix->zero();

  
  // phi is the quantity that will store the
  static FEBase *fe_base_local;
  static QBase *qbase_local;
  
  fe_base_local = NULL;
  qbase_local = NULL;
  
  static Elem* local_elem;
  
  local_elem = this->local_elem_for_DV_map[design_point_enum_ID];
  
  
  // initialize the element if it is not already initialized
  switch(domain) 
    {
    case ELEM_VOLUME_ENUM_ID:
      {
        // now init the fe and quadrature pointers
        switch (local_elem->dim())
          {
          case 1:
            {
              // nothing to be done here, since the quantities have 
              // been calculated explicitly
            }
            break;
            
          case 2:
          case 3:
            {
              fe_base_local = this->fe_base_map_for_DV[design_point_enum_ID][family];
              qbase_local = this->qbase_map_for_DV[design_point_enum_ID][family];
            }
            break;
            
          default:
            abort();
            break;
          }
        
      }
      break;
      
    case SIDE_ZERO_ENUM_ID:
    case SIDE_ONE_ENUM_ID:
    case SIDE_TWO_ENUM_ID:
    case SIDE_THREE_ENUM_ID:
    case SIDE_FOUR_ENUM_ID:
    case SIDE_FIVE_ENUM_ID:
      {
        unsigned int elem_side = 
        this->getSideNumberFromDomainEnum(domain);
        
        // now init the fe and quadrature pointers
        switch (local_elem->dim())
          {
          case 1:
            {
              // this does not need to be set, since it 
              // can be explicitly evaluated for a 1-D element
            }
            break;
            
          case 2:
          case 3:
            {
              this->initialize_element_side(elem_side,
                                            design_point_enum_ID);
              fe_base_local = this->fe_side_map_for_DV[design_point_enum_ID][family];
              qbase_local = this->qbase_side_map_for_DV[design_point_enum_ID][family];
            }
            break;
            
          default:
            abort();
            break;
          }
      }
      break;
      
    default:
      abort();
      break;
    }
  
  
  // phi_1 and phi_2 are the quantities that will store the
  // the shape functions or their derivatives for calculation purposes
  const std::vector<std::vector<Real> > *phi_1 = NULL, *phi_2 = NULL;
  const std::vector<Real>* JxW = NULL;
  unsigned int n_quad_points= 0;
  
  
  switch(local_elem->dim())
    {
    case 1:
      {
        const Point& point0 = local_elem->point(0);
        const Point& point1 = local_elem->point(1);
        
        static Point vec;
        vec = point1;
        vec -= point0;
        double length = vec.size();
        
        // set the entries in the matrix and return
        switch (qty_name)
          {
          case N_X_N_X_FACTOR_ENUM_ID:
            {
              switch (local_elem->type())
                {
                case EDGE2:
                  {
                    (*matrix)(0,0) = 1.0;   (*matrix)(0,1) = -1.0;
                    (*matrix)(1,0) = -1.0;   (*matrix)(1,1) = 1.0;
                    
                    switch(domain)
                      {
                      case ELEM_VOLUME_ENUM_ID:
                        {
                          matrix->scale(1.0 / length);
                        }
                        break;
                        
                      case SIDE_ZERO_ENUM_ID:
                      case SIDE_ONE_ENUM_ID:
                        {
                          matrix->scale(1.0  / length / length);
                        }
                        break;
                        
                      default:
                        abort();
                        break;
                      }
                  }
                  break;
                  
                case EDGE3:
                  {
                    // has not been implemented yet.
                    abort();
                  }
                  break;
                  
                default:
                  abort();
                  break;
                }
            }
            break;
            
          case N_X_N_FACTOR_ENUM_ID:
            {
              switch (local_elem->type())
                {
                case EDGE2:
                  {
                    switch(domain)
                      {
                      case ELEM_VOLUME_ENUM_ID:
                        {
                          (*matrix)(0,0) = -0.5;   (*matrix)(0,1) = -0.5;
                          (*matrix)(1,0) = 0.5;   (*matrix)(1,1) = 0.5;
                        }
                        break;
                        
                      case SIDE_ZERO_ENUM_ID:
                        { 
                          (*matrix)(0,0) = -1.0;   (*matrix)(0,1) = 0.0;
                          (*matrix)(1,0) = 1.0;   (*matrix)(1,1) = 0.0;
                          matrix->scale(1.0 / length);
                        }
                        break;
                        
                      case SIDE_ONE_ENUM_ID:
                        { 
                          (*matrix)(0,0) = 0.0;   (*matrix)(0,1) = -1.0;
                          (*matrix)(1,0) = 0.0;   (*matrix)(1,1) = 1.0;
                          matrix->scale(1.0 / length);
                        }
                        break;
                        
                      default:
                        abort();
                        break;
                      }
                  }
                  break;
                  
                case EDGE3:
                  {
                    // has not been implemented yet.
                    abort();
                  }
                  break;
                  
                default:
                  abort();
                  break;
                }
            }
            break;
            
          case N_N_FACTOR_ENUM_ID:
            {
              switch (local_elem->type())
                {
                case EDGE2:
                  {
                    switch(domain)
                      {
                      case ELEM_VOLUME_ENUM_ID:
                        {
                          (*matrix)(0,0) = 2.0/3.0; (*matrix)(0,1) = 1.0/3.0;
                          (*matrix)(1,0) = 1.0/3.0; (*matrix)(1,1) = 2.0/3.0;
                          matrix->scale(0.5 * length);
                        }
                        break;
                        
                      case SIDE_ZERO_ENUM_ID:
                        {
                          (*matrix)(0,0) = 1.0;
                        }
                        break;
                        
                      case SIDE_ONE_ENUM_ID:
                        {
                          (*matrix)(1,1) = 1.0;
                        }
                        break;
                        
                      default:
                        abort();
                        break;
                      }
                  }
                  break;
                  
                case EDGE3:
                  {
                    // has not been implemented yet.
                    abort();
                  }
                  break;
                  
                default:
                  abort();
                  break;
                }
            }
            break;
            
          case N_Y_N_Y_FACTOR_ENUM_ID:
          case N_Z_N_Z_FACTOR_ENUM_ID:
          case N_X_N_Y_FACTOR_ENUM_ID:
          case N_Y_N_Z_FACTOR_ENUM_ID:
          case N_Z_N_X_FACTOR_ENUM_ID:
          case N_Y_N_FACTOR_ENUM_ID:
          case N_Z_N_FACTOR_ENUM_ID:
          default:
            abort();
            break;
          }
        // now we return since for 1-D, nothing else needs to be done.
        return;
      }
      break;
      
    case 2:
    case 3:
      {
        // phi_1 and phi_2 are the quantities that will store the
        // the shape functions or their derivatives for calculation purposes
        JxW = &(fe_base_local->get_JxW());
        n_quad_points= qbase_local->n_points();
        
        switch (qty_name)
          {
          case N_X_N_X_FACTOR_ENUM_ID:
            {
              phi_1 = &(fe_base_local->get_dphidx());
              phi_2 = &(fe_base_local->get_dphidx());
            }
            break;
            
          case N_Y_N_Y_FACTOR_ENUM_ID:
            {
              phi_1 =  &(fe_base_local->get_dphidy());
              phi_2 =  &(fe_base_local->get_dphidy());
            }
            break;
            
          case N_Z_N_Z_FACTOR_ENUM_ID:
            {
              phi_1 =  &(fe_base_local->get_dphidz());
              phi_2 =  &(fe_base_local->get_dphidz());
            }
            break;
            
          case N_X_N_Y_FACTOR_ENUM_ID:
            {
              phi_1 =  &(fe_base_local->get_dphidx());
              phi_2 =  &(fe_base_local->get_dphidy());
            }
            break;
            
          case N_Y_N_Z_FACTOR_ENUM_ID:
            {
              phi_1 =  &(fe_base_local->get_dphidy());
              phi_2 =  &(fe_base_local->get_dphidz());
            }
            break;
            
          case N_Z_N_X_FACTOR_ENUM_ID:
            {
              phi_1 =  &(fe_base_local->get_dphidz());
              phi_2 =  &(fe_base_local->get_dphidx());
            }
            break;
            
          case N_X_N_FACTOR_ENUM_ID:
            {
              phi_1 = &(fe_base_local->get_dphidx());
              phi_2 = &(fe_base_local->get_phi());
            }
            break;
            
          case N_Y_N_FACTOR_ENUM_ID:
            {
              phi_1 = &(fe_base_local->get_dphidy());
              phi_2 = &(fe_base_local->get_phi());
            }
            break;
            
          case N_Z_N_FACTOR_ENUM_ID:
            {
              phi_1 = &(fe_base_local->get_dphidz());
              phi_2 = &(fe_base_local->get_phi());
            }
            break;
            
          case N_N_FACTOR_ENUM_ID:
            {
              phi_1 = &(fe_base_local->get_phi());
              phi_2 = &(fe_base_local->get_phi());
            }
            break;
            
          default:
            abort();
            break;
          }
      }
      break;
			
    default:
      abort();
      break;
    }
  
  for (unsigned int qp=0; qp<n_quad_points; qp++)
    for (unsigned int i=0; i<n_shape_funcs; i++)
      for (unsigned int j=0; j<n_shape_funcs; j++)
        (*matrix)(i,j) += (*JxW)[qp] * (*phi_1)[i][qp] * (*phi_2)[j][qp];
        
}




void 
FESystemElem::FESystemElemBase::calculateShapeFunctionFactors
(DenseVector<double>* vector,
 const unsigned int qty_name,
 const unsigned int design_point_enum_ID,
 const unsigned int domain,
 const FEFamily& family)
{
  
  assert (vector != NULL);

  // also, make sure that the family exists in the map. It is assumed that 
  // if this family exists in one map, then it should exist in all maps in this elem.
  Assert(this->fe_base_map_for_DV[design_point_enum_ID].find(family)
         != this->fe_base_map_for_DV[design_point_enum_ID].end(),
         ExcInternalError());

  
  
 if (!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID])
   this->initialize_element(design_point_enum_ID);
  
  // this needs to be called after the function has been initialized
  static unsigned int n_shape_funcs; 
  n_shape_funcs = 
    this->fe_base_map_for_DV[design_point_enum_ID][family]->n_shape_functions();  
  
  if (vector->size() != n_shape_funcs)
    vector->resize(n_shape_funcs);
  
	vector->zero();
  
  // phi is the quantity that will store the
  static FEBase *fe_base_local;
  static QBase *qbase_local;

  fe_base_local = NULL;
  qbase_local = NULL;

  static Elem* local_elem;
  
  local_elem = this->local_elem_for_DV_map[design_point_enum_ID];
  
  // initialize the element if it is not already initialized
  switch(domain) 
    {
    case ELEM_VOLUME_ENUM_ID:
      {
        // now init the fe and quadrature pointers
        switch (local_elem->dim())
          {
          case 1:
            {
              // nothing to be done here, since the quantities have 
              // been calculated explicitly
            }
            break;
            
          case 2:
          case 3:
            {
              fe_base_local = this->fe_base_map_for_DV[design_point_enum_ID][family];
              qbase_local = this->qbase_map_for_DV[design_point_enum_ID][family];
            }
            break;
            
          default:
            abort();
            break;
          }
        
      }
      break;
      
    case SIDE_ZERO_ENUM_ID:
    case SIDE_ONE_ENUM_ID:
    case SIDE_TWO_ENUM_ID:
    case SIDE_THREE_ENUM_ID:
    case SIDE_FOUR_ENUM_ID:
    case SIDE_FIVE_ENUM_ID:
      {
        unsigned int elem_side = 
        this->getSideNumberFromDomainEnum(domain);
        
        // now init the fe and quadrature pointers
        switch (local_elem->dim())
          {
          case 1:
            {
              // this does not need to be set, since it 
              // can be explicitly evaluated for a 1-D element
            }
            break;
            
          case 2:
          case 3:
            {
              this->initialize_element_side(elem_side,
                                            design_point_enum_ID);
              fe_base_local = this->fe_side_map_for_DV[design_point_enum_ID][family];
              qbase_local = this->qbase_side_map_for_DV[design_point_enum_ID][family];
            }
            break;
            
          default:
            abort();
            break;
          }
      }
      break;
      
    default:
      abort();
      break;
    }
  
  
  const std::vector<std::vector<Real> > *phi = NULL;
  const std::vector<Real>* JxW = NULL;
  unsigned int n_quad_points= 0;
  
  switch(local_elem->dim())
    {
    case 1:
      {
        const Point& point0 = local_elem->point(0);
        const Point& point1 = local_elem->point(1);
        
        static Point vec;
        vec = point1;
        vec -= point0;
        double length = vec.size();
        
        // set the entries in the matrix and return
        switch (qty_name)
          {
          case N_FACTOR_ENUM_ID:
            {
              switch (local_elem->type())
                {
                case EDGE2:
                  {
                    switch(domain)
                      {
                      case ELEM_VOLUME_ENUM_ID:
                        {
                          (*vector)(0) = 0.5 * length; 
                          (*vector)(1) = 0.5 * length; 
                        }
                        break;
                        
                      case SIDE_ZERO_ENUM_ID:
                        {
                          (*vector)(0) = 1.0;
                        }
                        break;
                        
                      case SIDE_ONE_ENUM_ID:
                        {
                          (*vector)(1) = 1.0;
                        }
                        break;
                        
                      default:
                        abort();
                        break;
                      }
                  }
                  break;
                  
                case EDGE3:
                  {
                    // has not been implemented yet.
                    abort();
                  }
                  break;
                  
                default:
                  abort();
                  break;
                }
            }
            break;
            
            
          default:
            abort();
            break;
          }
        // now we return since for 1-D, nothing else needs to be done.
        return;
      }
      break;
      
    case 2:
    case 3:
      {
        JxW = &(fe_base_local->get_JxW());
        n_quad_points= qbase_local->n_points();
        
        switch (qty_name)
          {
          case N_FACTOR_ENUM_ID:
            {
              phi = &(fe_base_local->get_phi());
            }
            break;
            
          default:
            abort();
            break;
          }
      }
      break;
			
    default:
      abort();
      break;
    }
  
  
  for (unsigned int qp=0; qp<n_quad_points; qp++)
    for (unsigned int i=0; i<n_shape_funcs; i++)
      {
      (*(vector))(i) += (*JxW)[qp] * (*phi)[i][qp];
      }
						
}






void
FESystemElem::FESystemElemBase::calculateSurfaceNormal(DenseVector<double>* vector, 
                                                       const unsigned int design_point_enum_ID,
                                                       const unsigned int domain)
{
  // the normal will be calculated for a side of the local_elem
  // get a pointer to the local elem, and 
  // initialize the element if it is not already initialized
  if (!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID])
      this->initialize_element(design_point_enum_ID);
  
  // this needs to be called after the function has been initialized
  switch (vector->size())
    {
    case 3:
      {
        // keep going
      }
      break;
      
    default:
      vector->resize(3);
      break;
    }
  
	vector->zero();
  
  static bool load_on_elem_face;
  load_on_elem_face = false;
    
  unsigned int side_num = 
    this->getSideNumberFromDomainEnum(domain);
  
  static Point normal;
  normal.zero();

  static Elem* local_elem;
  local_elem = this->local_elem_for_DV_map[design_point_enum_ID];

  Assert(side_num <= local_elem->n_sides(),
         ExcInternalError());
  
  switch (local_elem->n_sides() - side_num)
    {
    case 0:
      // then the load is on the element face of a quad element
      load_on_elem_face = true;
      break;
      
    default:
      {
        // nothing to be done, since the load is not on elem face
      }
    }
  
  
  switch (local_elem->dim())
    {
    case 1:
      {
        // calculate the vector along the length of the element from 0->1
        const Point& point0 = local_elem->point(0);
        const Point& point1 = local_elem->point(1);
        
        static Point vector01;
        vector01.assign(point1 - point0);
        normal.assign(vector01.unit());
        
        switch (side_num)
          {
          case 0:
            {
              // scale the vector by -1.0
              normal *= -1.0;
            }
            break;
            
          case 1:
            break;
            
          default:
            libmesh_error();
            break;
          }
      }
      break;
			
    case 2:
      {
        switch (load_on_elem_face)
          {
          case true:
            {
              // calculate the vector along 0->1 and 0->2. The cross product of these two will give the
              // surface normal
              const Point& point0 = local_elem->point(0);
              const Point& point1 = local_elem->point(1);
              const Point& point2 = local_elem->point(2);
              
              //calculate the position vectors from point 0
              static Point vector01, vector02, tmp_vector;
              vector01.assign(point1 - point0);
              vector02.assign(point2 - point0);
              
              
              // calculate the cross product and scale it by the side (for +ve or -ve)
              tmp_vector.assign(vector01.cross(vector02));
              
              normal.assign(tmp_vector.unit());
            }
            break;
            
          case false:
          default:
            {
              // calculate the vector along the length of this side
              // the surface normal will be the cross-product of this vector with the z-vector
              // calculate the vector along the length of the element from 0->1
              std::auto_ptr<Elem> side_elem(local_elem->build_side(side_num).release());
              
              const Point& point0 = side_elem->point(0);
              const Point& point1 = side_elem->point(1);
              
              static Point vector01, z_vector, tmp_vector;
              z_vector.zero();
              z_vector(2) =  1.0;
              vector01.assign(point1 - point0);
              tmp_vector.assign(vector01.cross(z_vector));
              normal.assign(tmp_vector.unit());
            }
            break;
          }
      }
      break;
      
    case 3:
      {
        // get the side of this element, and follow the same procedure as the one for load_on_face
        std::auto_ptr<Elem> side_elem(local_elem->build_side(side_num).release());
        
        const Point& point0 = side_elem->point(0);
        const Point& point1 = side_elem->point(1);
        const Point& point2 = side_elem->point(2);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, tmp_vector;
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        
        
        // calculate the cross product and scale it by the side (for +ve or -ve)
        tmp_vector.assign(vector01.cross(vector02));
        
        normal.assign(tmp_vector.unit());
      }
      break;
    }
  
  // now copy the normal vector components
  for (unsigned int i=0; i<3; i++)
    (*vector)(i)  = normal(i);
}





void
FESystemElem::FESystemElemBase::calculateFESystemElemQty(DenseMatrix<double>* qty, 
                                                         const unsigned int name, 
                                                         const unsigned int design_point_enum_ID,
                                                         const unsigned int domain,
                                                         const FEFamily& family)
{
  assert (qty != NULL);
  
  switch(name)
    {
    case  TRANSFORM_MATRIX_ENUM_ID:
      {
        this->calculate_T_matrix(qty, design_point_enum_ID);
      }
      break;
			
      
    case  N_N_FACTOR_ENUM_ID:
    case  N_X_N_X_FACTOR_ENUM_ID:
    case  N_Y_N_Y_FACTOR_ENUM_ID:
    case  N_Z_N_Z_FACTOR_ENUM_ID:
    case  N_X_N_Y_FACTOR_ENUM_ID:
    case  N_Y_N_Z_FACTOR_ENUM_ID:
    case  N_Z_N_X_FACTOR_ENUM_ID:
    case  N_X_N_FACTOR_ENUM_ID:
    case  N_Y_N_FACTOR_ENUM_ID:
    case  N_Z_N_FACTOR_ENUM_ID:
      {
        this->calculateShapeFunctionFactors(qty, name, design_point_enum_ID, domain, family);
      }
      break;
      
      
    default:
      abort();
      break;
    }
  
}



void
FESystemElem::FESystemElemBase::calculateFESystemElemQty(DenseVector<double>* quantity,
                                                         const unsigned int qty_name,
                                                         const unsigned int design_point_enum_ID,
                                                         const unsigned int domain,
                                                         const FEFamily& family)
{
  assert (quantity != NULL);
  
  switch(qty_name)
    {
    case  SURFACE_NORMAL_ENUM_ID:
      {
        this->calculateSurfaceNormal(quantity, design_point_enum_ID, domain);
      }
      break;
      
    case N_FACTOR_ENUM_ID:
      {
        this->calculateShapeFunctionFactors(quantity, qty_name, design_point_enum_ID,
                                            domain, family);
      }
      break;
      
    default:
      abort();
      break;
    }
}







std::string 
FESystemElem::FESystemElemBase::getStringNameForFESystemElemQty
(const unsigned int quantity, 
 const unsigned int design_point_enum_ID,
 const unsigned int domain,
 const FEFamily& family)
{
  std::string name;

  // the name is used only for valid FE names, since otherwise the quantitiy 
  // is independent of the FE type.
  if (family != INVALID_FE)
    name = Utility::enum_to_string<FEFamily>(family);
  name += "_";
  name += FESystemElem::FESystemElemQtyEnum::enumName(quantity);
  name += "_";
  name += FESystemElem::DesignPointElemEnum::enumName(design_point_enum_ID);
  name += "_";
  name += FESystemElem::IntegrationDomainEnum::enumName(domain);
  
  return name;
}




std::string 
FESystemElem::FESystemElemBase::getStringNameForFESystemElemQtySensitivity
(const unsigned int quantity,
 unsigned int DV_num,
 const unsigned int domain,
 const FEFamily& family)
{
  // get the name for this qty
  std::string qty_name, qty_sensitivity_name;
  // the sensitivities are always calculated at the base element. Hence, the name 
  // will be referred to the base elem
  qty_name = this->getStringNameForFESystemElemQty(quantity, 
                                                   FESystemElem::BASE_ELEM::num(),
                                                   domain, family);
	
  // create the string and append the DV number
  qty_sensitivity_name = "d";
  qty_sensitivity_name += qty_name;
  qty_sensitivity_name += "_dDV";
	
  std::ostringstream DV_number_to_string; 
  DV_number_to_string << DV_num;
  qty_sensitivity_name += DV_number_to_string.str();
	
  return qty_sensitivity_name;
}





void
FESystemElem::FESystemElemBase::initializeFEAndQuadratureDataStructures()
{
  Assert(!this->fe_quadrature_data_initialized,
         ExcInvalidState());
  
  
  // get the FETypes for this element
  std::vector<FEType> fetypes;
  std::map<FEFamily, std::pair<QuadratureType, Order> > quadratures;
  
  this->getFETypes(fetypes);
  this->getQuadratureRules(quadratures);
  
  // now initialize the data
  // first of all, create and initialize all the finite elements needed for this element
  
  std::vector<FEType>::const_iterator it, end;
  it = fetypes.begin();
  end = fetypes.end();
  
  std::map<FEFamily, std::pair<QuadratureType, Order> >::const_iterator 
    qend = quadratures.end();
  
  std::auto_ptr<FEBase> fe;
  std::auto_ptr<QBase> qbase;
  bool insert_success = false;
  
  // create and insert the data structures necessary
  for (; it != end; it++)
    {
    // get the quadrature rule definintion for this elem type
    std::map<FEFamily, std::pair<QuadratureType, Order> >::const_iterator 
    qit = quadratures.find((*it).family);
    Assert(qit != qend, ExcInternalError());
    
    // now iterate over all the maps fe and qbase maps for all design points
    // and add the data structures
    std::map<unsigned int, std::map<FEFamily, FEBase*> >::iterator fe_it, fe_end;

    // insert the fe first
    fe_it = this->fe_base_map_for_DV.begin();
    fe_end = this->fe_base_map_for_DV.end();
    
    // iterate over all design points and add the data structs
    for ( ; fe_it != fe_end; fe_it++)
      {
      // add the fe
      fe.reset(FEBase::build(this->dimension, *it).release());
      
      insert_success = fe_it->second.insert(std::map<FEFamily, FEBase*>::value_type
                                            (fe->get_family(), fe.get())).second;
      Assert(insert_success, 
             FESystemExceptions::ExcDuplicateID("FEFamily", fe->get_family()));

      // now add the quadrature rule
      qbase.reset(QBase::build(qit->second.first, this->dimension, qit->second.second).release());
      
      // the insert for this is not checked since it will be true if it passed the previous test
      // it is assumed that there will be a set of qrules per design point. Hence, the size of 
      // the fe and the qrule map should be same.
      this->qbase_map_for_DV[fe_it->first].insert(std::map<FEFamily, QBase*>::value_type
                                                  (fe->get_family(), qbase.get()));
      
      // attach the quad rule
      fe->attach_quadrature_rule(qbase.get());
      
      // release the pointers
      fe.release();
      qbase.release();
      
      // this is repeated for the data structs for the side elems.  
      // add the fe
      fe.reset(FEBase::build(this->dimension, *it).release());
      
      this->fe_side_map_for_DV[fe_it->first].insert(std::map<FEFamily, FEBase*>::value_type
                                                    (fe->get_family(), fe.get()));
      // now add the quadrature rule
      qbase.reset(QBase::build(qit->second.first,
                               (this->dimension - 1), qit->second.second).release());
      
      // the insert for this is not checked since it will be true if it passed the previous test
      // it is assumed that there will be a set of qrules per design point. Hence, the size of 
      // the fe and the qrule map should be same.
      this->qbase_side_map_for_DV[fe_it->first].insert(std::map<FEFamily, QBase*>::value_type
                                                       (fe->get_family(), qbase.get()));
      
      // attach the quad rule
      fe->attach_quadrature_rule(qbase.get());
      
      // release the pointers
      fe.release();
      qbase.release();
      }
    
//    this->JxW_map[fe->get_family()] = 
//      const_cast<std::vector<Real>*>(&(fe->get_JxW()));
//    this->phi_map[fe->get_family()] = 
//      const_cast<std::vector<std::vector<Real> >*>(&(fe->get_phi()));
//    this->dphi_map[fe->get_family()] = 
//      const_cast<std::vector<std::vector<RealGradient> >*>(&(fe->get_dphi()));
    }
  
  this->fe_quadrature_data_initialized = true;
}


double FESystemElem::FESystemElemBase::getElementMass(bool sensitivity)
{
  // then, check about the property card
  if (!this->property_card_initialized)
    this->initPropertyCard();
  
  double factor = 0.0, factor_sens = 0.0, vol = 0.0;
  DenseVector<double> vector;
  
  // get the mass factor from the element data card
  switch (this->dimension)
  {
    case 1:
    case 3:
      this->elem_property_card->getFactor(factor, MASS_FACTOR::num());
      break;
      
    case 2:
      this->elem_property_card->getFactor(factor, MEMBRANE_MASS_FACTOR::num());
      break;
      
    default:
      ExcInternalError();
  }

  // get the value of volume
  switch (this->dimension)
  {
    case 1:
    {
      static Elem* elem;
      elem = this->geometric_elems_for_DV_map[FESystemElem::BASE_ELEM::num()];
      
      // calculate the length of the element
      const Point& point0 = elem->point(0);
      const Point& point1 = elem->point(1);
      
      static double dx, dy, dz;
      dx = point0(0)-point1(0);
      dy = point0(1)-point1(1);
      dz = point0(2)-point1(2);
      vol = pow((dx*dx+dy*dy+dz*dz),0.5);
    }
      break;
      
    case 2:
    case 3:
    {
      this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &vector, 
                               FESystemElem::BASE_ELEM::num(),
                               FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
      for (unsigned int i=0; i < vector.size(); i++)
        vol += vector.el(i);
    }
      
    default:
      ExcInternalError();
  }

  // check for property sensitivity
  if (sensitivity && 
      this->sensitivity_parameter == DesignData::PROPERTY_PARAMETER::num())
    {
      if (!this->elem_property_card->checkElemAndMaterialCardGlobalParameterDependence
          (this->sensitivity_parameter_ID))
        return 0.0;
      
      switch (this->dimension)
      {
        case 1:
        case 3:
          this->elem_property_card->getFactorSensitivityForGlobalParameter
          (factor_sens, MASS_FACTOR::num(), this->sensitivity_parameter_ID);
          break;
          
        case 2:
          this->elem_property_card->getFactorSensitivityForGlobalParameter
          (factor_sens, MEMBRANE_MASS_FACTOR::num(), this->sensitivity_parameter_ID);
          break;
          
        default:
          ExcInternalError();
      }
    }
  
  // to get the jacobian factor, get the N_FACTOR integrated over the volume
  // and then add all the elements of the vector
  if (sensitivity)
    {
      switch (this->sensitivity_parameter)
      {
        case PROPERTY_PARAMETER_ENUM_ID:
          vol *= factor_sens;
          break;
          
        case SHAPE_PARAMETER_ENUM_ID:
        {
          // this will set the sensitivity value of volume of element
          switch (this->dimension)
          {
            case 1:
            {
              static Elem* perturbed_elem;
              perturbed_elem = 
              this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()];
              // calculate the length for the perturbed element, and find the sensitivity
              // calculate the length of the element
              const Point& point0 = perturbed_elem->point(0);
              const Point& point1 = perturbed_elem->point(1);
              
              static double dx, dy, dz;
              dx = point0(0)-point1(0);
              dy = point0(1)-point1(1);
              dz = point0(2)-point1(2);
              vol -= pow((dx*dx+dy*dy+dz*dz),0.5);
              vol /= -0.5;
            }
              break;
              
            case 2:
            case 3:
            {
              this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_FACTOR::num(), &vector,
                                                       FESystemElem::ELEM_VOLUME::num(),
                                                       LAGRANGE);
              vol = 0.0;
              for (unsigned int i=0; i < vector.size(); i++)
                vol += vector.el(i);
            }
              
            default:
              ExcInternalError();
          }

          // multiply with the mass factor to get the mass
          vol *= factor;
        }
          break;
          
        default:
          abort();
          break;
      }
    }
  else
    vol *= factor;
  
  return vol;
}


std::auto_ptr<FESystemElem::FESystemElemBase> 
FESystemElem::createFESystemElem(const unsigned int elem_kind_enum_ID,
                                 Discipline::AnalysisDisciplineBase& discipline_base)
{
  std::auto_ptr<FESystemElem::FESystemElemBase> elem(NULL);
  
  switch (discipline_base.getDisciplineEnumID())
  {
    case PISTON_THEORY_ENUM_ID:
      switch (elem_kind_enum_ID)
    {
      case STRUCTURAL_PLATE_DKT_ENUM_ID:
      case STRUCTURAL_TRI3_VON_KARMAN_ENUM_ID:
      {
        elem.reset(new FESystemElem::PistonTheoryTri3(discipline_base)); 
        return elem;
      }
        break;
    }

    case STRUCTURAL_DISCIPLINE_ENUM_ID:
      switch (elem_kind_enum_ID)
    {
      case STRUCTURAL_BAR_EDGE2_ENUM_ID:
        elem.reset(new FESystemElem::Bar(discipline_base));
        break;
        
      case STRUCTURAL_BEAM_EDGE2_ENUM_ID:
        elem.reset(new FESystemElem::Beam2(discipline_base));
        break;
        
      case STRUCTURAL_LINEAR_SPRING_EDGE2_ENUM_ID:
        elem.reset(new FESystemElem::LinearSpring(discipline_base));  
        break;
        
      case STRUCTURAL_MEMBRANE_QUAD4_ENUM_ID:
        elem.reset(new FESystemElem::MembraneQuad4(discipline_base));  
        break;
        
      case STRUCTURAL_MEMBRANE_TRI3_ENUM_ID:
        elem.reset(new FESystemElem::MembraneTri3(discipline_base));  
        break;
        
      case STRUCTURAL_PLATE_MITC4_QUAD4_ENUM_ID:
        elem.reset(new FESystemElem::PlateMITC4(discipline_base));  
        break;
        
      case STRUCTURAL_PLATE_DKT_ENUM_ID:
        elem.reset(new FESystemElem::PlateDKT(discipline_base));
        break;
        
      case STRUCTURAL_TRI3_VON_KARMAN_ENUM_ID:
        elem.reset(new FESystemElem::Tri3VonKarman(discipline_base));
        break;

      default: 
        Assert(false, ExcInternalError());
        break;
    }
      break;
      
    case THERMAL_DISCIPLINE_ENUM_ID:
      switch (elem_kind_enum_ID)
    {
      case THERMAL_CONDUCTION_EDGE2_ENUM_ID:
        elem.reset(new FESystemElem::Conduction_1D(discipline_base));  
        break;
        
      case THERMAL_CONDUCTION_HEX8_ENUM_ID:
        elem.reset(new FESystemElem::Conduction_Hex8(discipline_base));  
        break;
        
      case THERMAL_CONDUCTION_PRISM6_ENUM_ID:
        elem.reset(new FESystemElem::Conduction_Prism6(discipline_base));  
        break;
        
      case THERMAL_CONDUCTION_QUAD4_ENUM_ID:
        elem.reset(new FESystemElem::Conduction_Quad4(discipline_base));  
        break;
        
      case THERMAL_CONDUCTION_TRI3_ENUM_ID:
        elem.reset(new FESystemElem::Conduction_Tri3(discipline_base));  
        break;
        
      case THERMAL_FORCED_CONVECTION_EDGE2_ENUM_ID:
        elem.reset(new FESystemElem::ForcedConvection_1D(discipline_base));  
        break;

      default: 
        Assert(false, ExcInternalError());
        break;
    }
      break;
      
    default: 
      Assert(false, ExcInternalError());
  }
    
  return elem;
}





