// $Id: FEInterpolationElem.C,v 1.6 2006-10-23 23:42:35 manav Exp $

// C++ includes
#include <sstream>
#include <memory>

// FESystem includes
#include "Interpolation/FEInterpolationElem.h"
#include "Database/ElementDataStorage.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "Interpolation/FEInterpolation.h"
#include "Discipline/AnalysisDisciplineBase.h"

// libMesh includes
#include "geom/face_quad4.h"
#include "geom/edge_edge2.h"
#include "fe/fe.h"




FEInterpolationElem::FEInterpolationElem(const unsigned int dim, 
                                         const unsigned int elem_enum_ID,
                                         Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::FESystemElemBase(dim, elem_enum_ID, discipline)
{
  
}



FEInterpolationElem::~FEInterpolationElem()
{
	
}






void
FEInterpolationElem::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType& fe = fetypes.back();
  fe.order = FIRST;
  fe.family = LAGRANGE;
}





void
FEInterpolationElem::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
    quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                       value_type(LAGRANGE, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());
}






void FEInterpolationElem::getElementAssembledQty(const unsigned int name,
                                                 DenseMatrix<double>* qty)
{
  assert (qty != NULL);
  
  unsigned int n_nodes = this->getNNodes();
  
  if (qty->m() != (n_nodes) || qty->n() != (n_nodes))
    qty->resize(n_nodes, n_nodes);
  
  // calculate quantity and return it
  this->calculateAssembledQty(qty, name, FESystemElem::BASE_ELEM::num(), false);
}



void FEInterpolationElem::getElementAssembledQty(const unsigned int name,
                                                 DenseVector<double>* qty)
{
  assert (qty != NULL);
  unsigned int n_nodes = this->getNNodes();
  
  if (qty->size() != (n_nodes))
    qty->resize(n_nodes);
  
  // calculate quantity and return it
  this->calculateAssembledQty(qty, name, FESystemElem::BASE_ELEM::num(), false);
}




void FEInterpolationElem::getElementAssembledQtySensitivity(const unsigned int name,
                                                   DenseMatrix<double>* qty)
{
  assert (qty != NULL);
  
  unsigned int n_nodes = this->getNNodes();
  
  if (qty->m() != (n_nodes) || qty->n() != (n_nodes))
    qty->resize(n_nodes, n_nodes);
  
  // calculate quantity and return it
  this->calculateAssembledQty(qty, name, FESystemElem::BASE_ELEM::num(), true);
}	   





void FEInterpolationElem::getElementAssembledQtySensitivity(const unsigned int name,
                                                            DenseVector<double> *qty)
{
  assert (qty != NULL);
  
  unsigned int n_nodes = this->getNNodes();
  
  if (qty->size() != (n_nodes))
    qty->resize(n_nodes);
  
  // calculate quantity and return it
  this->calculateAssembledQty(qty, name, FESystemElem::BASE_ELEM::num(), true);
} 





void FEInterpolationElem::calculateAssembledQty(DenseMatrix<double>* quantity,
                                                const unsigned int qty_name,
                                                const unsigned int design_point_enum_ID,
                                                bool sensitivity_calc)
{
  assert (quantity != NULL);
  
  switch (qty_name)
    {
    case FE_INTERPOLATION_ELEM_K_MATRIX_ENUM_ID:
      {
        return this->calculate_K(quantity, design_point_enum_ID, sensitivity_calc);
      }
      break;
			
    default:
      abort();
    }
}



void 
FEInterpolationElem::calculateAssembledQty(DenseVector<double>* quantity,
                                           const unsigned int qty_name,
                                           const unsigned int design_point_enum_ID,
                                           bool sensitivity_calc)
{
  assert (quantity != NULL);
	
  switch (qty_name)
    {
    case FE_INTERPOLATION_ELEM_F_VECTOR_ENUM_ID:
      {
        return this->calculate_F(quantity, design_point_enum_ID, sensitivity_calc);
      }
      break;
			
    default:
      abort();
    }
}






//std::string 
//FEInterpolationElem::getStringNameForFEInterpolationElemQty(const unsigned int quantity, 
//                                                            const unsigned int domain)
//{
//  std::string name = FEInterpolationElemQtyEnum::enumName(quantity);
//	
//  // now add the side number to the name
//  unsigned int side_num = this->getSideNumberFromDomainEnum(domain);
//  std::ostringstream side_number_to_string; 
//  side_number_to_string << side_num;
//  name += "_side_";
//  name += side_number_to_string.str();
//  
//  return name;
//}
//
//
//
//
//std::string 
//FEInterpolationElem::getStringNameForFEInterpolationElemQtySensitivity
//(const unsigned int quantity,
// const unsigned int DV_num,
// const unsigned int domain)
//{
//  // get the name for this qty
//  std::string qty_name, qty_sensitivity_name;
//  qty_name = this->getStringNameForFEInterpolationElemQty(quantity, domain);
//	
//  // create the string and append the DV number
//  qty_sensitivity_name = "d";
//  qty_sensitivity_name += qty_name;
//  qty_sensitivity_name += "_dDV";
//	
//  std::ostringstream DV_number_to_string; 
//  DV_number_to_string << DV_num;
//  qty_sensitivity_name += DV_number_to_string.str();
//	
//  return qty_sensitivity_name;
//}
//





void
FEInterpolationElem::calculate_K(DenseMatrix<double>* matrix, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation)
{
  if (sensitivity_calculation == true)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          matrix->zero();
        }
        break;
				
      case SHAPE_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(),
                                   matrix, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), 
                                   LAGRANGE);
        }
        break;
				
      default:
        abort();
        break;
      }
    }
  else 
    {
    this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), matrix, 
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), 
                             LAGRANGE);
    }
}





void
FEInterpolationElem::calculate_F(DenseVector<double>* f_vector, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation)
{
  // calculation of the element load will involve the following procedure
  // -- create a vector of global coordinates of the Gauss points 
  //    of this element.
  // -- find out the element on the 'from_mesh' inside which this element
  // belongs.
  // -- find out the value of the dof to be interpolated at the point
  // -- once this has been done for all gauss points, caluculate the force
  //    vector
  
  // initialize the element if it is not already initialized
  if (this->local_elem_is_initialized_for_DV_map[design_point_enum_ID] == false)
    this->initialize_element(design_point_enum_ID);
  
  // get the phi for this elem, and the nodes for the base elem
  // multiply the phi and node location and get the localtion of the
  // element
	
  const std::vector<std::vector<Real> >& phi = 
    this->fe_base_map_for_DV[design_point_enum_ID][LAGRANGE]->get_phi();
  const std::vector<Real>& JxW =
    this->fe_base_map_for_DV[design_point_enum_ID][LAGRANGE]->get_JxW();
  
  
  std::vector<Point> gauss_coords(phi[0].size());
  
  unsigned int n_quad_points = phi[0].size();
  unsigned int n_nodes = phi.size(); 
  
  double x=0.0, y=0.0, z=0.0;
  std::vector<Elem*> source_elems;
  unsigned int mesh_ID = this->analysis_discipline.getAnalysisMeshID();
  
  FEInterpolation* fe_interpolation = 
    dynamic_cast<FEInterpolation*> (&(this->analysis_discipline.getAnalysisDriver()));
  
  for (unsigned int j=0; j<n_quad_points; j++) // iterate over each gauss point
    {
    x = 0.0, y=0.0, z=0.0;
    for (unsigned int i=0; i< n_nodes; i++)  // iterate over each node shape function
      {
      gauss_coords[j].add_scaled(this->geometric_elems_for_DV_map[design_point_enum_ID]->point(i), phi[i][j]);
      }
    
    // for each point, ask the interpolation discipline to find the element 
    // which contains the point
    source_elems.push_back(fe_interpolation->getElemContainingPoint
                           (gauss_coords[j],
                            this->elem_ID,
                            mesh_ID));
    }
  
  
  // now, for each quadrature point, get the local coordinates in the source
  // elem, calculate the value of the variable at the point and add it to the 
  // force vector
  Point inverse_coords;
  std::auto_ptr<FEBase> fe(NULL);
  
  switch (this->dimension)
    {
    case 1:
      {
        fe.reset(new FE<1,LAGRANGE>(FEType()));
      }
      break;
      
    case 2:
      {
        fe.reset(new FE<2,LAGRANGE>(FEType()));
      }
      break;
      
    case 3:
      {
        fe.reset(new FE<3,LAGRANGE>(FEType()));
      }
      break;
      
    default:
      abort();
    }
  
  DenseVector<double> source_dof_values, qp_values(n_quad_points);
  std::vector<Point> point_vec;
  
  for (unsigned int qp =0; qp < n_quad_points; qp++)
    {
    // get the shape function values in the source elem
    // this has to be initialized as a 2-D element. This is because libMesh needs the 
    // elements of n dim to live in an n-dim space for the method inverse map.
    // this needs to be looked into in greater detail to solve the problem.
    inverse_coords = FE<2,LAGRANGE>::inverse_map(source_elems[qp], 
                                                 gauss_coords[qp]);
    // calculate the value of the source variable at this 
    // point
    if (source_dof_values.size() != (source_elems[qp])->n_nodes())
      source_dof_values.resize((source_elems[qp])->n_nodes());
    source_dof_values.zero();
    
    fe_interpolation->getDofValuesForSourceElem(source_elems[qp],
                                                source_dof_values);
    point_vec.clear();
    point_vec.push_back(inverse_coords);
    fe->reinit(source_elems[qp], &point_vec);
    
    const std::vector<std::vector<Real> >& source_phi = fe->get_phi();
    for (unsigned int i=0; i < source_phi.size(); i++)
      qp_values(qp) += source_phi[i][0] * source_dof_values(i);
    }
  
  
  
  for (unsigned int qp=0; qp < n_quad_points; qp++)
    for (unsigned int i=0; i < phi.size(); i++)
      (*f_vector)(i) += JxW[qp] * phi[i][qp] * qp_values(qp);
  
}



std::auto_ptr<ElemPostProcessQty>
FEInterpolationElem::getElementPostProcessQty(std::vector<unsigned int> load_cases,
                                              std::vector<DesignData::DesignParameter*> DV_vector)
{
  // for this element, this does not make sense, hence, this method cannot be called.
  // consequently, an abort statement has been placed here.
  std::auto_ptr<ElemPostProcessQty> return_qty(NULL);
  
  abort();
  return return_qty;
}
