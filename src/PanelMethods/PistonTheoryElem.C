// $Id:$
/*
 *  PistonTheoryElem.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// FESystem includes
#include "PanelMethods/PistonTheoryElem.h"
#include "Discipline/PistonTheory.h"
#include "Loads/LoadCombination.h"


FESystemElem::PistonTheoryElem::PistonTheoryElem(const unsigned int dim, 
                                                 const unsigned int elem_enum_ID,
                                                 Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::FESystemElemBase(dim, elem_enum_ID, discipline),
fluid_flow_vector(3)
{
  
}



FESystemElem::PistonTheoryElem::~PistonTheoryElem()
{
  
}



void
FESystemElem::PistonTheoryElem::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType& fe = fetypes.back();
  fe.order = SECOND;
  fe.family = BCIZ;
  
}





void
FESystemElem::PistonTheoryElem::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
  quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                     value_type(BCIZ, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());
  
}




void
FESystemElem::PistonTheoryElem::getElementAssembledQty(const unsigned int name_enum_ID,
                                                       DenseMatrix<double>* qty)
{
  Assert(qty != NULL, ExcInternalError());
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static DenseMatrix<double> transformation_matrix(6*n_nodes, 6*n_nodes);
  
  if (!(qty->m() == 6*n_nodes && qty->n() == 6*n_nodes))
    {
      qty->resize(6*n_nodes, 6*n_nodes);
      transformation_matrix.resize(6*n_nodes, 6*n_nodes);
    }
  
  // get the transformation matrix
  transformation_matrix.zero();
  
  qty->zero();
  
  // calculate quantity and return it
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(), false);
  
  switch (this->dimension)
  {
    case 1:
    case 2:
    {
      this->getTransformationMatrix(&transformation_matrix, FESystemElem::BASE_ELEM::num(), false);
      
      // for a basic analysis, and property sensitivity, left multiply 
      // with the transformation matrix
      qty->left_multiply_transpose(transformation_matrix);
      qty->right_multiply(transformation_matrix);
    }
      break;
      
    default:
    {
      // keep going
    }
  }  
}



void
FESystemElem::PistonTheoryElem::getElementAssembledQty(const unsigned int name_enum_ID,
                                                       DenseVector<double>* qty)
{
  Assert(qty != NULL, ExcInternalError());
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static DenseVector<double> tmp_qty(6*n_nodes);
  static DenseMatrix<double> transformation_matrix;
  
  if (qty->size() != (6*n_nodes))
    {
      qty->resize(6*n_nodes);
      tmp_qty.resize(6*n_nodes);
      transformation_matrix.resize(6*n_nodes, 6*n_nodes);
    }
  
  // get the transformation matrix
  qty->zero();
  transformation_matrix.zero();
  tmp_qty.zero();
  
  // calculate quantity and return it
  this->calculateAssembledQty(&tmp_qty, name_enum_ID, FESystemElem::BASE_ELEM::num(), false);
  
  switch (this->dimension)
  {
    case 1:
    case 2:
    {
      this->getTransformationMatrix(&transformation_matrix, FESystemElem::BASE_ELEM::num(), false);
      transformation_matrix.left_multiply_vector(tmp_qty, *qty);
    }
      break;
      
    default:
    {
      // keep going
    }
  }
}



void
FESystemElem::PistonTheoryElem::getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                                  DenseMatrix<double>* qty)
{
  Assert(qty != NULL, ExcInternalError());
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static DenseMatrix<double> transformation_matrix(6*n_nodes, 6*n_nodes),
  transformation_matrix_sens(6*n_nodes, 6*n_nodes),
  basic_qty(6*n_nodes, 6*n_nodes);
  
  if (!(qty->m() == 6*n_nodes && qty->n() == 6*n_nodes))
    {
      qty->resize(6*n_nodes, 6*n_nodes);
      basic_qty.resize(6*n_nodes, 6*n_nodes);
      transformation_matrix.resize(6*n_nodes, 6*n_nodes);
      transformation_matrix_sens.resize(6*n_nodes, 6*n_nodes);
    }
  
  if (!this->property_card_initialized)
    this->initPropertyCard();
  
  // get the transformation matrix
  transformation_matrix.zero(); transformation_matrix_sens.zero();
  basic_qty.zero();
  
  qty->zero();
  
  // call the calculation routine
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(), true);
  
  switch (this->dimension)
  {
    case 1:
    case 2:
    {
      // get the transformation matrix    
      this->getTransformationMatrix(&transformation_matrix, 
                                    FESystemElem::BASE_ELEM::num(),false);
      
      // for a basic analysis, and property sensitivity, left multiply 
      // with the transformation matrix
      qty->left_multiply_transpose(transformation_matrix);
      qty->right_multiply(transformation_matrix);
      
      switch (this->sensitivity_parameter)
      {
        case SHAPE_PARAMETER_ENUM_ID:
        {
          // get the sensitivity of the T matrix, and the basic quantity       
          this->getTransformationMatrix(&transformation_matrix_sens,
                                        FESystemElem::BASE_ELEM::num(), true);
          this->calculateAssembledQty(&basic_qty, name_enum_ID, 
                                      FESystemElem::BASE_ELEM::num(), false);
          
          // the operation to be performed is d(transpose(T))/d_alpha * qty * T + 
          // transpose(T) * qty * d(T)/d_alpha. 
          basic_qty.left_multiply_transpose(transformation_matrix_sens);
          basic_qty.right_multiply(transformation_matrix);
          
          qty->add(1.0, basic_qty);
          qty->add_transpose(1.0, basic_qty);
        }
          break;
          
        default:
        {
          // keep going
        }
      }
    }
      break;
      
    default:
    {
      // keep going
    }
  }
}	   





void
FESystemElem::PistonTheoryElem::getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                                  DenseVector<double>* qty)
{
  Assert(qty != NULL, ExcInternalError());
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  static DenseMatrix<double> transformation_matrix(6*n_nodes, 6*n_nodes),
  transformation_matrix_sens(6*n_nodes, 6*n_nodes);
  static DenseVector<double> tmp_qty(6*n_nodes), basic_qty(6*n_nodes);
  
  if (qty->size() != (6*n_nodes))
    {
      qty->resize(6*n_nodes);
      tmp_qty.resize(6*n_nodes);
      basic_qty.resize(6*n_nodes);
      transformation_matrix.resize(6*n_nodes, 6*n_nodes);
      transformation_matrix_sens.resize(6*n_nodes, 6*n_nodes);
    }
  
  if (! this->property_card_initialized)
    this->initPropertyCard();
  
  // get the transformation matrix
  // create a quantity and resize it
  transformation_matrix.zero(); transformation_matrix_sens.zero();
  tmp_qty.zero(); basic_qty.zero();
  
  qty->zero();
  
  // call the calculation routine
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(), true);
  
  switch (this->dimension)
  {
    case 1:
    case 2:
    {
      // get the transformation matrix    
      this->getTransformationMatrix(&transformation_matrix, 
                                    FESystemElem::BASE_ELEM::num(),false);
      
      // for a basic analysis, and property sensitivity, left multiply 
      // with the transformation matrix
      transformation_matrix.left_multiply_vector(tmp_qty, *qty);
      
      switch (this->sensitivity_parameter)
      {
        case SHAPE_PARAMETER_ENUM_ID:
        {
          // get the sensitivity of the T matrix, and the basic quantity       
          this->getTransformationMatrix(&transformation_matrix_sens,
                                        FESystemElem::BASE_ELEM::num(), true);
          this->calculateAssembledQty(&basic_qty, name_enum_ID, 
                                      FESystemElem::BASE_ELEM::num(), false);
          
          tmp_qty.zero();
          transformation_matrix_sens.left_multiply_vector(basic_qty, tmp_qty);
          qty->add(1.0, tmp_qty);
        }
          break;
          
        default:
        {
          // keep going
        }
      }
    }
      break;
      
    default:
    {
      // keep going
    }
  }  
}	   






std::auto_ptr<ElemPostProcessQty> 
FESystemElem::PistonTheoryElem::getElementPostProcessQty(std::vector<unsigned int> load_cases, 
                                                         std::vector<DesignData::DesignParameter*> param)
{
  Assert(false, ExcInternalError());
}


void
FESystemElem::PistonTheoryElem::calculateAssembledQty(DenseMatrix<double>* qty,
                                                      const unsigned int qty_enum_ID,
                                                      const unsigned int design_point_enum_ID,
                                                      bool sensitivity_calculation)
{
  assert (qty != NULL);
  
  switch (qty_enum_ID)
  {
    case PISTON_THEORY_K_MATRIX_ENUM_ID:
      this->calculate_K(qty,design_point_enum_ID, sensitivity_calculation);
      break;
      
    case PISTON_THEORY_C1_MATRIX_ENUM_ID:
      this->calculate_C(qty,design_point_enum_ID, sensitivity_calculation);
      break;
      
    default:
      Assert(false, ExcInternalError());
      break;
  }
}



void 
FESystemElem::PistonTheoryElem::getTransformationMatrix(DenseMatrix<double>* matrix,
                                                        const unsigned int design_point_enum_ID,
                                                        bool shape_sensitivity)
{
  assert (matrix != NULL);
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  if (!(matrix->m() == 6*n_nodes && matrix->n() == 6*n_nodes))
    matrix->resize(n_nodes * 6, n_nodes * 6);
  
  matrix->zero();
  
  static DenseMatrix<double> T_mat;
  T_mat.zero();
  
  if (!shape_sensitivity)
    this->getFESystemElemQty(FESystemElem::TRANSFORM_MATRIX::num(), &T_mat, 
                             design_point_enum_ID, FESystemElem::ELEM_VOLUME::num(),
                             INVALID_FE);
  else
    {
      // make sure that the sensitivity is at the base elem
      Assert (design_point_enum_ID == FESystemElem::BASE_ELEM::num(), ExcInternalError());
      this->getFESystemElemQtyShapeSensitivity(FESystemElem::TRANSFORM_MATRIX::num(), &T_mat,
                                               FESystemElem::ELEM_VOLUME::num(), INVALID_FE);
    }
  
  for (unsigned k=0; k<2; k++)
    for (unsigned int i=0; i<(n_nodes*3); i++)
      for (unsigned int j =0; j<(n_nodes*3); j++)
        (*matrix)(i+(n_nodes*3*k),j+(n_nodes*3*k)) = T_mat(i,j);
}



void
FESystemElem::PistonTheoryElem::calculateAssembledQty(DenseVector<double>* qty,
                                                      const unsigned int qty_enum_ID,
                                                      const unsigned int design_point_enum_ID,
                                                      bool sensitivity_calculation)
{  
  // nothing to be calculated for this
  Assert(false, ExcInternalError());
}




void
FESystemElem::PistonTheoryElem::calculate_C(DenseMatrix<double>* matrix, 
                                            const unsigned int design_point_enum_ID,
                                            bool sensitivity_calculation)
{
  // nothing to be done for property sensitivity (property assumed to be materials, which are
  // not valid for this element kind)
  if (sensitivity_calculation && (this->sensitivity_parameter == PROPERTY_PARAMETER_ENUM_ID))
    return;

  unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static DenseMatrix<double> N_N(3*n_nodes, 3*n_nodes);
  N_N.zero(); 
  
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
  load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(PISTON_THEORY_SURFACE::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(PISTON_THEORY_SURFACE::num());
  
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();
  
  unsigned int side_num, n_sides;
  bool load_on_elem_face;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
  
  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
      side_num = load_iterator->second->getSurfaceID();
      
      // for 2-D elements, the element face itself can carry a surface flux load.
      load_on_elem_face = false;
      
      if (n_sides == side_num)
        load_on_elem_face = true;
            
      switch (this->dimension)
      {
        case 1:
          // piston theory not defined for 1-D elements
          Assert(false, ExcInternalError());
          break;
          
        case 2:
          // piston theory defined only on the face of the element, and not on the edge
          Assert(load_on_elem_face, ExcInternalError());
          break;
          
        case 3:
          // nothing specific to be done here. 
          break;
      }
      
      const unsigned int domain = this->getDomainEnumFromSideNumber(side_num);
            
      if (sensitivity_calculation)
        {
          // get the iterator for the load sensitivity
          load_sens_map_it = this->surface_load_sens->find(PISTON_THEORY_SURFACE::num());
          // if the load exists, its sensitivity should also be specified
          Assert(load_sens_map_it != this->surface_load_sens->end(), 
                 ExcInternalError());
          load_sens_iterator = load_sens_map_it->second.find(side_num);
          AssertThrow(load_sens_iterator != load_sens_map_it->second.end(),
                      ExcInternalError());
                    
          switch (this->sensitivity_parameter)
          {	      
            case SHAPE_PARAMETER_ENUM_ID:
              this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_N_FACTOR::num(), &N_N,
                                                       domain, BCIZ);
              
              break;
              
            case PROPERTY_PARAMETER_ENUM_ID:
            default:
              Assert(false, ExcInternalError());
              break;
          }
        }
      else
        this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N, 
                                 design_point_enum_ID, domain, BCIZ);
      
      // calculate the matrix
      // it is assumed here that rest of the code works fine, and so, 
      // m == n == n_nodes. Otherwise, it should be checked.
      for (unsigned int i=0; i<3*n_nodes; i++)
        for (unsigned int j=0; j<3*n_nodes; j++)
          (*matrix)(2*n_nodes+i,2*n_nodes+j) = N_N(i,j);
      matrix->scale(-1.0);

      // put a small diagonal number on dofs that do not participate
      for (unsigned int i=0; i<n_nodes;i++)
        {
          (*matrix)(i,i) = 1.0e-9;
          (*matrix)(n_nodes+i,n_nodes+i) = 1.0e-9;
          (*matrix)(5*n_nodes+i,5*n_nodes+i) = 1.0e-9;
        }
    }
}




void
FESystemElem::PistonTheoryElem::calculate_K(DenseMatrix<double>* matrix, 
                                            const unsigned int design_point_enum_ID,
                                            bool sensitivity_calculation)
{
  
  // nothing to be done for property sensitivity (property assumed to be materials, which are
  // not valid for this element kind)
  if (sensitivity_calculation && (this->sensitivity_parameter == PROPERTY_PARAMETER_ENUM_ID))
    return;
  
  unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static DenseMatrix<double> Nx_N(3*n_nodes, 3*n_nodes), Ny_N(3*n_nodes, 3*n_nodes),
  Nz_N(3*n_nodes, 3*n_nodes), transformation_matrix(n_nodes, n_nodes);
  Nx_N.zero(); Ny_N.zero(); Nz_N.zero(), transformation_matrix.zero(); 
  static DenseVector<double> fluid_vec(3);
  fluid_vec.zero();
  
  
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
  load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(PISTON_THEORY_SURFACE::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(PISTON_THEORY_SURFACE::num());
  
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();
  
  unsigned int side_num, n_sides;
  bool load_on_elem_face;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
  
  
  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
      side_num = load_iterator->second->getSurfaceID();
      
      // for 2-D elements, the element face itself can carry a surface flux load.
      load_on_elem_face = false;
      
      if (n_sides == side_num)
        load_on_elem_face = true;
      
      switch (this->dimension)
      {
        case 1:
          // piston theory not defined for 1-D elements
          Assert(false, ExcInternalError());
          break;
          
        case 2:
          // piston theory defined only on the face of the element, and not on the edge
          Assert(load_on_elem_face, ExcInternalError());
          break;
          
        case 3:
          // nothing specific to be done here. 
          break;
      }
      
      const unsigned int domain = this->getDomainEnumFromSideNumber(side_num);
      
      if (sensitivity_calculation)
        {
          // get the iterator for the load sensitivity
          load_sens_map_it = this->surface_load_sens->find(PISTON_THEORY_SURFACE::num());
          // if the load exists, its sensitivity should also be specified
          Assert(load_sens_map_it != this->surface_load_sens->end(), 
                 ExcInternalError());
          load_sens_iterator = load_sens_map_it->second.find(side_num);
          AssertThrow(load_sens_iterator != load_sens_map_it->second.end(),
                      ExcInternalError());
          
          switch (this->sensitivity_parameter)
          {	                    
            case SHAPE_PARAMETER_ENUM_ID:
            {
              Assert(false, ExcInternalError());
            }
              break;
              
            case PROPERTY_PARAMETER_ENUM_ID:
            default:
              Assert(false, ExcInternalError());
              break;
          }
        }
      else
        {
          this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &Nx_N, 
                                   design_point_enum_ID, domain, BCIZ);
          this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &Ny_N, 
                                   design_point_enum_ID, domain, BCIZ);
          this->getFESystemElemQty(FESystemElem::N_Z_N_FACTOR::num(), &Nz_N, 
                                   design_point_enum_ID, domain, BCIZ);

          // the fluid flow vector is copied into a vector of size that can be multiplied 
          // with the transformation matrix, which has a size of n_nodes*6
          DenseVector<double> vec1(n_nodes*6), vec2(n_nodes*6); 
          vec1.zero(); vec2.zero(); 
          for (unsigned int i=0; i<3; i++)
            vec1(i*3) = this->fluid_flow_vector(i);
            
            
          // get the transformation matrix
          this->getTransformationMatrix(&transformation_matrix, design_point_enum_ID, false);
          transformation_matrix.right_multiply_vector(vec1, vec2);
          
          for (unsigned int i=0; i<3; i++)
            fluid_vec(i) = vec2(i*3);
          
          Nx_N.scale(-fluid_vec(0));
          Ny_N.scale(-fluid_vec(1));
          Nz_N.scale(-fluid_vec(2));
        }
      
      // calculate the matrix
      // it is assumed here that rest of the code works fine, and so, 
      // m == n == n_nodes. Otherwise, it should be checked.
      for (unsigned int i=0; i<3*n_nodes; i++)
        for (unsigned int j=0; j<3*n_nodes; j++)
          (*matrix)(2*n_nodes+i,2*n_nodes+j) = Nx_N(j,i) + Ny_N(j,i) + Nz_N(j,i);
      
      // put a small diagonal number on dofs that do not participate
      for (unsigned int i=0; i<n_nodes;i++)
        {
          (*matrix)(i,i) = 1.0e-9;
          (*matrix)(n_nodes+i,n_nodes+i) = 1.0e-9;
          (*matrix)(5*n_nodes+i,5*n_nodes+i) = 1.0e-9;
        }

    }
}




void
FESystemElem::PistonTheoryElem::initPropertyCard()
{
  // nothing to be done here. 
}




void
FESystemElem::PistonTheoryElem::clearPropertyCardInitialization()
{
  // nothing to be done here
}





void
FESystemElem::PistonTheoryElem::setFluidFlowVector(const Point& vector)
{
  for (unsigned int i=0; i < 3; i++)
    this->fluid_flow_vector(i) = vector(i);
}





void
FESystemElem::PistonTheoryElem::clearLoadAndDofs()
{
  this->fluid_flow_vector.zero();
  FESystemElem::FESystemElemBase::clearLoadAndDofs();
}




