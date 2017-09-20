// $Id: FluidElem.C,v 1.1.2.1 2008-02-25 04:25:45 manav Exp $

// C++ includes
#include <sstream>


// FESystem includes
#include "FluidElems/fluid_elem.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Discipline/FluidAnalysis.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Database/ElementDataStorage.h"
#include "Database/GlobalDataStorage.h"
#include "Loads/LoadDatabase.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "Loads/LoadCombination.h"
#include "Properties/MaterialPropertyNameEnums.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic2DElemDataCard.h"
#include "Properties/Isotropic1DElemDataCard.h"



FESystemElem::FluidElem::FluidElem
(const unsigned int dim, 
 const unsigned int elem_enum_ID, 
 Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::FESystemElemBase(dim, elem_enum_ID, discipline)
{
  
}



FESystemElem::FluidElem::~FluidElem()
{
	
}







void
FESystemElem::FluidElem::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType& fe = fetypes.back();
  fe.order = FIRST;
  fe.family = LAGRANGE;
}





void
FESystemElem::FluidElem::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
    quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                       value_type(LAGRANGE, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());
}





void
FESystemElem::FluidElem::calculateAssembledQty(DenseMatrix<double>* quantity,
                                                 const unsigned int qty_name,
                                                 const unsigned int design_point_enum_ID,
                                                 bool sensitivity_calc)
{
  assert (quantity != NULL);
  
  // make sure that if sensitivity is being requested, it is at the base point
  if(sensitivity_calc)
    Assert (design_point_enum_ID == FESystemElem::BASE_ELEM::num(), 
            ExcInternalError());
  
  // check if the local element has been initialized or not. If not, 
  // initialize it, also, resize the quantity
  if (!this->local_elem_is_initialized_for_DV_map.find(design_point_enum_ID)->second)
    this->initialize_element(design_point_enum_ID);
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  if (!(quantity->m() == n_nodes && quantity->n() == n_nodes))
    {
    quantity->resize(n_nodes, n_nodes);
    quantity->zero();
    }
  
  if (!this->property_card_initialized)
    this->initPropertyCard();
  
  switch (qty_name)
    {
    case FLUID_INVISCID_JACOBIAN_MATRIX_ENUM_ID:
        this->calculateInviscidJacobian(quantity, design_point_enum_ID, sensitivity_calc);
      break;
			
    case FLUID_VISCOUS_JACOBIAN_MATRIX_ENUM_ID:
        this->calculateViscousJacobian(quantity, design_point_enum_ID, sensitivity_calc);
      break;
			
    default:
      Assert(false, ExcInternalError());
    }
}



void
FESystemElem::FluidElem::calculateAssembledQty(DenseVector<double>* quantity,
                                                 const unsigned int qty_name,
                                                 const unsigned int design_point_enum_ID,
                                                 bool sensitivity_calc)
{
  assert (quantity != NULL);

  // make sure that if sensitivity is being requested, it is at the base point
  if(sensitivity_calc)
    Assert (design_point_enum_ID == FESystemElem::BASE_ELEM::num(), 
            ExcInternalError());
  
  // check if the local element has been initialized or not. If not, 
  // initialize it, also, resize the quantity
  if (!this->local_elem_is_initialized_for_DV_map.find(design_point_enum_ID)->second)
    this->initialize_element(design_point_enum_ID);
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  if (quantity->size() != n_nodes)
    {
    quantity->resize(n_nodes);
    quantity->zero();
    }
  
  if (!this->property_card_initialized)
    this->initPropertyCard();
  

  switch (qty_name)
    {
    case FLUID_F_VECTOR_ENUM_ID:
      this->calculate_F(quantity,	 design_point_enum_ID, sensitivity_calc);
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
}





void
FESystemElem::FluidElem::calculateInviscidJacobian(DenseMatrix<double>* matrix, 
                                                   const unsigned int design_point_enum_ID,
                                                   bool sensitivity_calculation)
{
  static double factor, factor_sens; 
  factor = 0.0, factor_sens = 0.0; 
  
  Assert(!sensitivity_calculation, ExcInternalError());
  
//  if (sensitivity_calculation && 
//      this->sensitivity_parameter == DesignData::PROPERTY_PARAMETER::num())
//    {
//    if (!this->elem_property_card->checkElemAndMaterialCardGlobalParameterDependence
//        (this->sensitivity_parameter_ID))
//      return;
//    
//    this->elem_property_card->getFactorSensitivityForGlobalParameter
//      (factor_sens, FLUID_CONDUCTANCE_FACTOR::num(), this->sensitivity_parameter_ID);
//    }
  
  static DenseMatrix<double> mat;
  
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          switch ( this->dimension)
            {
            case 3:
              {
                mat.zero();
                // get the factors and add them to the matrix dN/dx * dN/dx
                this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &mat, 
                                         design_point_enum_ID,
                                         FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 2:
              {
                // dN/dy * dN/dy
                mat.zero();
                this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &mat, 
                                         design_point_enum_ID,
                                         FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 1:
              {
                mat.zero();
                // get the factors and add them to the matrix dN/dx * dN/dx
                this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &mat, 
                                         design_point_enum_ID,
                                         FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              break;
              
            default:
              abort();
            }
          
          matrix->scale(factor_sens);
        }
        break;
				
      case SHAPE_PARAMETER_ENUM_ID:
        {
          switch ( this->dimension)
            {
            case 3:
              {
                // get the factors and add them to the matrix dN/dx * dN/dx
                mat.zero();
                this->getFESystemElemQtyShapeSensitivity
                  (FESystemElem::N_Z_N_Z_FACTOR::num(), &mat,
                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 2:
              {  
                // dN/dy * dN/dy
                mat.zero();
                this->getFESystemElemQtyShapeSensitivity
                  (FESystemElem::N_Y_N_Y_FACTOR::num(), &mat, 
                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 1:
              {
                // get the factors and add them to the matrix dN/dx * dN/dx
                mat.zero();
                this->getFESystemElemQtyShapeSensitivity
                  (FESystemElem::N_X_N_X_FACTOR::num(), &mat,
                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              break;
              
            default:
              abort();
            }
          
          matrix->scale(factor);
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else 
    {
    switch ( this->dimension)
      {
      case 3:
        {
          mat.zero();
          // get the factors and add them to the matrix dN/dx * dN/dx
          this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &mat, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          matrix->add(1.0, mat);
        }
        
      case 2:
        {
          // dN/dy * dN/dy
          mat.zero();
          this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &mat, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          matrix->add(1.0, mat);
        }
        
      case 1:
        {
          mat.zero();
          // get the factors and add them to the matrix dN/dx * dN/dx
          this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &mat, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          matrix->add(1.0, mat);
        }
        break;
        
      default:
        abort();
      }
    
    matrix->scale(factor);
    }
}




void
FESystemElem::FluidElem::calculate_C_Jac(DenseMatrix<double>* matrix,
                                           const unsigned int design_point_enum_ID)
{
  if (!this->analysis_discipline.checkPropertyDependenceOnTemperature())
    return;
  
  static double factor; 
  factor = 0.0; 
  
  this->elem_property_card->getFactorSensitivityForLocalParameter
    (factor, FLUID_CAPACITANCE_FACTOR::num(), Property::TEMPERATURE::num());
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  // it is assumed that the matrices and vectors remain the same size as 
  // created here. Hence, they are not resized again.
  static DenseMatrix<double> mat1(n_nodes, n_nodes);
  static DenseVector<double> vec2(n_nodes);
  if ((vec2.size() != n_nodes))
    {
    vec2.resize(n_nodes);
    mat1.resize(n_nodes, n_nodes); 
    }

  mat1.zero();
  vec2.zero();  
  
  this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &mat1, 
                           design_point_enum_ID,
                           FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
  mat1.scale(factor / (1.0 * n_nodes));
  
  mat1.right_multiply_vector(this->dof_values_vec[1], vec2);
  
  // now, this vector has to be placed inside the jacobian matrix
  static double value = 0.0;
  for (unsigned int i=0; i<n_nodes; i++)
    {
    value = vec2(i);
    for (unsigned int j=0; j < n_nodes; j++)
      (*matrix)(i,j) = value;
    }

  // diagonalize the matrix since the non-diagonal capacitance is not working well
  for (unsigned int i=0; i < matrix->m(); i++)
    for (unsigned int j=0; j < matrix->m(); j++)
      {
      if (i==j)
        continue;
      (*matrix)(i,i) += (*matrix)(i,j);
      (*matrix)(i,j) = 0.0;
      }
}



void
FESystemElem::FluidElem::calculate_K_c(DenseMatrix<double>* matrix, 
                                         const unsigned int design_point_enum_ID,
                                         bool sensitivity_calculation)
{
  static double factor, factor_sens; 
  factor = 0.0, factor_sens = 0.0; 
  
  this->elem_property_card->getFactor(factor, FLUID_CONDUCTANCE_FACTOR::num());
  if (sensitivity_calculation && 
      this->sensitivity_parameter == DesignData::PROPERTY_PARAMETER::num())
    {
    if (!this->elem_property_card->checkElemAndMaterialCardGlobalParameterDependence
        (this->sensitivity_parameter_ID))
      return;
    
    this->elem_property_card->getFactorSensitivityForGlobalParameter
      (factor_sens, FLUID_CONDUCTANCE_FACTOR::num(), this->sensitivity_parameter_ID);
    }
  
  static DenseMatrix<double> mat;

  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          switch ( this->dimension)
            {
            case 3:
              {
                mat.zero();
                // get the factors and add them to the matrix dN/dx * dN/dx
                this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &mat, 
                                         design_point_enum_ID,
                                         FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 2:
              {
                // dN/dy * dN/dy
                mat.zero();
                this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &mat, 
                                         design_point_enum_ID,
                                         FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
              case 1:
              {
                mat.zero();
                // get the factors and add them to the matrix dN/dx * dN/dx
                this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &mat, 
                                         design_point_enum_ID,
                                         FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              break;
             
              default:
                abort();
            }
          
          matrix->scale(factor_sens);
        }
        break;
				
      case SHAPE_PARAMETER_ENUM_ID:
        {
          switch ( this->dimension)
            {
            case 3:
              {
                // get the factors and add them to the matrix dN/dx * dN/dx
                mat.zero();
                this->getFESystemElemQtyShapeSensitivity
                  (FESystemElem::N_Z_N_Z_FACTOR::num(), &mat,
                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 2:
              {  
                // dN/dy * dN/dy
                mat.zero();
                this->getFESystemElemQtyShapeSensitivity
                  (FESystemElem::N_Y_N_Y_FACTOR::num(), &mat, 
                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              
            case 1:
              {
                // get the factors and add them to the matrix dN/dx * dN/dx
                mat.zero();
                this->getFESystemElemQtyShapeSensitivity
                  (FESystemElem::N_X_N_X_FACTOR::num(), &mat,
                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
                matrix->add(1.0, mat);
              }
              break;
              
            default:
              abort();
            }
          
          matrix->scale(factor);
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else 
    {
    switch ( this->dimension)
      {
      case 3:
        {
          mat.zero();
          // get the factors and add them to the matrix dN/dx * dN/dx
          this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &mat, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          matrix->add(1.0, mat);
        }
        
      case 2:
        {
          // dN/dy * dN/dy
          mat.zero();
          this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &mat, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          matrix->add(1.0, mat);
        }
        
      case 1:
        {
          mat.zero();
          // get the factors and add them to the matrix dN/dx * dN/dx
          this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &mat, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          matrix->add(1.0, mat);
        }
        break;
        
      default:
        abort();
      }
    
    matrix->scale(factor);
    }
}




void
FESystemElem::FluidElem::calculate_K_c_Jac(DenseMatrix<double>* matrix,
                                             const unsigned int design_point_enum_ID)
{
  if (!this->analysis_discipline.checkPropertyDependenceOnTemperature())
    return;
  
  static double factor; 
  factor = 0.0; 
  
  this->elem_property_card->getFactorSensitivityForLocalParameter
    (factor, FLUID_CONDUCTANCE_FACTOR::num(), Property::TEMPERATURE::num());

  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  // it is assumed that the matrices and vectors remain the same size as 
  // created here. Hence, they are not resized again.
  static DenseMatrix<double> mat1(n_nodes, n_nodes), mat2(n_nodes, n_nodes);
  static DenseVector<double> vec2(n_nodes);
  if ((vec2.size() != n_nodes))
    {
    vec2.resize(n_nodes);
    mat1.resize(n_nodes, n_nodes); mat2.resize(n_nodes, n_nodes);
    }

  vec2.zero();  
  
  switch (this->dimension)  // get the factors and add them to the matrix dN/dx * dN/dx
    {
    case 3:
      {
        mat1.zero(); mat2.zero();
        this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &mat1, 
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        mat2.add(1.0, mat1);
      }
      
    case 2:
      {
        // dN/dy * dN/dy
        mat1.zero(); mat2.zero();
        mat1.zero();
        this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &mat1, 
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        mat2.add(1.0, mat1);
      }
      
    case 1:
      {
        mat1.zero(); mat2.zero();
        this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &mat1,
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        mat2.add(1.0, mat1);
      }
      break;
      
    default:
      abort();
    }
  
  mat2.scale(factor / (1.0 * n_nodes));
  
  mat2.right_multiply_vector(this->dof_values_vec[0], vec2);
  
  // now, this vector has to be placed inside the jacobian matrix
  static double value = 0.0;
  for (unsigned int i=0; i<n_nodes; i++)
    {
    value = vec2(i);
    for (unsigned int j=0; j < n_nodes; j++)
      (*matrix)(i,j) = value;
    }
}




void
FESystemElem::FluidElem::calculate_K_h(DenseMatrix<double>* matrix, 
                                         const unsigned int design_point_enum_ID,
                                         bool sensitivity_calculation)
{
  // initialize the data structures for loads and their sensitivities
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
    load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(SURFACE_CONVECTION_HEAT_LOAD::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(SURFACE_CONVECTION_HEAT_LOAD::num());
	
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();
	  
  unsigned int side_num;
  bool load_on_elem_face;
	 unsigned int n_sides;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
  
  double load_val, load_val_sens;
  load_val = 0.0; load_val_sens = 0.0;

  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
    side_num = load_iterator->second->getSurfaceID();
		
    // for 2-D elements, the element face itself can carry a surface flux load.
    load_on_elem_face = false;
		
    if (n_sides == side_num)
      load_on_elem_face = true;
		
    double factor, factor_sens;
    factor = 0.0; factor_sens = 0.0;
    
    switch (this->dimension)
      {
      case 1:
        {
          this->elem_property_card->getPropertyValue(AREA::num(), factor);
          if (sensitivity_calculation)
            this->elem_property_card->getPropertyDerivativeForGlobalParameter
              (AREA::num(), this->sensitivity_parameter_ID, factor_sens);
        }
        break;
				
      case 2:
        {
          // the integration will be different if the load is on the element face
          if (load_on_elem_face)
            {
            factor =  1.0;
            factor_sens = 0.0;
            }
          else 
            {
            this->elem_property_card->getPropertyValue(THICKNESS_2D_ELEM::num(), factor);
            if (sensitivity_calculation)
              this->elem_property_card->getPropertyDerivativeForGlobalParameter
                (THICKNESS_2D_ELEM::num(), this->sensitivity_parameter_ID, factor_sens);
            }
        }
        break;
        
      case 3:
        {
          factor = 1.0;
          factor_sens = 0.0;
        }
        break;
      }
		
    const unsigned int domain = this->getDomainEnumFromSideNumber(side_num);
    
    load_val = 
      dynamic_cast<const Loads::SurfaceConvectionLoadCombination*>(load_iterator->second)->getConvectionCoeff();
    
    if (sensitivity_calculation)
      {
      load_sens_map_it = this->surface_load_sens->find(SURFACE_CONVECTION_HEAT_LOAD::num());
      // if the load exists, its sensitivity should also be specified
      Assert(load_sens_map_it != this->surface_load_sens->end(),
             ExcInternalError());
      load_sens_iterator = load_sens_map_it->second.find(side_num);
      AssertThrow(load_sens_iterator != load_sens_map_it->second.end(),
             ExcInternalError());
      
      load_val_sens = dynamic_cast<const SurfaceConvectionLoad*>
        (load_sens_iterator->second)->getConvectionCoeff();

      
      switch (this->sensitivity_parameter)
        {
        case PROPERTY_PARAMETER_ENUM_ID:
          {
            this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), matrix, 
                                     design_point_enum_ID,
                                     domain, LAGRANGE);
            matrix->scale(factor * load_val_sens + factor_sens * load_val);
          }
          break;
          
        case SHAPE_PARAMETER_ENUM_ID:
          {
            static DenseMatrix<double> qty, qty_sens;
            qty.zero(); qty_sens.zero();
            
            this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &qty,
                                     design_point_enum_ID,
                                     domain, LAGRANGE);
            this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_N_FACTOR::num(), 
                                                     &qty_sens, domain, LAGRANGE);
            
            qty.scale(factor * load_val_sens);
            qty_sens.scale(factor * load_val);
            
            matrix->add(1.0, qty);
            matrix->add(1.0, qty_sens);
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
                               domain, LAGRANGE);
      matrix->scale(factor * load_val);
      }
    }
}





void
FESystemElem::FluidElem::calculate_F_qsurf(DenseVector<double>* vector, 
                                             const unsigned int design_point_enum_ID,
                                             bool sensitivity_calculation)
{
  // initialize the data structures for loads and their sensitivities
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
  load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(SURFACE_HEAT_LOAD::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(SURFACE_HEAT_LOAD::num());
		
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();

  unsigned int side_num, n_sides;
  bool load_on_elem_face;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
	
  double load_val, load_val_sens;
  load_val = 0.0; load_val_sens = 0.0;

  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
    side_num = load_iterator->second->getSurfaceID();
		
    // for 2-D elements, the element face itself can carry a surface flux load.
    load_on_elem_face = false;
		
    if (n_sides == side_num)
      load_on_elem_face = true;
		
    double factor, factor_sens;
    factor = 0.0; factor_sens = 0.0;
    
    switch (this->dimension)
      {
      case 1:
        {
          this->elem_property_card->getPropertyValue(AREA::num(), factor);
          if (sensitivity_calculation)
            this->elem_property_card->getPropertyDerivativeForGlobalParameter
            (AREA::num(), this->sensitivity_parameter_ID, factor_sens);
        }
        break;
				
      case 2:
        {
          // the integration will be different if the load is on the element face
          if (load_on_elem_face)
            {
            factor =  1.0;
            factor_sens = 0.0;
            }
          else 
            {
            this->elem_property_card->getPropertyValue(THICKNESS_2D_ELEM::num(), factor);
            if (sensitivity_calculation)
            this->elem_property_card->getPropertyDerivativeForGlobalParameter
              (THICKNESS_2D_ELEM::num(), this->sensitivity_parameter_ID, factor_sens);
            }
        }
        break;
        
      case 3:
        {
          factor = 1.0;
          factor_sens = 0.0;
        }
        break;
      }
    
    const unsigned int domain = this->getDomainEnumFromSideNumber(side_num);
		
    load_val = 
      dynamic_cast<const Loads::ScalarSurfaceLoadCombination*>(load_iterator->second)->getValue();
    
    if (sensitivity_calculation)
      {
      load_sens_map_it = this->surface_load_sens->find(SURFACE_HEAT_LOAD::num());
      // if the load exists, its sensitivity should also be specified
      Assert(load_sens_map_it != this->surface_load_sens->end(),
             ExcInternalError());
      load_sens_iterator = load_sens_map_it->second.find(side_num);
      AssertThrow(load_sens_iterator != load_sens_map_it->second.end(),
                  ExcInternalError());
      
      load_val_sens = dynamic_cast<const Loads::ScalarSurfaceLoadCombination*>
        (load_sens_iterator->second)->getValue();
      
      
      switch (this->sensitivity_parameter)
        {
        case PROPERTY_PARAMETER_ENUM_ID:
          {
            this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), vector, 
                                     design_point_enum_ID,
                                     domain, LAGRANGE);
            vector->scale(factor * load_val_sens + factor_sens * load_val);
          }
          break;
          
        case SHAPE_PARAMETER_ENUM_ID:
          {
            static DenseVector<double> qty, qty_sens;
            qty.zero(); qty_sens.zero();
            
            this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &qty, 
                                     design_point_enum_ID,
                                     domain, LAGRANGE);
            this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_FACTOR::num(),
                                                     &qty_sens, domain, LAGRANGE);
            
            qty.scale(factor * load_val_sens);
            qty_sens.scale(factor * load_val);
            
            vector->add(1.0, qty);
            vector->add(1.0, qty_sens);
          }
          break;
          
        default:
          abort();
          break;
        }
      }
    else
      {
      this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), vector, 
                               design_point_enum_ID,
                               domain, LAGRANGE);
      vector->scale(factor * load_val);
      }
    }
}





void
FESystemElem::FluidElem::calculate_F_Qvol(DenseVector<double>* vector, 
                                            const unsigned int design_point_enum_ID,
                                            bool sensitivity_calculation)
{
  // initialize the data structures for loads and their sensitivities
  std::map<unsigned int, const Loads::VolumeLoadCombination*>::const_iterator load_map_it, load_sens_map_it;
    
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->volume_loads->count(VOLUME_HEAT_LOAD::num()) == 0)
    return;
  else
    load_map_it = this->volume_loads->find(VOLUME_HEAT_LOAD::num());

  double factor, factor_sens;
  factor = 0.0; factor_sens = 0.0;
    
  switch (this->dimension)
    {
    case 1:
      {
        this->elem_property_card->getPropertyValue(AREA::num(), factor);
        if (sensitivity_calculation)
          this->elem_property_card->getPropertyDerivativeForGlobalParameter
            (AREA::num(), this->sensitivity_parameter_ID, factor_sens);
      }
      break;
      
    case 2:
      {
        this->elem_property_card->getPropertyValue(THICKNESS_2D_ELEM::num(), factor);
        if (sensitivity_calculation)
          this->elem_property_card->getPropertyDerivativeForGlobalParameter
            (THICKNESS_2D_ELEM::num(), this->sensitivity_parameter_ID, factor_sens);
      }
      break;
      
    case 3:
      {
        factor = 1.0;
        factor_sens = 0.0;
      }
      break;
    }
  
  double load_val, load_val_sens;
  load_val = 0.0; load_val_sens = 0.0;

  load_val = dynamic_cast<const Loads::ScalarVolumeLoadCombination*>(load_map_it->second)->getValue();
  
  // now, scale the matrix with the right factors to take into 
  // account the sensitivity if needed
  if (sensitivity_calculation)
    {
    load_sens_map_it = this->volume_load_sens->find(VOLUME_HEAT_LOAD::num()); 
    // if the load exists, its sensitivity should also be specified 
    Assert(load_sens_map_it != this->volume_load_sens->end(),
           ExcInternalError());
    
    load_val_sens = dynamic_cast<const Loads::ScalarVolumeLoadCombination*>(load_sens_map_it->second)->getValue();
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), vector, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          vector->scale(factor * load_val_sens + factor_sens * load_val);
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          static DenseVector<double> qty, qty_sens;
          qty.zero(); qty_sens.zero();
          
          this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &qty, 
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_FACTOR::num(), &qty_sens,
                                                   FESystemElem::ELEM_VOLUME::num(),
                                                   LAGRANGE);
          
          qty.scale(factor * load_val_sens);
          qty_sens.scale(factor * load_val);
          
          vector->add(1.0, qty);
          vector->add(1.0, qty_sens);
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else
    {
    this->getFESystemElemQty( FESystemElem::N_FACTOR::num(), vector, 
                              design_point_enum_ID,
                              FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    vector->scale(factor * load_val);
    }
}





void
FESystemElem::FluidElem::calculate_F_h(DenseVector<double>* vector, 
                                         const unsigned int design_point_enum_ID,
                                         bool sensitivity_calculation)
{
  // initialize the data structures for loads and their sensitivities
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
    load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(SURFACE_CONVECTION_HEAT_LOAD::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(SURFACE_CONVECTION_HEAT_LOAD::num());
  
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();
	
  unsigned int side_num;
  bool load_on_elem_face;
  unsigned int n_sides;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
   
  double coeff_val, coeff_val_sens, temp_val, temp_val_sens;
  coeff_val = 0.0; coeff_val_sens = 0.0; temp_val = 0.0; temp_val_sens = 0.0;
   
  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
    side_num = load_iterator->second->getSurfaceID();
		
    // for 2-D elements, the element face itself can carry a surface flux load.
    load_on_elem_face = false;
		
    if (n_sides == side_num)
      load_on_elem_face = true;
    
    double factor, factor_sens;
    factor = 0.0; factor_sens = 0.0;
    
    switch (this->dimension)
      {
      case 1:
        {
          this->elem_property_card->getPropertyValue(AREA::num(), factor);
          if (sensitivity_calculation)
            this->elem_property_card->getPropertyDerivativeForGlobalParameter
            (AREA::num(), this->sensitivity_parameter_ID, factor_sens);
        }
        break;
				
      case 2:
        {
          // the integration will be different if the load is on the element face
          if (load_on_elem_face)
            {
            factor =  1.0;
            factor_sens = 0.0;
            }
          else 
            {
            this->elem_property_card->getPropertyValue(THICKNESS_2D_ELEM::num(), factor);
            if (sensitivity_calculation)
            this->elem_property_card->getPropertyDerivativeForGlobalParameter
              (THICKNESS_2D_ELEM::num(), this->sensitivity_parameter_ID, factor_sens);
            }
        }
        break;
        
      case 3:
        {
          factor = 1.0;
          factor_sens = 0.0;
        }
        break;
      }
        
    const unsigned int domain = this->getDomainEnumFromSideNumber(side_num);
    
    coeff_val = 
      dynamic_cast<const Loads::SurfaceConvectionLoadCombination*>(load_iterator->second)->getConvectionCoeff();
    temp_val = 
      dynamic_cast<const Loads::SurfaceConvectionLoadCombination*>(load_iterator->second)->getAmbientTemperature();
    
    // now, scale the matrix with the right factors to take into 
    // account the sensitivity if needed
    if (sensitivity_calculation)
      {
      load_sens_map_it = this->surface_load_sens->find(SURFACE_CONVECTION_HEAT_LOAD::num());
      // if the load exists, its sensitivity should also be specified
      Assert(load_sens_map_it != this->surface_load_sens->end(),
             ExcInternalError());
      load_sens_iterator = load_sens_map_it->second.find(side_num);
      AssertThrow(load_sens_iterator != load_sens_map_it->second.end(),
                  ExcInternalError());
      
      coeff_val_sens = dynamic_cast<const Loads::SurfaceConvectionLoadCombination*>
        (load_sens_iterator->second)->getConvectionCoeff();
      temp_val_sens = dynamic_cast<const Loads::SurfaceConvectionLoadCombination*>
        (load_sens_iterator->second)->getAmbientTemperature();

      
      switch (this->sensitivity_parameter)
        {
        case PROPERTY_PARAMETER_ENUM_ID:
          {
            this->getFESystemElemQty( FESystemElem::N_FACTOR::num(), vector,
                                      design_point_enum_ID,
                                      domain, LAGRANGE);
            vector->scale(factor * (coeff_val_sens *temp_val + 
                                    coeff_val * temp_val_sens) + 
                          factor_sens * coeff_val * temp_val);
          }
          break;
          
        case SHAPE_PARAMETER_ENUM_ID:
          {
            static DenseVector<double> qty, qty_sens;
            qty.zero(); qty_sens.zero();
            
            this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &qty,
                                     design_point_enum_ID,
                                     domain, LAGRANGE);
            this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_FACTOR::num(),
                                                     &qty_sens, domain, LAGRANGE);
            
            qty.scale(factor * (coeff_val_sens *temp_val + 
                                coeff_val * temp_val_sens));
            qty_sens.scale(factor * coeff_val * temp_val);
            
            vector->add(1.0, qty);
            vector->add(1.0, qty_sens);
          }
          break;
          
        default:
          abort();
          break;
        }
      }
    else
      {
      this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), vector,
                               design_point_enum_ID,
                               domain, LAGRANGE);
      vector->scale(factor * coeff_val * temp_val);
      }
    }
}





void
FESystemElem::FluidElem::calculate_F_sigma(DenseVector<double>* vector, 
                                             const unsigned int design_point_enum_ID,
                                             bool sensitivity_calculation)
{
  // initialize the data structures for loads and their sensitivities
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
  load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(SURFACE_RADIATION_HEAT_LOAD::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(SURFACE_RADIATION_HEAT_LOAD::num());
		
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();
  
  unsigned int side_num, n_sides;
  bool load_on_elem_face;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
	
  double load_val, load_val_sens;
  load_val = 0.0; load_val_sens = 0.0;
	
  static const double absolute_temp =
    this->analysis_discipline.getFESystemController().analysis_case->getRealParameter("ABSOLUTE_TEMP");
  static const double sb_constant =
    this->analysis_discipline.getFESystemController().analysis_case->getRealParameter("S_B_CONSTANT");
		
  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
    side_num = load_iterator->second->getSurfaceID();
		
    // for 2-D elements, the element face itself can carry a surface flux load.
    load_on_elem_face = false;
		
    if (n_sides == side_num)
      load_on_elem_face = true;
		
    double factor, factor_sens;
    factor = 0.0; factor_sens = 0.0;
    
    switch (this->dimension)
      {
      case 1:
        {
          this->elem_property_card->getFactor(factor, FLUID_EMITTED_LOAD_FACTOR::num());
          this->elem_property_card->getFactorSensitivityForGlobalParameter
            (factor_sens, FLUID_EMITTED_LOAD_FACTOR::num(), this->sensitivity_parameter_ID);
        }
        break;
				
      case 2:
        {
          // the integration will be different if the load is on the element face
          if (load_on_elem_face)
            {
            this->elem_property_card->getPropertyValueFromMaterialCard(EMISSIVITY::num(),
                                                                       factor);
            
            this->elem_property_card->getPropertyValueDerivativeForGlobalParameterFromMaterialCard
              (EMISSIVITY::num(), this->sensitivity_parameter_ID, factor_sens);
            }
          else 
            {
            this->elem_property_card->getFactor(factor, FLUID_EMITTED_LOAD_FACTOR::num());
            this->elem_property_card->getFactorSensitivityForGlobalParameter
              (factor_sens, FLUID_EMITTED_LOAD_FACTOR::num(), this->sensitivity_parameter_ID);
            }
        }
        break;
        
      case 3:
        {
          this->elem_property_card->getPropertyValueFromMaterialCard(EMISSIVITY::num(),
                                                                     factor);
          this->elem_property_card->getPropertyValueDerivativeForGlobalParameterFromMaterialCard
            (EMISSIVITY::num(), this->sensitivity_parameter_ID, factor_sens);
        }
        break;
      }

    factor *= sb_constant;
    factor_sens *= sb_constant;		
    
 
    const unsigned int domain = this->getDomainEnumFromSideNumber(side_num);
    
    static DenseVector<double> qty_face, qty_sigma, qty_sigma_sens, qty_face_sens;
    qty_face.zero(); qty_sigma.zero(); qty_sigma_sens.zero(); qty_face_sens.zero();

    load_val = 
      dynamic_cast<const Loads::SurfaceRadiationLoadCombination*>(load_iterator->second)->getTemperature();

    // now, scale the matrix with the right factors to take into 
    // account the sensitivity if needed
    if (sensitivity_calculation)
      {
      load_sens_map_it = this->surface_load_sens->find(FLUID_EMITTED_LOAD_FACTOR::num());
      // if the load exists, its sensitivity should also be specified
      Assert(load_sens_map_it != this->surface_load_sens->end(),
             ExcInternalError());
      load_sens_iterator = load_sens_map_it->second.find(side_num);
      AssertThrow(load_sens_iterator != load_sens_map_it->second.end(),
                  ExcInternalError());
      
      load_val_sens = dynamic_cast<const Loads::SurfaceRadiationLoadCombination*>
        (load_sens_iterator->second)->getTemperature();
      
      
      switch (this->sensitivity_parameter)
        {
        case PROPERTY_PARAMETER_ENUM_ID:
          {
            this->getFluidElemQty(FESystemElem::FLUID_F_SIGMA_FACTOR::num(), 
                                    &qty_sigma,  design_point_enum_ID,
                                    domain, false);
            
            this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &qty_face, 
                                     design_point_enum_ID,
                                     domain, LAGRANGE);
            
            
            vector->add(-factor_sens, qty_sigma);
            
            vector->add((factor_sens * pow((load_val + absolute_temp), 4)+
                         4.0 * factor * pow((load_val + absolute_temp), 3) * load_val_sens), 
                        qty_face);
          }
          break;
          
        case SHAPE_PARAMETER_ENUM_ID:
          {
            this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &qty_face, 
                                     design_point_enum_ID, domain, LAGRANGE);
            
            this->getFluidElemQty
              (FESystemElem::FLUID_F_SIGMA_FACTOR::num(), &qty_sigma_sens, 
               design_point_enum_ID, domain, true);
            
            this->getFESystemElemQtyShapeSensitivity
              (FESystemElem::N_FACTOR::num(), &qty_face_sens, domain, LAGRANGE);
            
            
            vector->add(-factor, qty_sigma_sens);
            vector->add(4.0 * factor * pow((load_val + absolute_temp), 3) * load_val_sens, 
                        qty_face);
            vector->add(factor * pow((load_val + absolute_temp), 4),
                        qty_face_sens);
          }
          break;
          
        default:
          abort();
          break;
        }
      }
    else
      {
      this->getFluidElemQty(FESystemElem::FLUID_F_SIGMA_FACTOR::num(), &qty_sigma,
                              design_point_enum_ID, domain, false);
      
      this->getFESystemElemQty(FESystemElem::N_FACTOR::num(), &qty_face, 
                               design_point_enum_ID, domain, LAGRANGE);
      
      vector->add(-factor, qty_sigma);
      vector->add(factor * pow((load_val + absolute_temp), 4),
                  qty_face);
      }
    }
  
}





void
FESystemElem::FluidElem::calculate_F_sigma_Jac(DenseMatrix<double>* matrix,
                                                 const unsigned int design_point_enum_ID)
{
  // initialize the data structures for loads and their sensitivities
  std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> >::const_iterator 
  load_map_it, load_sens_map_it;
  std::map<unsigned int, const Loads::SurfaceLoadCombination*>::const_iterator 
  load_iterator, end_iterator,		// iterators for the loads
  load_sens_iterator;		// iterators for the load sensitivities
  
  // if no loads exist return from the function with an arguement false
  // this implies that the quantity is zero and need not be processed
  if (this->surface_loads->count(SURFACE_RADIATION_HEAT_LOAD::num()) == 0)
    return;
  else
    load_map_it = this->surface_loads->find(SURFACE_RADIATION_HEAT_LOAD::num());
		
  // obtain the two iterators as separate values
  load_iterator = load_map_it->second.begin();
  end_iterator = load_map_it->second.end();
  
  unsigned int side_num, n_sides;
  bool load_on_elem_face;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();
	
  double load_val, load_val_sens;
  load_val = 0.0; load_val_sens = 0.0;
	
  static const double absolute_temp =
    this->analysis_discipline.getFESystemController().analysis_case->getRealParameter("ABSOLUTE_TEMP");
  static const double sb_constant =
    this->analysis_discipline.getFESystemController().analysis_case->getRealParameter("S_B_CONSTANT");
  
  // process the load and update the element load vector.
  for (; load_iterator != end_iterator; load_iterator++)
    {
    side_num = load_iterator->second->getSurfaceID();
		
    // for 2-D elements, the element face itself can carry a surface flux load.
    load_on_elem_face = false;
		
    FEBase *fe_base_local = NULL;
    QBase *qbase_local = NULL;
    const std::vector<std::vector<Real> > *phi = NULL;
    const std::vector<Real>* JxW = NULL;

    if (n_sides == side_num)
      {
      load_on_elem_face = true;
      fe_base_local = this->fe_base_map_for_DV[design_point_enum_ID][LAGRANGE];
      qbase_local = this->qbase_map_for_DV[design_point_enum_ID][LAGRANGE];
      }
    else if (this->dimension > 1)
      {
      this->initialize_element_side(side_num, design_point_enum_ID);
      fe_base_local = this->fe_side_map_for_DV[design_point_enum_ID][LAGRANGE];
      qbase_local = this->qbase_side_map_for_DV[design_point_enum_ID][LAGRANGE];
      }
    
    if (this->dimension > 1)
      {
      JxW = &(fe_base_local->get_JxW());
      phi = &(fe_base_local->get_phi());
      }
		
    double factor;
    factor = 0.0; 
    
    switch (this->dimension)
      {
      case 1:
        {
          this->elem_property_card->getFactor(factor, FLUID_EMITTED_LOAD_FACTOR::num());
        }
        break;
				
      case 2:
        {
          // the integration will be different if the load is on the element face
          if (load_on_elem_face)
            this->elem_property_card->getPropertyValueFromMaterialCard
              (EMISSIVITY::num(), factor);
          else 
            this->elem_property_card->getFactor(factor, FLUID_EMITTED_LOAD_FACTOR::num());
        }
        break;
        
      case 3:
        {
          this->elem_property_card->getPropertyValueFromMaterialCard(EMISSIVITY::num(),
                                                                     factor);
        }
        break;
      }
    
    factor *= (4.0 * sb_constant);
    
    switch (this->dimension)
      {
      case 1:
        {
          switch (side_num)
            {
            case 0:
              (*matrix)(0,0) += -factor * pow((this->dof_values_vec[0](0)+absolute_temp), 3);
              break;
            case 1:
              (*matrix)(1,1) += -factor * pow((this->dof_values_vec[0](1)+absolute_temp), 3);
              break;
            default:
              error();
              break;
            }
        }
        break;
				
      case 2:
      case 3:
        {
          const unsigned int n_q_points = qbase_local->n_points();
          
          for (unsigned int qp=0; qp<n_q_points; qp++)
            {
            // calculate the temperature at this quadrature point
            double temp = 0.0;
            for (unsigned int i=0; i<(*phi).size(); i++)
              temp += (*phi)[i][qp]* this->dof_values_vec[0](i);
            
            for (unsigned int i=0; i<(*phi).size(); i++)
              for (unsigned int j=0; j<(*phi).size(); j++)
                (*matrix)(i,j) += - (*JxW)[qp]* factor * pow((temp+absolute_temp),3)*
                  ((*phi)[i][qp]* (*phi)[j][qp]);
            }
        }
        break;
      }
    }	
  
  
  // now, to add the effect of temperature dependent property, the following needs to be done. 
  // the whole operation is simplified by directly using the radiative load
  
  // create a local static vector to be used for this computation
  unsigned int n_nodes;
  n_nodes = this->getNNodes();
  static DenseVector<double> vec1(n_nodes);
  double eps, eps_temp_sens;
  
  if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
    {
    vec1.zero();
    eps = 0.0; eps_temp_sens = 0.0;
    this->getElementAssembledQty(FESystemElem::FLUID_F_EMITTED_RAD_VECTOR::num(),
                                 &vec1);
    this->elem_property_card->getPropertyValueFromMaterialCard(EMISSIVITY::num(), eps);
    this->elem_property_card->getPropertyValueDerivativeForLocalParameterFromMaterialCard
      (EMISSIVITY::num(), Property::TEMPERATURE::num(), eps_temp_sens);
    
    // this scaling is based on the assumption that the property is calculated on the 
    // average temperature per elem.
    vec1.scale(eps_temp_sens / (1.0 * n_nodes) / eps);
    
    // now add this to the vector
    double val = 0.0;
    for (unsigned int i=0; i<n_nodes; i++)
      {
      val = vec1(i);
      for (unsigned int j=0; j<n_nodes; j++)
        (*matrix)(i,j) += val;
      }
    }
}




void 
FESystemElem::FluidElem::calculateEmittedRadiationLoadFactor
(DenseVector<double>* vector,
 const unsigned int design_point_enum_ID,
 const unsigned int domain)
{
  static const double absolute_temp =
  this->analysis_discipline.getFESystemController().analysis_case->getRealParameter("ABSOLUTE_TEMP");

  // for 2-D elements, the element face itself can carry a surface flux load.
  bool load_on_elem_face;
  unsigned int side_num;
  unsigned int n_sides;
  n_sides = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second->n_sides();

  side_num = this->getSideNumberFromDomainEnum(domain);
  load_on_elem_face = false;
  
  FEBase *fe_base_local = NULL;
  QBase *qbase_local = NULL;
  const std::vector<std::vector<Real> > *phi = NULL;
  const std::vector<Real>* JxW = NULL;
  
  if (n_sides == side_num)
    {
    load_on_elem_face = true;
    fe_base_local = this->fe_base_map_for_DV[design_point_enum_ID][LAGRANGE];
    qbase_local = this->qbase_map_for_DV[design_point_enum_ID][LAGRANGE];
    }
  else if (this->dimension > 1)
    {
    this->initialize_element_side(side_num, design_point_enum_ID);
    fe_base_local = this->fe_side_map_for_DV[design_point_enum_ID][LAGRANGE];
    qbase_local = this->qbase_side_map_for_DV[design_point_enum_ID][LAGRANGE];
    }
  
  if (this->dimension > 1)
    {
    JxW = &(fe_base_local->get_JxW());
    phi = &(fe_base_local->get_phi());
    }

  DenseVector<double>& vec_ref = *vector;
  
  switch (this->dimension)
    {
    case 1:
      {
        switch (side_num)
          {
          case 0:
            vec_ref(0) += pow((this->dof_values_vec[0](0)+absolute_temp), 4);
            break;
          case 1:
            vec_ref(1) += pow((this->dof_values_vec[0](1)+absolute_temp), 4);
            break;
          default:
            error();
            break;
          }
      }
      break;
      
    case 2:
    case 3:
      {
        const unsigned int n_q_points = qbase_local->n_points();
        double temp_power=0.0, temp = 0.0;
        for (unsigned int qp=0; qp<n_q_points; qp++)
          {
          // calculate the temperature at this quadrature point
          temp = 0.0;
          for (unsigned int i=0; i<(*phi).size(); i++)
            temp += (*phi)[i][qp]* this->dof_values_vec[0](i);
          temp_power = pow(temp+absolute_temp, 4);
          
          for (unsigned int i=0; i<(*phi).size(); i++)
            vec_ref(i) += (*JxW)[qp]* (*phi)[i][qp] * temp_power;
          }
      }        
      break;
    }
}




std::auto_ptr<ElemPostProcessQty>
FESystemElem::FluidElem::getElementPostProcessQty
(std::vector<unsigned int> load_cases,
 std::vector<DesignData::DesignParameter*> DV_vector)
{
  (void) load_cases;
  (void) DV_vector;

  std::auto_ptr<ElemPostProcessQty> return_qty;
  
  
  return return_qty;
}




void 
FESystemElem::FluidElem::calculateFluidElemQty(DenseMatrix<double>* quantity,
                                                   const unsigned int qty_name,
                                                   const unsigned int design_point_enum_ID,
                                                   const unsigned int domain)
{
  // unused parameter here
  (void) domain;

  assert (quantity != NULL);
  
  // check if the local element has been initialized or not. If not, 
  // initialize it, also, resize the quantity
  if (!this->local_elem_is_initialized_for_DV_map.find(design_point_enum_ID)->second)
    this->initialize_element(design_point_enum_ID);
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  if (!(quantity->m() == n_nodes && quantity->n() == n_nodes))
    {
    quantity->resize(n_nodes, n_nodes);
    quantity->zero();
    }
  
  switch(qty_name)
    {
    case FLUID_F_SIGMA_FACTOR_ENUM_ID:
    default:
      abort();
      break;
    }
}




void 
FESystemElem::FluidElem::calculateFluidElemQty(DenseVector<double>* quantity,
                                                   const unsigned int qty_name,
                                                   const unsigned int design_point_enum_ID,
                                                   const unsigned int domain)
{
  assert (quantity != NULL);
	
  // check if the local element has been initialized or not. If not, 
  // initialize it, also, resize the quantity
  if (!this->local_elem_is_initialized_for_DV_map.find(design_point_enum_ID)->second)
    this->initialize_element(design_point_enum_ID);
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  if ((quantity->size() != n_nodes))
    {
    quantity->resize(n_nodes);
    quantity->zero();
    }
  
  switch(qty_name)
    {
    case FLUID_F_SIGMA_FACTOR_ENUM_ID:
      {
        this->calculateEmittedRadiationLoadFactor(quantity, design_point_enum_ID, domain);
      }
      break;
			
    default:
      abort();
      break;
    }
}



void
FESystemElem::FluidElem::initPropertyCard()
{
  Assert(!this->property_card_initialized, ExcInternalError());

  // if the properties are temperature dependent, then initialize the
  // element property cards at the element average temperature
  if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
    {
    // get the load case number and the dof values
    static unsigned int n_nodes;
    n_nodes = this->getNNodes();
    static std::map<unsigned int, double> local_param_map;
    
    // now calculate the average
    static double temp;
    temp = 0.0;
    
    for (unsigned int i=0; i < n_nodes; i++)
      temp += this->dof_values_vec[0](i);
    
    temp /= (1.0 * n_nodes);
    
    local_param_map[Property::TEMPERATURE::num()] = temp; 
    this->elem_property_card->reinitElemAndMaterialCardForLocalParameters(&local_param_map);
    }
  else
    this->elem_property_card->reinitElemAndMaterialCardForLocalParameters();
  
  this->property_card_initialized = true;
}


void
FESystemElem::FluidElem::clearPropertyCardInitialization()
{
  static std::vector<unsigned int> param_vec;
  static std::vector<unsigned int> no_param_vec;
  
  if (param_vec.size() != 0)
   param_vec.push_back(Property::TEMPERATURE::num());
  
  // if the properties are temperature dependent, then initialize the
  // element property cards at the element average temperature
  if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
    this->elem_property_card->partialClearElemAndMaterialCardLocalParameterInitialization
      (param_vec);
  else
    this->elem_property_card->partialClearElemAndMaterialCardLocalParameterInitialization
      (no_param_vec);
  
  this->property_card_initialized = false;
}

