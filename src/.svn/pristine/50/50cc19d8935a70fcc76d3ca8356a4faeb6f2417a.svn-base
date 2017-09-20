// $Id: bar.C,v 1.25.6.1 2007-03-14 22:05:03 manav Exp $

// C++ includes


// FESystem includes
#include "StructuralElems/bar.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic1DElemDataCard.h"
#include "FESystem/AnalysisCase.h"

// libMesh includes
#include "geom/edge_edge2.h"

FESystemElem::Bar::Bar(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::StructuralElem(1, FESystemElem::STRUCTURAL_BAR_EDGE2::num(), discipline)
{
  
}



FESystemElem::Bar::~Bar()
{
	
}


void
FESystemElem::Bar::calculate_K_G(DenseMatrix<double>* matrix, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation)
{
  // params not used here
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;
    
  matrix->zero();
}



void
FESystemElem::Bar::getFETypes(std::vector<FEType>& fetypes)
{
  // this element calculates the matrices as closed form expressions.
  // Hence, no finite element is needed
  fetypes.clear();
}





void
FESystemElem::Bar::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  // this element does not use and FE, hence, no FE are needed.
  quadratures.clear();
}




void
FESystemElem::Bar::calculate_M(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation)
{
  static Elem* elem;
  elem = this->geometric_elems_for_DV_map[design_point_enum_ID];
  
  // calculate the length of the element
  const Point& point0 = elem->point(0);
  const Point& point1 = elem->point(1);
  
  static double length, dx, dy, dz;
  dx = point0(0)-point1(0);
  dy = point0(1)-point1(1);
  dz = point0(2)-point1(2);
  length = pow((dx*dx+dy*dy+dz*dz),0.5);
	
	
  static double factor = 0.0, factor_sens = 0.0, factor_temp_sens = 0.0; 
  factor = 0.0; factor_sens = 0.0, factor_temp_sens = 0.0;
  this->elem_property_card->getFactor(factor, MASS_FACTOR::num());
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (factor_sens, MASS_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (factor_temp_sens, MASS_FACTOR::num(), Property::TEMPERATURE::num());
        break;
      }
    }
  
  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(2);
  nodal_temp_sens.zero();
  
  // set the elements for the matrix.
  (*matrix)(0,0) = (2.0/6.0); // u1,u1
  (*matrix)(0,1) = (1.0/6.0); // u1,u2
  (*matrix)(1,0) = (1.0/6.0); // u2,u1
  (*matrix)(1,1) = (2.0/6.0); // u2,u2
	
	
  if (sensitivity_calculation)
    {
    static double avg_factor_temp_sens = 0.0;
    avg_factor_temp_sens = 0.0;

    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
      avg_factor_temp_sens = factor_temp_sens * 0.5 * (nodal_temp_sens(0) + nodal_temp_sens(1));
      }
    
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          matrix->scale((factor_sens + avg_factor_temp_sens) * length);
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          static Elem* perturbed_elem;
          perturbed_elem = 
            this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()];
          // calculate the length for the perturbed element, and find the sensitivity
          // calculate the length of the element
          const Point& point0 = perturbed_elem->point(0);
          const Point& point1 = perturbed_elem->point(1);
          
          static double perturbed_length;
          dx = point0(0)-point1(0);
          dy = point0(1)-point1(1);
          dz = point0(2)-point1(2);
          perturbed_length = pow((dx*dx+dy*dy+dz*dz),0.5);
          
          matrix->scale(avg_factor_temp_sens * length +
                        factor * (perturbed_length - length)/ this->perturbation);
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else
    matrix->scale(factor * length);
  
}




void
FESystemElem::Bar::calculate_K(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation)
{
  static Elem* elem;
  elem = this->geometric_elems_for_DV_map[design_point_enum_ID];

  // calculate the length of the element
  const Point& point0 = elem->point(0);
  const Point& point1 = elem->point(1);
	
  static double length, dx, dy, dz;
  {
    double dx = point0(0)-point1(0);
    double dy = point0(1)-point1(1);
    double dz = point0(2)-point1(2);
    length = pow((dx*dx+dy*dy+dz*dz),0.5);
  }
	
	
   static double factor = 0.0, factor_sens = 0.0, factor_temp_sens = 0.0; 
  factor = 0.0; factor_sens = 0.0, factor_temp_sens = 0.0;

  this->elem_property_card->getFactor(factor, EA_FACTOR::num());
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (factor_sens, EA_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (factor_sens, EA_FACTOR::num(), Property::TEMPERATURE::num());
        break;
      }
    }
  
  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(2);
  nodal_temp_sens.zero();

  // set the matrix elements
  (*matrix)(0,0) = 1.0; // u1,u1
  (*matrix)(0,1) = -1.0; // u1,u2
  (*matrix)(1,0) = -1.0; // u2,u1
  (*matrix)(1,1) = 1.0; // u2,u2
	
	
  if  (sensitivity_calculation)
    {
    static double avg_factor_temp_sens = 0.0;
    avg_factor_temp_sens = 0.0;
    
    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
      avg_factor_temp_sens = factor_temp_sens * 0.5 * (nodal_temp_sens(0) + nodal_temp_sens(1));
      }
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          matrix->scale((factor_sens + avg_factor_temp_sens)/ length);
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          static Elem* perturbed_elem;
          perturbed_elem = 
            this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()];
          // calculate the length for the perturbed element, and find the sensitivity
          // calculate the length of the element
          const Point& point0 = perturbed_elem->point(0);
          const Point& point1 = perturbed_elem->point(1);
          
          static double perturbed_length;
          dx = point0(0)-point1(0);
          dy = point0(1)-point1(1);
          dz = point0(2)-point1(2);
          perturbed_length = pow((dx*dx+dy*dy+dz*dz),0.5);
          
          matrix->scale(avg_factor_temp_sens / length -
                        (factor / (length*length))* 
                        (perturbed_length - length)/ this->perturbation);
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else
    {
    matrix->scale(factor / length);
    }
}





void
FESystemElem::Bar::calculate_F_T(DenseVector<double>* vector, 
                                const unsigned int design_point_enum_ID,
                                bool sensitivity_calculation)
{
  // param not used here
  (void) design_point_enum_ID;

  // the load vector is calculated based on the following expression
  // 
  //					| -1    -1	| |T1|					|1 |
  // F_T =  (E*A*alpha/2.0)*	|		| |   | + (E*A*alpha* T_ref)	|  |
  //					| 1  1		| |T2|					|-1|
  //
  //
  
  static double ref_temp = this->analysis_discipline.getFESystemController().
  analysis_case->getRealParameter("REFERENCE_TEMP");
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  // initialize the data structures for the load
  static DenseVector<double> nodal_temp(n_nodes), nodal_temp_sens(n_nodes);
  
  if (nodal_temp.size() != n_nodes)
    {
    nodal_temp.resize(n_nodes);
    nodal_temp_sens.resize(n_nodes);
    }
  
  nodal_temp.zero(); nodal_temp_sens.zero();
  
  this->extractNodalTemperatureVectorFromLoads(nodal_temp, false);
  
  static double factor = 0.0, factor_sens = 0.0, factor_temp_sens = 0.0;
  factor = 0.0; factor_sens = 0.0; factor_temp_sens = 0.0;
  this->elem_property_card->getFactor(factor, THERMAL_EXPANSION_FACTOR::num());
  
  if (sensitivity_calculation)
    {
    // get the nodal temperature vector sensitivity
    this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (factor_sens, THERMAL_EXPANSION_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (factor_temp_sens, THERMAL_EXPANSION_FACTOR::num(), Property::TEMPERATURE::num());
        break;
      }
    }
  
  // now, scale the matrix with the right factors to take into 
  // account the sensitivity if needed
  if (sensitivity_calculation)
    {
    // this load vector does not have any shape effect, hence, the
    // shape and property sensitivity are both similar, where 
    // the factor_sens will be zero for shape sensitivity
    // also, the reference temperature is assumed to be constant
    (*vector)(0) = 
    - (factor_sens + factor_temp_sens * 0.5 * (nodal_temp_sens(0) + nodal_temp_sens(1))) * 
    0.5 * (nodal_temp(0) + nodal_temp(1) - 2.0 * ref_temp) +
    factor * 0.5 * (nodal_temp_sens(0) + nodal_temp_sens(1));
    (*vector)(1) = -(*vector)(0);
    }
  else
    {
    (*vector)(0) = - factor * 0.5 * (nodal_temp(0) + nodal_temp(1)) + factor * ref_temp; // u1
    (*vector)(1) =  -(*vector)(0); // u2
    }
}



std::auto_ptr<ElemPostProcessQty > 
FESystemElem::Bar::getElementPostProcessQty(std::vector<unsigned int> load_cases, 
                                            std::vector<DesignData::DesignParameter*> dv_vector)
{
  (void) load_cases;
  (void) dv_vector;
//   std::auto_ptr<ElemPostProcessQty> return_qty(new ElemPostProcessQty(this->elem_ID));
  
//   TensorValue<double> tensor;
  
//   unsigned int n_nodes = this->getNNodes();
  
//   // iterate over all the load cases, to calculate the strains and stresses
	
//   std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
//   std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
	
//   std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
//   dv_begin = dv_vector.begin();
//   dv_end = dv_vector.end();
	
	
//   double strain = 0.0, strain_sens = 0.0,
//     stress = 0.0, stress_sens = 0.0,
//     mechanical_strain = 0.0, mechanical_strain_sens = 0.0,
//     material_factor = 0.0, material_factor_sens = 0.0,
//     temp_strain_material_factor = 0.0, temp_strain_material_factor_sens = 0.0,
//     ref_temp = 0.0,
//     length = 0.0,
//     length_factor = 0.0, perturbed_length_factor = 0.0;
  
  
//   std::auto_ptr<DenseMatrix<double> > transform_mat(new DenseMatrix<double>), 
//     transform_mat_sens(new DenseMatrix<double>);
  
//   this->getStructuralT_matrix(transform_mat.get(), FESystemElem::BASE_ELEM::num(), false);
  
  
//   // also, create dof_vectors
//   DenseVector<double> dof_values(6*n_nodes), dof_value_sens(6*n_nodes),
//     local_dof(6*n_nodes), local_dof_sens(6*n_nodes),scratch_vec(6*n_nodes),
//     nodal_temp(n_nodes), nodal_temp_sens(n_nodes);
  
//   Elem *base_elem, *perturbed_elem;
//   base_elem = this->geometric_elems_for_DV_map[FESystemElem::BASE_ELEM::num()];
//   perturbed_elem = this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()];
  
//   // calculate the length of the element and set the factors 
//   const Point& point0 = base_elem->point(0);
//   const Point& point1 = base_elem->point(1);
	
//   {
//     double dx=0.0, dy=0.0, dz=0.0;
//     dx = point0(0)-point1(0);
//     dy = point0(1)-point1(1);
//     dz = point0(2)-point1(2);
//     length = pow((dx*dx+dy*dy+dz*dz),0.5);
//   }
  
//   length_factor = 1.0/ length;
  
//   this->elem_property_card->getPropertyValueFromMaterialCard(YOUNGS_MODULUS::num(),
//                                                              material_factor);
//   this->elem_property_card->getPropertyValueFromMaterialCard(ALPHA_EXPANSION::num(),
//                                                              temp_strain_material_factor);
//   this->elem_property_card->getPropertyValueFromMaterialCard(TEMP_REF::num(),
//                                                              ref_temp);

//   double temp = 0.0, temp_sens = 0.0;
  
//   // iterate over each load case, and ask solver to solve for it
//   for (; load_case_it != load_case_end; load_case_it++)
//     {
    
//     // get the DOF vector for this elem
//      this->analysis_discipline.getElemDofValues(base_elem, 
//                                                 dof_values,
//                                                 *load_case_it);
    
//     // transform these dofs to the local coordinate system
//     local_dof.zero();
//     transform_mat->right_multiply_vector(dof_values, local_dof);
    
//     this->getNodalTemperatureVector(nodal_temp, *load_case_it, false );
    
//     // calculate the factor and the strain and stress
//     strain = (local_dof(1) - local_dof(0)) * length_factor;
//     temp = (nodal_temp(0) + nodal_temp(1) ) * 0.5 - ref_temp ; // temperature at the middle of element
//     mechanical_strain = (strain - temp_strain_material_factor * temp);
//     stress =   material_factor * mechanical_strain;
    
//     // add the two tensors
//     tensor.zero();
//     tensor(0,0) = strain;
//     return_qty->addStrainTensor(tensor, *load_case_it);
    
//     tensor.zero();
//     tensor(0,0) = stress;
//     return_qty->addStressTensor(tensor, *load_case_it);
    
//     tensor.zero();
//     tensor(0,0) = mechanical_strain;
//     return_qty->addMechanicalStrainTensor(tensor, *load_case_it);
    
//     unsigned int dv_ID = 0;
//     // now calculate the sensitivities
//     dv_it = dv_begin;
//     for (; dv_it != dv_end; dv_it++)
//       {
//       dv_ID = (*dv_it)->getID();
//       this->clearSensitivityInitialization();
      
//       // get the dof value sensitivity for this case
//       this->analysis_discipline.getElemDofValues(base_elem,
//                                                  dof_value_sens,
//                                                  *load_case_it,
//                                                  true,
//                                                  dv_ID);
      
//       // get the load sensitivity vectors
//       this->getNodalTemperatureVector(nodal_temp_sens, *load_case_it, true, dv_ID);
      
//       material_factor_sens = 0.0;
//       temp_strain_material_factor_sens = 0.0;
            
//       switch ((*dv_it)->getParameterTypeEnumID())
//         {
//         case PROPERTY_PARAMETER_ENUM_ID:
//           {
//             this->reinitForPropertySensitivity(dv_ID);
            
//             perturbed_length_factor = 0.0;
            
//             // if the property ID is the same as this elements property, then set the 
//             // value of the property sensitivity. Else the value be zero
//             if (this->elem_property_card->checkGlobalParameterDependence(dv_ID))
//               {
//               this->elem_property_card->getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//               (YOUNGS_MODULUS::num(), dv_ID, material_factor_sens);
//               this->elem_property_card->getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//                 (ALPHA_EXPANSION::num(), dv_ID, temp_strain_material_factor_sens);
//               }
//             else 
//               {
//               material_factor_sens = 0.0;
//               temp_strain_material_factor_sens = 0.0;
//               }
            
//             // also, calculate the sensitivity of the dof vector. The transformation matrix sensitivity 
//             // will be zero for this case, since shape sensitivity will be zero
//             transform_mat->right_multiply_vector(dof_value_sens, local_dof_sens);
//           }
//           break;
          
//         case SHAPE_PARAMETER_ENUM_ID:
//           {
//             Elem* pert_elem = NULL;
//             pert_elem = this->analysis_discipline.getPerturbedElemForShapeParameter(this->elem_ID,
//                                                                                     *dv_it);
//             this->reinitForShapeSensitivity(dv_ID, pert_elem, 
//                                             (*dv_it)->getPerturbationStepSize());
            
//             // set the factor_sens to zero
//             material_factor_sens = 0.0;
//             temp_strain_material_factor_sens = 0.0;
            
//             // get the transformation matrix sensitivity and calculate the local_dof_sens
//             this->getStructuralT_matrix(transform_mat_sens.get(), FESystemElem::BASE_ELEM::num(), true);
            
//             local_dof_sens.zero(); scratch_vec.zero();
            
//             transform_mat_sens->right_multiply_vector(dof_values, local_dof_sens);
//             transform_mat->right_multiply_vector(dof_value_sens, scratch_vec);
//             local_dof_sens.add(1.0, scratch_vec);
            
//             // calculate the length for the perturbed element, and find the sensitivity
//             // calculate the length of the element
//             // first init the perturbed elem
//             const Point& point0 = perturbed_elem->point(0);
//             const Point& point1 = perturbed_elem->point(1);
	  				
//             double perturbed_length = 0.0;
//             {
//               double dx=0.0, dy=0.0, dz=0.0;
//               dx = point0(0)-point1(0);
//               dy = point0(1)-point1(1);
//               dz = point0(2)-point1(2);
//               perturbed_length = pow((dx*dx+dy*dy+dz*dz),0.5);
//             }
	  				
//             perturbed_length_factor = -1.0 / (length * length) * 
//               (perturbed_length - length)/ this->perturbation;
//           }
//           break;
          
//         default:
//           abort();
//           break;
//         }
      
//       // calculate the strain and stress
//       strain_sens = length_factor * (local_dof_sens(1)-local_dof_sens(0)) + 
//         perturbed_length_factor * (local_dof(1) - local_dof(0));
//       temp_sens = (nodal_temp(0) + nodal_temp(1) ) * 0.5;
//       mechanical_strain_sens = (strain_sens - temp_strain_material_factor_sens * temp - 
//                                 temp_strain_material_factor * temp_sens);
//       stress_sens = material_factor * mechanical_strain_sens + material_factor_sens * mechanical_strain;
      
//       // add the two tensors
//       tensor.zero();
//       tensor(0,0) = strain_sens;
//       return_qty->addStrainTensor(tensor, *load_case_it, dv_ID);
      
//       tensor.zero();
//       tensor(0,0) = stress_sens;
//       return_qty->addStressTensor(tensor, *load_case_it, dv_ID);
      
//       tensor.zero();
//       tensor(0,0) = mechanical_strain_sens;
//       return_qty->addMechanicalStrainTensor(tensor, *load_case_it, dv_ID);
//       }
//     }  
  
//   // finally return the post process qty
//   return return_qty;
  
}


