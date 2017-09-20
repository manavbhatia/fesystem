// $Id: membrane.C,v 1.22.6.1 2007-03-14 22:05:03 manav Exp $

// C++ includes


// FESystem includes
#include "StructuralElems/membrane.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic2DElemDataCard.h"
#include "FESystem/AnalysisCase.h"

// libMesh includes





FESystemElem::Membrane::Membrane(const unsigned int dim,
                                 const unsigned int elem_enum_ID,
                                 Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::StructuralElem(dim, elem_enum_ID, discipline)
{
  
}





FESystemElem::Membrane::~Membrane()
{
	
}





void
FESystemElem::Membrane::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType& fe = fetypes.back();
  fe.order = FIRST;
  fe.family = LAGRANGE;
}





void
FESystemElem::Membrane::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
    quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                       value_type(LAGRANGE, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());
}



void
FESystemElem::Membrane::calculate_K_G(DenseMatrix<double>* matrix, 
                                const unsigned int design_point_enum_ID,
                                   bool sensitivity_calculation)
{
  // param not used here
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;

  abort();
}


void
FESystemElem::Membrane::calculate_M(DenseMatrix<double>* matrix, 
                                const unsigned int design_point_enum_ID,
                                    bool sensitivity_calculation)
{
  static double factor, factor_sens, factor_temp_sens;
  factor = 0.0; factor_sens = 0.0; factor_temp_sens = 0.0;
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  this->elem_property_card->getFactor(factor, MEMBRANE_MASS_FACTOR::num());
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (factor_sens, MEMBRANE_MASS_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (factor_temp_sens, MEMBRANE_MASS_FACTOR::num(), Property::TEMPERATURE::num());
        break;
      }
    }


  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(n_nodes);
  nodal_temp_sens.zero();
  static double avg_factor_temp_sens, avg_temp_sens, N_N_factor, N_N_sens_factor;
  avg_factor_temp_sens = 0.0; avg_temp_sens = 0.0;
  N_N_factor = 0.0; N_N_sens_factor = 0.0;

  static DenseMatrix<double> N_N(n_nodes, n_nodes), N_N_sens(n_nodes, n_nodes);
  N_N.zero(); N_N_sens.zero();
	
  if (sensitivity_calculation)
    {
    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
      for (unsigned int i=0; i < n_nodes; i++)
        avg_temp_sens += nodal_temp_sens(i);
      avg_temp_sens /= (1.0 * n_nodes);
      
      avg_factor_temp_sens = factor_temp_sens * avg_temp_sens;
      }
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          
          N_N_factor = factor_sens + avg_factor_temp_sens;
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_N_FACTOR::num(), &N_N_sens,
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          
          N_N_sens_factor = factor;
          
          if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
            {
            this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            
            N_N_factor = avg_factor_temp_sens;
            }
        }
        break;
        
      default:
        Assert(false, ExcInternalError());
      }
    }
  else
    {
    this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    
    N_N_factor = factor;
    }
  
  
  // now perform the additions
  N_N.scale(N_N_factor);
  N_N_sens.scale(N_N_sens_factor);
  
  // calculate the stiffness matrix
  // it is assumed here that rest of the code works fine, and so, 
  // m == n == n_nodes. Otherwise, it should be checked.
  for (unsigned int i=0; i<n_nodes; i++)
    for (unsigned int j=0; j<n_nodes; j++)
      {
      (*matrix)(i,j) = N_N(i,j) + N_N_sens(i,j);
      (*matrix)(n_nodes+i,n_nodes+j) = N_N(i,j) + N_N_sens(i,j);
      (*matrix)(2*n_nodes+i,2*n_nodes+j) = N_N(i,j) + N_N_sens(i,j);
      }
}







void
FESystemElem::Membrane::calculate_K(DenseMatrix<double>* matrix, 
                                    const unsigned int design_point_enum_ID,
                                    bool sensitivity_calculation)
{
  
  static DenseMatrix<double> material_mat(3,3), material_mat_sens(3,3),
  material_mat_temp_sens(3,3);
  material_mat.zero();   material_mat_sens.zero();   material_mat_temp_sens.zero();
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  this->elem_property_card->getFactor(material_mat, STIFFNESS_A_MATRIX_FACTOR::num());
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (material_mat_sens, STIFFNESS_A_MATRIX_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (material_mat_temp_sens, STIFFNESS_A_MATRIX_FACTOR::num(), 
             Property::TEMPERATURE::num());
        break;
      }
    }
   
  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(n_nodes);
  nodal_temp_sens.zero();
  
  static DenseMatrix<double> Nx_Nx(n_nodes,n_nodes), Ny_Ny(n_nodes,n_nodes),
    Nx_Ny(n_nodes,n_nodes);
  static DenseMatrix<double> Nx_Nx_sens(n_nodes,n_nodes), Ny_Ny_sens(n_nodes,n_nodes),
    Nx_Ny_sens(n_nodes,n_nodes);
  static DenseMatrix<double> N_N_factor(3,3), N_N_sens_factor(3,3);
  
	Nx_Nx.zero(); Ny_Ny.zero(); Nx_Ny.zero(); 
	Nx_Nx_sens.zero(); Ny_Ny_sens.zero(); Nx_Ny_sens.zero(); 
	N_N_factor.zero(), N_N_sens_factor.zero();
	
  if (sensitivity_calculation)
    {
    
    static double avg_temp_sens = 0.0;
    avg_temp_sens = 0.0;
    
    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
      for (unsigned int i=0; i < n_nodes; i++)
        avg_temp_sens += nodal_temp_sens(i);
      avg_temp_sens /= (1.0 * n_nodes);
      
      material_mat_temp_sens.scale(avg_temp_sens);
      }
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);

          N_N_factor = material_mat_sens;
          N_N_factor += material_mat_temp_sens;
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx_sens, 
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny_sens, 
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny_sens, 
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          
          N_N_sens_factor = material_mat;
          
          if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
            {
            this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            
            N_N_factor = material_mat_temp_sens;
            }
        }
        break;
        
      default:
        Assert(false, ExcInternalError());
      }
    }
  else
    {
    this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    N_N_factor = material_mat;
    }
	
  // calculate the stiffness matrix
  // it is assumed here that all dimensions for the matrices are correct, hence, is not being
  // checked. 
  static double factor1, factor2, factor3, factor1_sens, factor2_sens, factor3_sens;
  factor1 = N_N_factor(0,0);
  factor2 = N_N_factor(0,1);
  factor3 = N_N_factor(2,2);
  
  factor1_sens = N_N_sens_factor(0,0);
  factor2_sens = N_N_sens_factor(0,1);
  factor3_sens = N_N_sens_factor(2,2);
  
  for (unsigned int i=0; i<n_nodes; i++)
    for (unsigned int j=0; j<n_nodes; j++)
      {
      (*matrix)(i,j) = factor1*Nx_Nx(i,j) + factor3 * Ny_Ny(i,j) + 
      factor1_sens*Nx_Nx_sens(i,j) + factor3_sens * Ny_Ny_sens(i,j);
      
      (*matrix)(i,n_nodes+j) = factor2*Nx_Ny(i,j) + factor3 * Nx_Ny(j,i) + 
        factor2_sens*Nx_Ny_sens(i,j) + factor3_sens* Nx_Ny_sens(j,i);
      
      (*matrix)(n_nodes+i,j) = factor2*Nx_Ny(j,i) + factor3*(Nx_Ny(i,j)) + 
        factor2_sens*Nx_Ny_sens(j,i) + factor3_sens*Nx_Ny_sens(i,j);

      (*matrix)(n_nodes+i,n_nodes+j) = factor1*Ny_Ny(i,j) + factor3*Nx_Nx(i,j) + 
        factor1_sens*Ny_Ny_sens(i,j) + factor3_sens*Nx_Nx_sens(i,j);
      }
}








void
FESystemElem::Membrane::calculate_F_T(DenseVector<double>* vector, 
                                      const unsigned int design_point_enum_ID,
                                      bool sensitivity_calculation)
{
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
  
  // get the loads for this element
  this->extractNodalTemperatureVectorFromLoads(nodal_temp, false);
  
  static double factor, factor_sens, factor_temp_sens;
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
  
  static DenseMatrix<double> nx_n_matrix(n_nodes,n_nodes), 
    ny_n_matrix(n_nodes,n_nodes), 
    nx_n_matrix_sens(n_nodes,n_nodes), 
    ny_n_matrix_sens(n_nodes,n_nodes);
  
  nx_n_matrix.zero(); ny_n_matrix.zero(); 
  nx_n_matrix_sens.zero(); ny_n_matrix_sens.zero(); 

  static double avg_factor_temp_sens, avg_temp_sens;
  avg_factor_temp_sens = 0.0; avg_temp_sens = 0.0;

  // now, scale the matrix with the right factors to take into 
  // account the sensitivity if needed
  switch (sensitivity_calculation)
    {

    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      for (unsigned int i=0; i < n_nodes; i++)
        avg_temp_sens += nodal_temp_sens(i);
      avg_temp_sens /= (1.0 * n_nodes);
      
      avg_factor_temp_sens = factor_temp_sens * avg_temp_sens;
      }
    
    case true:
      {
        switch (this->sensitivity_parameter)
          {
          case PROPERTY_PARAMETER_ENUM_ID:
            {
              this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              
              double sum1 = 0.0, sum2 = 0.0;
              for (unsigned int i=0; i<n_nodes; i++)
                {
                sum1 = 0.0; sum2 = 0.0;
                for (unsigned int j=0; j<n_nodes; j++)
                  {
                  sum1 += nx_n_matrix(i,j) * 
                  ((factor_sens + avg_factor_temp_sens)* (nodal_temp(j)-ref_temp) +
                   factor * nodal_temp_sens(j));
                  sum2 += ny_n_matrix(i,j) * 
                    ((factor_sens + avg_factor_temp_sens)* (nodal_temp(j)-ref_temp) + 
                     factor * nodal_temp_sens(j));
                  }
                (*vector)(i) = sum1;
                (*vector)(n_nodes+i) = sum2;
                }
            }
            break;
            
          case SHAPE_PARAMETER_ENUM_ID:
            {
              this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_FACTOR::num(),
                                                       &nx_n_matrix_sens,
                                                       FESystemElem::ELEM_VOLUME::num(),
                                                       LAGRANGE);
              this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Y_N_FACTOR::num(), 
                                                       &ny_n_matrix_sens,
                                                       FESystemElem::ELEM_VOLUME::num(),
                                                       LAGRANGE);
              
              double sum1 = 0.0, sum2 = 0.0;
              for (unsigned int i=0; i<n_nodes; i++)
                {
                sum1 = 0.0; sum2 = 0.0; 
                for (unsigned int j=0; j<n_nodes; j++)
                  {
                  sum1 += factor* (nx_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) + 
                                   nx_n_matrix(i,j) * nodal_temp_sens(j)) +
                  avg_factor_temp_sens * nx_n_matrix(i,j) * (nodal_temp(j)-ref_temp);
                  sum2 += factor* ( ny_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) + 
                                    ny_n_matrix(i,j) * nodal_temp_sens(j)) +
                    avg_factor_temp_sens * ny_n_matrix(i,j) * (nodal_temp(j)-ref_temp);	
                  }
                (*vector)(i) = sum1;
                (*vector)(n_nodes + i) = sum2;
                }
            }
            break;
            
          default:
            Assert(false, ExcInternalError());
          }
      }
      break;
      
    case false:
    default:
      {
        this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix, 
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix, 
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        
        double sum1 = 0.0, sum2 = 0.0;
        for (unsigned int i=0; i<n_nodes; i++)
          {
          sum1 = 0.0; sum2 = 0.0;
          for (unsigned int j=0; j<n_nodes; j++)
            {
            sum1 += nx_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
            sum2 += ny_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
            }
          (*vector)(i) = sum1;
          (*vector)(n_nodes + i) = sum2;
          }
      }
      break;
    }
}



std::auto_ptr<ElemPostProcessQty> 
FESystemElem::Membrane::getElementPostProcessQty(std::vector<unsigned int> load_cases, 
                                                 std::vector<DesignData::DesignParameter*> dv_vector)
{
  (void) load_cases;
  (void) dv_vector;

//   std::auto_ptr<ElemPostProcessQty> return_qty(new ElemPostProcessQty(this->elem_ID));
  
//   TensorValue<double> tensor;
  
//   unsigned int n_nodes = this->base_elem->n_nodes();
  
//   static DenseMatrix<double> B_strain(3, 6*n_nodes) , B_strain_sens(3, 6*n_nodes);
//   B_strain.zero(); B_strain_sens.zero();
  
//   // create rest of the necessary data structures
//   static DenseVector<double>  strain(3), strain_sens(3),
//     stress(3), stress_sens(3), strain_plus_thermal(3), 
//     strain_plus_thermal_sens(3), strain_scratch(3);
  
//   strain.zero(); strain_sens.zero();
//   stress.zero(); stress_sens.zero(); strain_plus_thermal.zero(); 
//   strain_plus_thermal_sens.zero(); strain_scratch.zero();
  
//   static DenseMatrix<double> material_factor(3,3), material_factor_sens(3,3);
//   material_factor.zero(); 
//   material_factor_sens.zero();
  
//   double temp_strain_material_factor = 0.0, temp_strain_material_factor_sens = 0.0,
//     ref_temp = 0.0;
    
    
//     // point at which this element will be initialized
//     // for now, this is only the centroid of the element
//     std::vector<Point> point_vec;
//     point_vec.push_back(Point(0.0,0.0,0.0));
    
//     // now initialize the element at the centroid of the element where all the 
//     // data structures will be calculated
//     this->initialize_element(&point_vec);
    
//     // get the necessary data out of this element 
//     const std::vector< std::vector<Real> >& dphi_dx = this->fe_base_map[LAGRANGE]->get_dphidx();
//     const std::vector< std::vector<Real> >& dphi_dy = this->fe_base_map[LAGRANGE]->get_dphidy();
//     const std::vector<std::vector<Real> >& phi = this->fe_base_map[LAGRANGE]->get_phi();
    
//     // calculate B_strain matrix. Since the calculation is being performed at a single point, 
//     // the vector of shape factor derivatives should have size = 1
//     assert (dphi_dx[0].size() == 1);
//     for (unsigned int i=0; i < n_nodes; i++)
//       {
//       B_strain(0,i) = dphi_dx[i][0];  // epsilon_xx
//       B_strain(1,n_nodes+i) = dphi_dy[i][0];  // epsilon_yy
//       B_strain(2,i) = dphi_dy[i][0];  // epsilon  xy 
//       B_strain(2,n_nodes+i) = dphi_dx[i][0];  // epsilon xy
//       }
    
    
//     // iterate over all the load cases, to calculate the strains and stresses
//     std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
//     std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
    
//     std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
//     dv_begin = dv_vector.begin();
//     dv_end = dv_vector.end();
    
    
//     static DenseMatrix<double> transform_mat, transform_mat_sens;
//     transform_mat.zero(); transform_mat_sens.zero();
    
//     this->getStructuralT_matrix(&transform_mat, false);
    
    
//     // also, create dof_vectors
//     static DenseVector<double>  dof_values(6*n_nodes), dof_value_sens(6*n_nodes),
//       local_dof(6*n_nodes), local_dof_sens(6*n_nodes), scratch_vec(6*n_nodes),
//       nodal_temp(n_nodes),nodal_temp_sens(n_nodes),
//       dphi_dx_sens(n_nodes), dphi_dy_sens(n_nodes);
    
//     dof_values.zero(); dof_value_sens.zero();
//     local_dof.zero(); local_dof_sens.zero(); 
//     nodal_temp.zero(); nodal_temp_sens.zero();
//     dphi_dx_sens.zero(); dphi_dy_sens.zero(); 
//     scratch_vec.zero();
    
//     // initialize the material factors
//     this->elem_property_card->getPropertyValueFromMaterialCard(ALPHA_EXPANSION::num(),
//                                                                temp_strain_material_factor);
//     this->elem_property_card->getPropertyValueFromMaterialCard(TEMP_REF::num(),
//                                                                ref_temp);

//     this->elem_property_card->getFactor(material_factor, STRESS_STRAIN_FACTOR::num());

//     double temp = 0.0, temp_sens = 0.0;
    
//     // iterate over each load case, and ask solver to solve for it
//     for (; load_case_it != load_case_end; load_case_it++)
//       {
      
//       // get the DOF vector for this elem
//       this->analysis_discipline.getElemDofValues(this->base_elem, 
//                                                  dof_values,
//                                                  *load_case_it);
      
//       // transform these dofs to the local coordinate system
//       local_dof.zero();
//       transform_mat.right_multiply_vector(dof_values, local_dof);
      
//       // get the loads for this element
//       this->getNodalTemperatureVector(nodal_temp, *load_case_it);
      
//       // calculate the factor and the strain and stress
//       B_strain.right_multiply_vector(local_dof, strain);
      
//       // add the thermal strain component and then calculate the stress
//       for (unsigned int i=0; i<phi.size(); i++)
//         temp += phi[i][0] * nodal_temp(i);
      
//       strain_plus_thermal.zero();
//       strain_plus_thermal(0) = strain(0) - temp_strain_material_factor * (temp - ref_temp);
//       strain_plus_thermal(1) = strain(1) - temp_strain_material_factor * (temp - ref_temp);
//       strain_plus_thermal(2) = strain(2);
      
//       // now multiply to calculate the stresses
//       material_factor.right_multiply_vector(strain_plus_thermal, stress);
      
      
//       // add the two tensors
//       tensor.zero();
//       tensor(0,0) = strain(0);  // xx
//       tensor(1,1) = strain(1);  // yy
//       tensor(0,1) = strain(2);  // xy
//       tensor(1,0) = strain(2);  // yx
//       return_qty->addStrainTensor(tensor, *load_case_it);
      
//       tensor.zero();
//       tensor(0,0) = stress(0);  // xx
//       tensor(1,1) = stress(1);  // yy
//       tensor(0,1) = stress(2);  // xy
//       tensor(1,0) = stress(2);  // yx
//       return_qty->addStressTensor(tensor, *load_case_it);
      
//       tensor.zero();
//       tensor(0,0) = strain_plus_thermal(0);  // xx
//       tensor(1,1) = strain_plus_thermal(1);  // yy
//       tensor(0,1) = strain_plus_thermal(2);  // xy
//       tensor(1,0) = strain_plus_thermal(2);  // yx
//       return_qty->addMechanicalStrainTensor(tensor, *load_case_it);
      
      
//       unsigned int dv_ID = 0;
//       // now calculate the sensitivities
//       dv_it = dv_begin;
//       for (; dv_it != dv_end; dv_it++)
//         {
//         dv_ID = (*dv_it)->getID();
//         this->clearSensitivityInitialization();
        
//         // get the dof value sensitivity for this case
//         this->analysis_discipline.getElemDofValues(this->base_elem,
//                                                    dof_value_sens,
//                                                    *load_case_it,
//                                                    true,
//                                                    dv_ID);
        
//         // get the load sensitivity vectors
//         this->getNodalTemperatureVector(nodal_temp_sens, *load_case_it, true, dv_ID);
        
//         material_factor_sens.zero();
//         temp_strain_material_factor_sens = 0.0;
        
//         switch ((*dv_it)->getParameterTypeEnumID())
//           {
//           case PROPERTY_PARAMETER_ENUM_ID:
//             {
//               this->reinitForPropertySensitivity(dv_ID);

//               // for a property DV, the strain operator sensitivity will be zero
//               B_strain_sens.zero();
              
//               // if the property ID is the same as this elements property, then set the 
//               // value of the property sensitivity. Else the value be zero
//               if ( this->elem_property_card->checkGlobalParameterDependence(dv_ID))
//                 {
//                 this->elem_property_card->getFactorSensitivityForGlobalParameter
//                 (material_factor_sens, STRESS_STRAIN_FACTOR::num(), dv_ID);
//                 this->elem_property_card->
//                 getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//                 (ALPHA_EXPANSION::num(), dv_ID, temp_strain_material_factor_sens);
//                 }
//               else 
//                 {
//                 material_factor_sens.zero();
                
//                 temp_strain_material_factor_sens = 0.0;
//                 }
              
//               // also, calculate the sensitivity of the dof vector. The transformation matrix sensitivity 
//               // will be zero for this case, since shape sensitivity will be zero
//               transform_mat.right_multiply_vector(dof_value_sens, local_dof_sens);
//             }
//             break;
            
//           case SHAPE_PARAMETER_ENUM_ID:
//             {
//               Elem* pert_elem = NULL;
//               pert_elem = this->analysis_discipline.getPerturbedElemForShapeParameter(this->elem_ID,
//                                                                                       *dv_it);
//               this->reinitForShapeSensitivity(dv_ID, pert_elem, 
//                                               (*dv_it)->getPerturbationStepSize());
              
//               // set the factor_sens to zero
//               material_factor_sens.zero();
//               temp_strain_material_factor_sens = 0.0;
              
//               // get the transformation matrix sensitivity and calculate the local_dof_sens
//               this->getStructuralT_matrix(&transform_mat_sens, true);
              
//               local_dof_sens.zero(); scratch_vec.zero();
//               transform_mat_sens.right_multiply_vector(dof_values, local_dof_sens);
              
//               scratch_vec.zero();
//               transform_mat.right_multiply_vector(dof_value_sens, scratch_vec);
              
//               local_dof_sens.add(1.0, scratch_vec);

//               // now init the fe_sens to the perturbed DV and calculate the shape sensitivity of the
//               // B matrix
//               // then, reinit it at the base elem, at the element centroid, where the 
//               // element strains and stresses will be evaluated
//               this->attachNewElem(this->perturbed_elem);
//               this->initialize_element(&point_vec);
              
//               // get the necessary data out of this element 
//               const std::vector< std::vector<Real> >& dphi_dx_perturbed = 
//                 this->fe_base_map[LAGRANGE]->get_dphidx();
//               const std::vector< std::vector<Real> >& dphi_dy_perturbed = 
//                 this->fe_base_map[LAGRANGE]->get_dphidy();
              
//               dphi_dx_sens.zero();
//               dphi_dy_sens.zero();
              
//               for (unsigned int i=0; i<n_nodes; i++)
//                 {
//                 dphi_dx_sens(i) = dphi_dx_perturbed[i][0] - dphi_dx[i][0];
//                 dphi_dy_sens(i) = dphi_dy_perturbed[i][0] - dphi_dy[i][0];
//                 }
              
//               dphi_dx_sens.scale(1.0/this->perturbation);
//               dphi_dy_sens.scale(1.0/this->perturbation);
              
//               // calculate B_strain matrix. Since the calculation is being performed at a single point, 
//               // the vector of shape factor derivatives should have size = 1
              
//               for (unsigned int i=0; i < n_nodes; i++)
//                 {
//                 B_strain_sens(0,i) = dphi_dx_sens(i);  // epsilon_xx
//                 B_strain_sens(1,n_nodes+i) = dphi_dy_sens(i);  // epsilon_yy
//                 B_strain_sens(2,i) = dphi_dy_sens(i);  // epsilon  xy 
//                 B_strain_sens(2,n_nodes+i) = dphi_dx_sens(i);  // epsilon xy
//                 }
              
//             }
//             break;
            
//           default:
//             abort();
//             break;
//           }
        
//         // calculate the strain and stress
//         // calculate the factor and the strain and stress
//         strain_scratch.zero();
//         B_strain.right_multiply_vector(local_dof_sens, strain_sens);
//         B_strain_sens.right_multiply_vector(local_dof, strain_scratch);
//         strain_sens.add(1.0, strain_scratch);
        
//         // add the thermal strain component and then calculate the stress
//         temp_sens = 0.0;
//         for (unsigned int i=0; i<phi.size(); i++)
//           temp_sens += phi[i][0] * nodal_temp_sens(i);
        
//         // add the temperture sensitivity effects to the strain sensitivity
//         strain_plus_thermal_sens.zero();
//         strain_plus_thermal_sens(0) = strain_sens(0) - temp_strain_material_factor * temp_sens -
//           temp_strain_material_factor_sens * (temp - ref_temp);
//         strain_plus_thermal_sens(1) = strain_sens(1) - temp_strain_material_factor * temp_sens -
//           temp_strain_material_factor_sens * (temp - ref_temp);
//         strain_plus_thermal_sens(2) = strain_sens(2);
        
//         // now multiply to calculate the stresses
//         strain_scratch.zero();
//         material_factor_sens.right_multiply_vector(strain_plus_thermal, stress_sens);
//         material_factor.right_multiply_vector(strain_plus_thermal_sens, strain_scratch);
//         stress_sens.add(1.0, strain_scratch);
        
//         // add the two tensors
//         tensor.zero();
//         tensor(0,0) = strain_sens(0);  // xx
//         tensor(1,1) = strain_sens(1);  // yy
//         tensor(0,1) = strain_sens(2);  // xy
//         tensor(1,0) = strain_sens(2);  // yx
//         return_qty->addStrainTensor(tensor, *load_case_it, dv_ID);
        
//         tensor.zero();
//         tensor(0,0) = stress_sens(0);  // xx
//         tensor(1,1) = stress_sens(1);  // yy
//         tensor(0,1) = stress_sens(2);  // xy
//         tensor(1,0) = stress_sens(2);  // yx
//         return_qty->addStressTensor(tensor, *load_case_it, dv_ID);
        
//         tensor.zero();
//         tensor(0,0) = strain_plus_thermal_sens(0);  // xx
//         tensor(1,1) = strain_plus_thermal_sens(1);  // yy
//         tensor(0,1) = strain_plus_thermal_sens(2);  // xy
//         tensor(1,0) = strain_plus_thermal_sens(2);  // yx
//         return_qty->addMechanicalStrainTensor(tensor, *load_case_it,  dv_ID);
        
//         }
//       }  
    
//     return return_qty;
}
