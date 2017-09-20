// $Id: brick_hex8.C,v 1.24.6.1 2007-03-14 22:05:03 manav Exp $

// C++ includes


// FESystem includes
#include "StructuralElems/brick_hex8.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic3DElemDataCard.h"
#include "FESystem/AnalysisCase.h"

// libMesh includes
#include "geom/cell_hex8.h"

FESystemElem::BrickHex8::BrickHex8(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::StructuralElem(3, FESystemElem::STRUCTURAL_BRICK_HEX8::num(), discipline)
{
  
}





FESystemElem::BrickHex8::~BrickHex8()
{
	
}



void
FESystemElem::BrickHex8::calculate_K_G(DenseMatrix<double>* matrix, 
                                       const unsigned int design_point_enum_ID,
                                       bool sensitivity_calculation)
{
  // params not used here
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;

  abort();
}




void
FESystemElem::BrickHex8::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType& fe = fetypes.back();
  fe.order = FIRST;
  fe.family = LAGRANGE;
}





void
FESystemElem::BrickHex8::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
  quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                     value_type(LAGRANGE, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());
}





void
FESystemElem::BrickHex8::calculate_M(DenseMatrix<double>* matrix, 
                                     const unsigned int design_point_enum_ID,
                                     bool sensitivity_calculation)
{
  static double factor, factor_sens, factor_temp_sens;
  factor = 0.0; factor_sens = 0.0; factor_temp_sens = 0.0;
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
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
  static DenseVector<double> nodal_temp_sens(n_nodes);
  nodal_temp_sens.zero();
  static double avg_factor_temp_sens, avg_temp_sens, N_N_factor, N_N_sens_factor;
  avg_factor_temp_sens = 0.0; avg_temp_sens = 0.0;
  N_N_factor = 0.0; N_N_sens_factor = 0.0;

  static DenseMatrix<double> N_N, N_N_sens;
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
FESystemElem::BrickHex8::calculate_K(DenseMatrix<double>* matrix, 
                                     const unsigned int design_point_enum_ID,
                                     bool sensitivity_calculation)
{
  static DenseMatrix<double> material_mat(6,6), material_mat_sens(6,6),
  material_mat_temp_sens(6,6);
  material_mat.zero();   material_mat_sens.zero();   material_mat_temp_sens.zero();
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  this->elem_property_card->getFactor(material_mat, SOLID_STIFFNESS_MATRIX_FACTOR::num());
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (material_mat_sens, SOLID_STIFFNESS_MATRIX_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (material_mat_temp_sens, SOLID_STIFFNESS_MATRIX_FACTOR::num(), 
             Property::TEMPERATURE::num());
        break;
      }
    }

  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(n_nodes);
  nodal_temp_sens.zero();
  
  static DenseMatrix<double> Nx_Nx(n_nodes,n_nodes), Ny_Ny(n_nodes,n_nodes),
    Nz_Nz(n_nodes,n_nodes), Nx_Ny(n_nodes,n_nodes), Ny_Nz(n_nodes,n_nodes),
    Nz_Nx(n_nodes,n_nodes);
  static DenseMatrix<double> Nx_Nx_sens(n_nodes,n_nodes), Ny_Ny_sens(n_nodes,n_nodes),
    Nz_Nz_sens(n_nodes,n_nodes), Nx_Ny_sens(n_nodes,n_nodes), Ny_Nz_sens(n_nodes,n_nodes),
    Nz_Nx_sens(n_nodes,n_nodes);
  static DenseMatrix<double> N_N_factor(6,6), N_N_sens_factor(6,6);
  
	Nx_Nx.zero(); Ny_Ny.zero(); Nz_Nz.zero(); Nx_Ny.zero(); 
  Ny_Nz.zero(); Nz_Nz.zero();
	Nx_Nx_sens.zero(); Ny_Ny_sens.zero(); Nz_Nz_sens.zero(); Nx_Ny_sens.zero(); 
  Ny_Nz_sens.zero(); Nz_Nz_sens.zero();
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
          this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &Nz_Nz,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty(FESystemElem::N_Y_N_Z_FACTOR::num(), &Ny_Nz,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty(FESystemElem::N_Z_N_X_FACTOR::num(), &Nz_Nx,
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
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Z_N_Z_FACTOR::num(), &Nz_Nz_sens,
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny_sens, 
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Y_N_Z_FACTOR::num(), &Ny_Nz_sens,
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Z_N_X_FACTOR::num(), &Nz_Nx_sens,
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
            this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &Nz_Nz,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            this->getFESystemElemQty(FESystemElem::N_Y_N_Z_FACTOR::num(), &Ny_Nz,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
            this->getFESystemElemQty(FESystemElem::N_Z_N_X_FACTOR::num(), &Nz_Nx,
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
    this->getFESystemElemQty(FESystemElem::N_Z_N_Z_FACTOR::num(), &Nz_Nz,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_Y_N_Z_FACTOR::num(), &Ny_Nz,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_Z_N_X_FACTOR::num(), &Nz_Nx,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    }	
	
  // calculate the stiffness matrix
  // it is assumed here that all dimensions for the matrices are correct, hence, is not being
  // checked. 
  static double factor1, factor2, factor3, factor1_sens, factor2_sens, factor3_sens;
  factor1 = N_N_factor(0,0);
  factor2 = N_N_factor(0,1);
  factor3 = N_N_factor(3,3);

  factor1_sens = N_N_sens_factor(0,0);
  factor2_sens = N_N_sens_factor(0,1);
  factor3_sens = N_N_sens_factor(3,3);
  
  for (unsigned int i=0; i<n_nodes; i++)
    for (unsigned int j=0; j<n_nodes; j++)
      {
      (*matrix)(i,j) = factor1*Nx_Nx(i,j) + factor3*(Ny_Ny(i,j)+Nz_Nz(i,j)) + 
      factor1_sens*Nx_Nx_sens(i,j) + factor3_sens*(Ny_Ny_sens(i,j)+Nz_Nz_sens(i,j));
      
      (*matrix)(i,n_nodes+j) = factor2*Nx_Ny(i,j) + factor3*(Nx_Ny(j,i)) + 
        factor2_sens*Nx_Ny_sens(i,j) + factor3_sens*(Nx_Ny_sens(j,i));
      (*matrix)(i,2*n_nodes+j) = factor2*Nz_Nx(j,i) + factor3*(Nz_Nx(i,j)) + 
        factor2_sens*Nz_Nx_sens(j,i) + factor3_sens*(Nz_Nx_sens(i,j));
      
      (*matrix)(n_nodes+i,j) = factor2*Nx_Ny(j,i) + factor3*(Nx_Ny(i,j)) + 
        factor2_sens*Nx_Ny_sens(j,i) + factor3_sens*(Nx_Ny_sens(i,j));
      (*matrix)(n_nodes+i,n_nodes+j) = factor1*Ny_Ny(i,j) + factor3*(Nx_Nx(i,j)+Nz_Nz(i,j)) + 
        factor1_sens*Ny_Ny_sens(i,j) + factor3_sens*(Nx_Nx_sens(i,j)+Nz_Nz_sens(i,j));
      (*matrix)(n_nodes+i,2*n_nodes+j) = factor2*Ny_Nz(i,j) + factor3*(Ny_Nz(j,i)) + 
        factor2_sens*Ny_Nz_sens(i,j) + factor3_sens*(Ny_Nz_sens(j,i));
      
      (*matrix)(2*n_nodes+i,j) = factor2*Nz_Nx(i,j) +  factor3*(Nz_Nx(j,i)) + 
        factor2_sens*Nz_Nx_sens(i,j) +  factor3_sens*(Nz_Nx_sens(j,i));
      (*matrix)(2*n_nodes+i,n_nodes+j) = factor2*Ny_Nz(j,i) + factor3*(Ny_Nz(i,j)) + 
        factor2_sens*Ny_Nz_sens(j,i) + factor3_sens*(Ny_Nz_sens(i,j));
      (*matrix)(2*n_nodes+i,2*n_nodes+j) = factor1*Nz_Nz(i,j) + factor3*(Ny_Ny(i,j)+Nx_Nx(i,j)) + 
        factor1_sens*Nz_Nz_sens(i,j) + factor3_sens*(Ny_Ny_sens(i,j)+Nx_Nx_sens(i,j));
      }
}





void
FESystemElem::BrickHex8::calculate_F_T(DenseVector<double>* vector, 
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
    nz_n_matrix(n_nodes, n_nodes),
    nx_n_matrix_sens(n_nodes,n_nodes), 
    ny_n_matrix_sens(n_nodes,n_nodes),
    nz_n_matrix_sens(n_nodes, n_nodes);
  
  nx_n_matrix.zero();   ny_n_matrix.zero();   nz_n_matrix.zero(); 
  nx_n_matrix_sens.zero();   ny_n_matrix_sens.zero();   nz_n_matrix_sens.zero(); 

  static double avg_factor_temp_sens, avg_temp_sens;
  avg_factor_temp_sens = 0.0; avg_temp_sens = 0.0;
  
  // now, scale the matrix with the right factors to take into 
  // account the sensitivity if needed
  if (sensitivity_calculation)
    {

    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      for (unsigned int i=0; i < n_nodes; i++)
        avg_temp_sens += nodal_temp_sens(i);
      avg_temp_sens /= (1.0 * n_nodes);
      
      avg_factor_temp_sens = factor_temp_sens * avg_temp_sens;
      }
    
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
          this->getFESystemElemQty(FESystemElem::N_Z_N_FACTOR::num(), &nz_n_matrix,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          
          double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
          for (unsigned int i=0; i<n_nodes; i++)
            {
            sum1 = 0.0; sum2 = 0.0; sum3 = 0.0;
            for (unsigned int j=0; j<n_nodes; j++)
              {
              sum1 += nx_n_matrix(i,j) * 
              ((factor_sens + avg_factor_temp_sens)* (nodal_temp(j)-ref_temp) +
               factor * nodal_temp_sens(j));
              sum2 += ny_n_matrix(i,j) * 
                ((factor_sens + avg_factor_temp_sens)* (nodal_temp(j)-ref_temp) + 
                 factor * nodal_temp_sens(j));
              sum3 += nz_n_matrix(i,j) * 
                ((factor_sens + avg_factor_temp_sens)* (nodal_temp(j)-ref_temp) + 
                 factor * nodal_temp_sens(j));
              }
            (*vector)(i) = sum1;
            (*vector)(n_nodes+i) = sum2;
            (*vector)(2*n_nodes + i) = sum3;
            }
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQty( FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix,
                                    design_point_enum_ID,
                                    FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty( FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix,
                                    design_point_enum_ID,
                                    FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty( FESystemElem::N_Z_N_FACTOR::num(), &nz_n_matrix,
                                    design_point_enum_ID,
                                    FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity
            ( FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix_sens,
              FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity
            ( FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix_sens,
              FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity
            ( FESystemElem::N_Z_N_FACTOR::num(), &nz_n_matrix_sens,
              FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          
          double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
          for (unsigned int i=0; i<n_nodes; i++)
            {
            sum1 = 0.0; sum2 = 0.0; sum3= 0.0;
            for (unsigned int j=0; j<n_nodes; j++)
              {
              sum1 += factor* (nx_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) + 
                               nx_n_matrix(i,j) * nodal_temp_sens(j)) +
              avg_factor_temp_sens * nx_n_matrix(i,j) * (nodal_temp(j)-ref_temp);
              sum2 += factor* ( ny_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) + 
                                ny_n_matrix(i,j) * nodal_temp_sens(j)) +
                avg_factor_temp_sens * ny_n_matrix(i,j) * (nodal_temp(j)-ref_temp);	
              sum3 += factor* ( nz_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) +
                                nz_n_matrix(i,j) * nodal_temp_sens(j)) +
                avg_factor_temp_sens * nz_n_matrix(i,j) * (nodal_temp(j)-ref_temp);
              }
            (*vector)(i) = sum1;
            (*vector)(n_nodes + i) = sum2;
            (*vector)(2*n_nodes + i) = sum3;
            }
        }
        break;
        
      default:
        Assert(false, ExcInternalError());
      }
    }
  else
    {
    this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix, 
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix, 
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    this->getFESystemElemQty(FESystemElem::N_Z_N_FACTOR::num(), &nz_n_matrix,
                             design_point_enum_ID,
                             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
    
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    for (unsigned int i=0; i<n_nodes; i++)
      {
      sum1 = 0.0; sum2 = 0.0; sum3 = 0.0;
      for (unsigned int j=0; j<n_nodes; j++)
        {
        sum1 += nx_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
        sum2 += ny_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
        sum3 += nz_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
        }
      (*vector)(i) = sum1;
      (*vector)(n_nodes + i) = sum2;
      (*vector)(2*n_nodes + i) = sum3;
      }
    
    }
  
}



std::auto_ptr<ElemPostProcessQty> 
FESystemElem::BrickHex8::getElementPostProcessQty(std::vector<unsigned int> load_cases, 
                                                  std::vector<DesignData::DesignParameter*> dv_vector)
{
  (void) load_cases;
  (void) dv_vector;
//  // this element does not have any strain or stress. Hence, it will just return a zero vector
//  std::auto_ptr<ElemPostProcessQty> return_qty(new ElemPostProcessQty(this->elem_ID));
//  
//  TensorValue<double> tensor;
//  
//  // iterate over all the load cases, to calculate the strains and stresses
//	
//  std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
//  std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
//	
//  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
//  dv_begin = dv_vector.begin();
//  dv_end = dv_vector.end();
//  
//  // iterate over each load case, and ask solver to solve for it
//  for (; load_case_it != load_case_end; load_case_it++)
//    {
//    return_qty->addStrainTensor(tensor, *load_case_it);
//    return_qty->addStressTensor(tensor, *load_case_it);
//    
//    dv_it = dv_begin;
//    for (; dv_it != dv_end; dv_it++)
//      {
//      return_qty->addStrainTensor(tensor, *load_case_it, (*dv_it)->getID());
//      return_qty->addStressTensor(tensor, *load_case_it, (*dv_it)->getID());
//      }
//    }
//  
//  return return_qty;  
}
