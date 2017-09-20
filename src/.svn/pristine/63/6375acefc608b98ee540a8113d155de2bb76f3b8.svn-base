// $Id: PlateDKBatoz.C,v 1.9.6.1 2007-03-14 22:05:03 manav Exp $

// C++ includes

// FESystem includes
#include "StructuralElems/PlateDKBatoz.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic2DElemDataCard.h"
#include "FESystem/AnalysisCase.h"


// libMesh includes



FESystemElem::PlateDKBatoz::PlateDKBatoz(const unsigned int dim, 
                                         const unsigned int elem_enum_ID,
                                         Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::StructuralElem(dim, elem_enum_ID, discipline)
{
  
}





FESystemElem::PlateDKBatoz::~PlateDKBatoz()
{
	
}




void
FESystemElem::PlateDKBatoz::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType* fe = &(fetypes.back());
  fe->order = FIRST;
  fe->family = LAGRANGE;
  
  fetypes.push_back(FEType());
  fe = &(fetypes.back());
  fe->order = SECOND;
  fe->family = BCIZ;

  fetypes.push_back(FEType());
  fe = &(fetypes.back());
  fe->order = SECOND;
  fe->family = BATOZ;
}





void
FESystemElem::PlateDKBatoz::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
    quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                       value_type(LAGRANGE, std::make_pair(QGAUSS, NINTH))).second;
  
  Assert(insert_success, ExcInternalError());
  
  insert_success = 
  quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                     value_type(BCIZ, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());

  insert_success = 
    quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                       value_type(BATOZ, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());

}




void
FESystemElem::PlateDKBatoz::calculate_M(DenseMatrix<double>* matrix, 
                                        const unsigned int design_point_enum_ID,
                                      bool sensitivity_calculation)
{
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static double mem_factor, mem_factor_sens, mem_factor_temp_sens,
    plate_factor, plate_factor_sens, plate_factor_temp_sens;
  mem_factor = 0.0; mem_factor_sens = 0.0; mem_factor_temp_sens = 0.0;
  plate_factor = 0.0; plate_factor_sens = 0.0; plate_factor_temp_sens = 0.0;  

  this->elem_property_card->getFactor(mem_factor, MEMBRANE_MASS_FACTOR::num());
  this->elem_property_card->getFactor(plate_factor, PLATE_MASS_FACTOR::num());

  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          this->elem_property_card->getFactorSensitivityForGlobalParameter
          (mem_factor_sens, MEMBRANE_MASS_FACTOR::num(), this->sensitivity_parameter_ID);
          this->elem_property_card->getFactorSensitivityForGlobalParameter
            (plate_factor_sens, PLATE_MASS_FACTOR::num(), this->sensitivity_parameter_ID);
        }
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          {
          this->elem_property_card->getFactorSensitivityForLocalParameter
          (mem_factor_temp_sens, MEMBRANE_MASS_FACTOR::num(), Property::TEMPERATURE::num());
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (plate_factor_temp_sens, PLATE_MASS_FACTOR::num(), Property::TEMPERATURE::num());
          }
        break;
      }
    }


  
  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(n_nodes);
  nodal_temp_sens.zero();
  
  static double avg_mem_factor_temp_sens, avg_temp_sens, N_N_mem_factor, N_N_mem_sens_factor,
    avg_plate_factor_temp_sens, N_N_plate_factor, N_N_plate_sens_factor;

  avg_mem_factor_temp_sens = 0.0; avg_mem_factor_temp_sens = 0.0; avg_temp_sens = 0.0;
  avg_plate_factor_temp_sens = 0.0; avg_plate_factor_temp_sens = 0.0;

  N_N_mem_factor = 0.0; N_N_mem_sens_factor = 0.0;
  N_N_plate_factor = 0.0; N_N_plate_sens_factor = 0.0;
  
  static DenseMatrix<double> N_N(n_nodes, n_nodes), N_N_sens(n_nodes, n_nodes), 
  N_N_bend(n_nodes*3,n_nodes*3), N_N_bend_sens(n_nodes*3,n_nodes*3);
  N_N.zero(); N_N_sens.zero();
  N_N_bend.zero(); N_N_bend_sens.zero();
  
  
  if (sensitivity_calculation)
    {
    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
      for (unsigned int i=0; i < n_nodes; i++)
        avg_temp_sens += nodal_temp_sens(i);
      avg_temp_sens /= (1.0 * n_nodes);
      
      avg_mem_factor_temp_sens = mem_factor_temp_sens * avg_temp_sens;
      avg_plate_factor_temp_sens = plate_factor_temp_sens * avg_temp_sens;
      }
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N_bend,
                                   design_point_enum_ID,
                                   FESystemElem::ELEM_VOLUME::num(), BCIZ);
          
          N_N_mem_factor = mem_factor_sens + avg_mem_factor_temp_sens;
          N_N_plate_factor = plate_factor_sens + avg_plate_factor_temp_sens;
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_N_FACTOR::num(), &N_N_sens,
                                                   FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_N_FACTOR::num(), &N_N_bend_sens,
                                                   FESystemElem::ELEM_VOLUME::num(), BCIZ);
          
          N_N_mem_sens_factor = mem_factor;
          N_N_plate_sens_factor = plate_factor;
          
          if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
            {
            this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N,
                                     design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N_bend,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), BCIZ);
              
            N_N_mem_factor = avg_mem_factor_temp_sens;
            N_N_plate_factor = avg_plate_factor_temp_sens;
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
      this->getFESystemElemQty(FESystemElem::N_N_FACTOR::num(), &N_N_bend,
                               design_point_enum_ID,
                               FESystemElem::ELEM_VOLUME::num(), BCIZ);
      
    N_N_mem_factor = mem_factor;
    N_N_plate_factor = plate_factor;
    }
  
  
  // calculate the mass matrix
  static double val1;
  for (unsigned int i=0; i<n_nodes; i++)
    for (unsigned int j=0; j<n_nodes; j++)
      {
        val1 = N_N_mem_factor * N_N(i,j) + N_N_mem_sens_factor * N_N_sens(i,j);
        (*matrix)(i,j) = val1;
        (*matrix)(n_nodes+i,n_nodes+j) = val1;
      }

  //std::cout << _N_bend;
  
  for (unsigned int i=0; i<3*n_nodes; i++)
    for (unsigned int j=0; j<3*n_nodes; j++)
      (*matrix)(2*n_nodes+i,2*n_nodes+j) = N_N_mem_factor * N_N_bend(i,j) + N_N_mem_sens_factor * N_N_bend_sens(i,j);
  
  for (unsigned int i=0; i<n_nodes;i++)
    (*matrix)(5*n_nodes+i,5*n_nodes+i) = 1.0e-6;
}




void
FESystemElem::PlateDKBatoz::calculate_K_G(DenseMatrix<double>* matrix, 
                                          const unsigned int design_point_enum_ID,
                                   bool sensitivity_calculation)
{
  // params not used here
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;

  Assert(false, ExcInternalError());
}



void
FESystemElem::PlateDKBatoz::calculate_K(DenseMatrix<double>* matrix, 
                                        const unsigned int design_point_enum_ID,
                                      bool sensitivity_calculation)
{
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  static DenseMatrix<double> mem_material_mat(3,3), mem_material_mat_sens(3,3),
    mem_material_mat_temp_sens(3,3);
  static DenseMatrix<double> plate_material_mat(3,3), plate_material_mat_sens(3,3),
    plate_material_mat_temp_sens(3,3);
  mem_material_mat.zero(); mem_material_mat_sens.zero(); mem_material_mat_temp_sens.zero();
  plate_material_mat.zero(); plate_material_mat_sens.zero(); plate_material_mat_temp_sens.zero();
  

  this->elem_property_card->getFactor(mem_material_mat, STIFFNESS_A_MATRIX_FACTOR::num());
  this->elem_property_card->getFactor(plate_material_mat, STIFFNESS_D_MATRIX_FACTOR::num());
  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
        this->elem_property_card->getFactorSensitivityForGlobalParameter
          (mem_material_mat_sens, STIFFNESS_A_MATRIX_FACTOR::num(), 
           this->sensitivity_parameter_ID);
        this->elem_property_card->getFactorSensitivityForGlobalParameter
          (plate_material_mat_sens, STIFFNESS_D_MATRIX_FACTOR::num(), 
           this->sensitivity_parameter_ID);
        }
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          {
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (mem_material_mat_temp_sens, STIFFNESS_A_MATRIX_FACTOR::num(), 
             Property::TEMPERATURE::num());
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (plate_material_mat_temp_sens, STIFFNESS_D_MATRIX_FACTOR::num(), 
             Property::TEMPERATURE::num());
          }
        break;
      }
    }

  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(n_nodes);
  nodal_temp_sens.zero();

  // following factors are needed for the membrane effects
  static DenseMatrix<double> Nx_Nx(n_nodes,n_nodes), Ny_Ny(n_nodes,n_nodes),
    Nx_Ny(n_nodes,n_nodes);
  static DenseMatrix<double> Nx_Nx_sens(n_nodes,n_nodes), Ny_Ny_sens(n_nodes,n_nodes),
    Nx_Ny_sens(n_nodes,n_nodes);
  static DenseMatrix<double> N_N_mem_factor(3,3), N_N_sens_mem_factor(3,3),
    N_N_plate_factor(3,3), N_N_sens_plate_factor(3,3);
  
	Nx_Nx.zero(); Ny_Ny.zero(); Nx_Ny.zero(); 
	Nx_Nx_sens.zero(); Ny_Ny_sens.zero(); Nx_Ny_sens.zero(); 
  N_N_mem_factor.zero(); N_N_sens_mem_factor.zero();
  N_N_plate_factor.zero(); N_N_sens_plate_factor.zero();
  
  // following factors are needed for plate bending effects
  static DenseMatrix<double> Hxx_Hxx(n_nodes*3, n_nodes*3), Hxx_Hyy(n_nodes*3, n_nodes*3),
    Hyy_Hyy(n_nodes*3, n_nodes*3), Hxy_Hxy(n_nodes*3, n_nodes*3), 
    Hyx_Hxy(n_nodes*3, n_nodes*3), Hyx_Hyx(n_nodes*3, n_nodes*3);
  static DenseMatrix<double> Hxx_Hxx_sens(n_nodes*3, n_nodes*3), 
    Hxx_Hyy_sens(n_nodes*3, n_nodes*3),
    Hyy_Hyy_sens(n_nodes*3, n_nodes*3), Hxy_Hxy_sens(n_nodes*3, n_nodes*3), 
    Hyx_Hxy_sens(n_nodes*3, n_nodes*3), Hyx_Hyx_sens(n_nodes*3, n_nodes*3);
  Hxx_Hxx.zero(); Hxx_Hyy.zero(); Hyy_Hyy.zero(); Hxy_Hxy.zero();
  Hyx_Hxy.zero(); Hyx_Hyx.zero();
  Hxx_Hxx_sens.zero(); Hxx_Hyy_sens.zero(); Hyy_Hyy_sens.zero(); Hxy_Hxy_sens.zero();
  Hyx_Hxy_sens.zero(); Hyx_Hyx_sens.zero();
	
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
      
      mem_material_mat_temp_sens.scale(avg_temp_sens);
      plate_material_mat_temp_sens.scale(avg_temp_sens);
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
          
          FESystemElem::getPlateDKBatozQty(this,
                                           FESystemElem::PLATE_DK_BATOZ_HXX_HXX_FACTOR::num(), 
                                           &Hxx_Hxx, design_point_enum_ID,
                                           FESystemElem::ELEM_VOLUME::num(), false);
          FESystemElem::getPlateDKBatozQty(this,
                                           FESystemElem::PLATE_DK_BATOZ_HYY_HYY_FACTOR::num(), 
                                           &Hyy_Hyy, design_point_enum_ID,
                                           FESystemElem::ELEM_VOLUME::num(), false);
          FESystemElem::getPlateDKBatozQty(this,
                                           FESystemElem::PLATE_DK_BATOZ_HXX_HYY_FACTOR::num(), 
                                           &Hxx_Hyy, design_point_enum_ID,
                                           FESystemElem::ELEM_VOLUME::num(), false);
          FESystemElem::getPlateDKBatozQty(this,
                                           FESystemElem::PLATE_DK_BATOZ_HXY_HXY_FACTOR::num(), 
                                           &Hxy_Hxy, design_point_enum_ID,
                                           FESystemElem::ELEM_VOLUME::num(), false);
          FESystemElem::getPlateDKBatozQty(this,
                                           FESystemElem::PLATE_DK_BATOZ_HYX_HXY_FACTOR::num(), 
                                           &Hyx_Hxy, design_point_enum_ID,
                                           FESystemElem::ELEM_VOLUME::num(), false);
          FESystemElem::getPlateDKBatozQty(this,
                                           FESystemElem::PLATE_DK_BATOZ_HYX_HYX_FACTOR::num(), 
                                           &Hyx_Hyx, design_point_enum_ID,
                                           FESystemElem::ELEM_VOLUME::num(), false);
          
          N_N_mem_factor = mem_material_mat_sens;
          N_N_mem_factor += mem_material_mat_temp_sens;

          N_N_plate_factor = plate_material_mat_sens;
          N_N_plate_factor += plate_material_mat_temp_sens;
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          this->getFESystemElemQtyShapeSensitivity
          (FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx_sens,
           FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity
            (FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny_sens,
             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          this->getFESystemElemQtyShapeSensitivity
            (FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny_sens,
             FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
          
        FESystemElem::getPlateDKBatozQty
          (this,FESystemElem::PLATE_DK_BATOZ_HXX_HXX_FACTOR::num(), 
           &Hxx_Hxx_sens, 
           design_point_enum_ID,
           FESystemElem::ELEM_VOLUME::num(), true);
        FESystemElem::getPlateDKBatozQty
          (this, FESystemElem::PLATE_DK_BATOZ_HYY_HYY_FACTOR::num(), 
           &Hyy_Hyy_sens, 
           design_point_enum_ID,
           FESystemElem::ELEM_VOLUME::num(), true);
        FESystemElem::getPlateDKBatozQty
          (this,FESystemElem::PLATE_DK_BATOZ_HXX_HYY_FACTOR::num(), 
           &Hxx_Hyy_sens,
           design_point_enum_ID,
           FESystemElem::ELEM_VOLUME::num(), true);
        FESystemElem::getPlateDKBatozQty
          (this,FESystemElem::PLATE_DK_BATOZ_HXY_HXY_FACTOR::num(), 
           &Hxy_Hxy_sens,
           design_point_enum_ID,
           FESystemElem::ELEM_VOLUME::num(), true);
        FESystemElem::getPlateDKBatozQty
          (this,FESystemElem::PLATE_DK_BATOZ_HYX_HXY_FACTOR::num(), 
           &Hyx_Hxy_sens,
           design_point_enum_ID,
           FESystemElem::ELEM_VOLUME::num(), true);
        FESystemElem::getPlateDKBatozQty
          (this,FESystemElem::PLATE_DK_BATOZ_HYX_HYX_FACTOR::num(), 
           &Hyx_Hyx_sens,
           design_point_enum_ID,
           FESystemElem::ELEM_VOLUME::num(), true);

          N_N_sens_mem_factor = mem_material_mat;
          N_N_sens_plate_factor = plate_material_mat;

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
            
            FESystemElem::getPlateDKBatozQty(this,
                                             FESystemElem::PLATE_DK_BATOZ_HXX_HXX_FACTOR::num(), 
                                             &Hxx_Hxx, design_point_enum_ID,
                                             FESystemElem::ELEM_VOLUME::num(), false);
            FESystemElem::getPlateDKBatozQty(this,
                                             FESystemElem::PLATE_DK_BATOZ_HYY_HYY_FACTOR::num(), 
                                             &Hyy_Hyy, design_point_enum_ID,
                                             FESystemElem::ELEM_VOLUME::num(), false);
            FESystemElem::getPlateDKBatozQty(this,
                                             FESystemElem::PLATE_DK_BATOZ_HXX_HYY_FACTOR::num(), 
                                             &Hxx_Hyy, design_point_enum_ID,
                                             FESystemElem::ELEM_VOLUME::num(), false);
            FESystemElem::getPlateDKBatozQty(this,
                                             FESystemElem::PLATE_DK_BATOZ_HXY_HXY_FACTOR::num(), 
                                             &Hxy_Hxy, design_point_enum_ID,
                                             FESystemElem::ELEM_VOLUME::num(), false);
            FESystemElem::getPlateDKBatozQty(this,
                                             FESystemElem::PLATE_DK_BATOZ_HYX_HXY_FACTOR::num(), 
                                             &Hyx_Hxy, design_point_enum_ID,
                                             FESystemElem::ELEM_VOLUME::num(), false);
            FESystemElem::getPlateDKBatozQty(this,
                                             FESystemElem::PLATE_DK_BATOZ_HYX_HYX_FACTOR::num(), 
                                             &Hyx_Hyx, design_point_enum_ID,
                                             FESystemElem::ELEM_VOLUME::num(), false);
            N_N_mem_factor = mem_material_mat_temp_sens;
            N_N_plate_factor = plate_material_mat_temp_sens;
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
    
    FESystemElem::getPlateDKBatozQty(this,FESystemElem::PLATE_DK_BATOZ_HXX_HXX_FACTOR::num(), 
                                     &Hxx_Hxx, design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), false);
    FESystemElem::getPlateDKBatozQty(this,FESystemElem::PLATE_DK_BATOZ_HYY_HYY_FACTOR::num(), 
                                     &Hyy_Hyy, design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), false);
    FESystemElem::getPlateDKBatozQty(this,FESystemElem::PLATE_DK_BATOZ_HXX_HYY_FACTOR::num(), 
                                     &Hxx_Hyy, design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), false);
    FESystemElem::getPlateDKBatozQty(this,FESystemElem::PLATE_DK_BATOZ_HXY_HXY_FACTOR::num(), 
                                     &Hxy_Hxy, design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), false);
    FESystemElem::getPlateDKBatozQty(this,FESystemElem::PLATE_DK_BATOZ_HYX_HXY_FACTOR::num(), 
                                     &Hyx_Hxy, design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), false);
    FESystemElem::getPlateDKBatozQty(this,FESystemElem::PLATE_DK_BATOZ_HYX_HYX_FACTOR::num(), 
                                     &Hyx_Hyx, design_point_enum_ID,
                                     FESystemElem::ELEM_VOLUME::num(), false);
    N_N_mem_factor = mem_material_mat;
    N_N_plate_factor = plate_material_mat;
    }
  
	
	// calculate the stiffness matrix. Add the membrane stiffness effects for u, v dofs
  static double factor1, factor2, factor3, factor1_sens, factor2_sens, factor3_sens;
  factor1 = N_N_mem_factor(0,0);
  factor2 = N_N_mem_factor(0,1);
  factor3 = N_N_mem_factor(2,2);
  
  factor1_sens = N_N_sens_mem_factor(0,0);
  factor2_sens = N_N_sens_mem_factor(0,1);
  factor3_sens = N_N_sens_mem_factor(2,2);
  
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
      
    // now add the stiffness effects of the plate bending for the w, thetax, thetay dofs
    for (unsigned int i=0; i<n_nodes*3; i++)
      for (unsigned int j=0; j<n_nodes*3; j++)
        {
        (*matrix)(2*n_nodes+i,2*n_nodes+j) += 
        N_N_plate_factor(0,0)* (Hxx_Hxx(i,j) + Hyy_Hyy(i,j)) +
        N_N_plate_factor(0,1)* (Hxx_Hyy(i,j) + Hxx_Hyy(j,i)) +
        N_N_plate_factor(2,2)* (Hxy_Hxy(i,j) + Hyx_Hxy(i,j) + Hyx_Hxy(j,i) + Hyx_Hyx(i,j)) + 
        N_N_sens_plate_factor(0,0)* (Hxx_Hxx_sens(i,j) + Hyy_Hyy_sens(i,j)) +
        N_N_sens_plate_factor(0,1)* (Hxx_Hyy_sens(i,j) + Hxx_Hyy_sens(j,i)) +
        N_N_sens_plate_factor(2,2)* (Hxy_Hxy_sens(i,j) + Hyx_Hxy_sens(i,j) + 
                                     Hyx_Hxy_sens(j,i) + Hyx_Hyx_sens(i,j));
        }
        
    // finally, place a unity at the diagonal for the theta_z rows, to remove singularity 
    // of the element
    // this will be changed in future, since the value of this fictitious stiffness is dependent on the 
    // value of the other diagonal terms 
    for (unsigned int i=0; i<n_nodes; i++)
      (*matrix)(5*n_nodes+i,5*n_nodes+i) = 1.0;
}




void
FESystemElem::PlateDKBatoz::calculate_F_T(DenseVector<double>* vector, 
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
FESystemElem::PlateDKBatoz::getElementPostProcessQty
(std::vector<unsigned int> load_cases, 
 std::vector<DesignData::DesignParameter*> dv_vector)
{
  (void) load_cases;
  (void) dv_vector;
//  std::auto_ptr<ElemPostProcessQty> return_qty(new ElemPostProcessQty(this->elem_ID));
//  
//  TensorValue<double> tensor;
//  
//  // create an fe that will be used for calculationi of the element strains and
//  // stresses. The vector fe_sens is for shape sensitivity
////  std::auto_ptr<FEBase> fe(new FE<2,LAGRANGE>(FEType())),
////    fe_sens (new FE<2,LAGRANGE>(FEType()));
//    
//  
//  unsigned int n_nodes = this->base_elem->n_nodes();
//  
//  static DenseMatrix<double> B_strain(3, 6*n_nodes) , B_strain_sens(3, 6*n_nodes);
//  B_strain.zero(); B_strain_sens.zero();
//  
//  // create rest of the necessary data structures
//  static DenseVector<double>  strain(3), strain_sens(3),
//    stress(3), stress_sens(3), strain_plus_thermal(3), 
//    strain_plus_thermal_sens(3), strain_scratch(3);
//  
//  strain.zero(); strain_sens.zero();
//  stress.zero(); stress_sens.zero(); strain_plus_thermal.zero(); 
//  strain_plus_thermal_sens.zero(); strain_scratch.zero();
//  
//  static DenseMatrix<double> material_factor(3,3), material_factor_sens(3,3);
//  material_factor.zero(); 
//  material_factor_sens.zero();
//  
//  double temp_strain_material_factor = 0.0, temp_strain_material_factor_sens = 0.0,
//    ref_temp = 0.0, thickness = 0.0, thickness_sens = 0.0;
//    
//  // point at which this element will be initialized
//  // for now, this is only the centroid of the element
//  std::vector<Point> point_vec;
//  point_vec.push_back(Point(0.0,0.0,0.0));
//  
//  // then, reinit it at the base elem, at the element centroid, where the 
//  // element strains and stresses will be evaluated
//  this->initialize_element(&point_vec);
//  
//  // get the necessary data out of this element 
//  const std::vector< std::vector<Real> >& dphi_dx = this->fe_base->get_dphidx();
//  const std::vector< std::vector<Real> >& dphi_dy = this->fe_base->get_dphidy();
//  const std::vector<std::vector<Real> >& phi = this->fe_base->get_phi();
//  
//  // calculate B_strain matrix. Since the calculation is being performed at a single point, 
//  // the vector of shape factor derivatives should have size = 1
//  assert (dphi_dx[0].size() == 1);
//  for (unsigned int i=0; i < n_nodes; i++)
//    {
//    B_strain(0,i) = dphi_dx[i][0];  // epsilon_xx
//    B_strain(0,4*n_nodes+i) = 0.5 * thickness * dphi_dx[i][0]; // bending contribution
//    B_strain(1,n_nodes+i) = dphi_dy[i][0];  // epsilon_yy
//    B_strain(1,3*n_nodes+i) = - 0.5 * thickness * dphi_dy[i][0]; // bending contribution
//    B_strain(2,i) = dphi_dy[i][0];  // epsilon  xy 
//    B_strain(2,3*n_nodes+i) = - 0.5 * thickness * dphi_dx[i][0]; // bending contribution
//    B_strain(2,n_nodes+i) = dphi_dx[i][0];  // epsilon xy
//    B_strain(2,4*n_nodes+i) = 0.5 * thickness * dphi_dy[i][0]; // bending contribution
//    }
//  
//  
//  // iterate over all the load cases, to calculate the strains and stresses
//  std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
//  std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
//	
//  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
//  dv_begin = dv_vector.begin();
//  dv_end = dv_vector.end();
//  
//  
//  static DenseMatrix<double> transform_mat, transform_mat_sens;
//  transform_mat.zero(); transform_mat_sens.zero();
//  
//  this->getStructuralT_matrix(&transform_mat, false);
//  
//  
//  // also, create dof_vectors
//  static DenseVector<double>  dof_values(6*n_nodes), dof_value_sens(6*n_nodes),
//    local_dof(6*n_nodes), local_dof_sens(6*n_nodes), scratch_vec(6*n_nodes),
//    nodal_temp(n_nodes),nodal_temp_sens(n_nodes),
//    dphi_dx_sens(n_nodes), dphi_dy_sens(n_nodes);
//  
//  dof_values.zero(); dof_value_sens.zero();
//  local_dof.zero(); local_dof_sens.zero(); 
//  nodal_temp.zero(); nodal_temp_sens.zero();
//  dphi_dx_sens.zero(); dphi_dy_sens.zero(); 
//  scratch_vec.zero();
//
//  
//  // initialize the material factors
//  this->elem_property_card->getPropertyValueFromMaterialCard(ALPHA_EXPANSION::num(),
//                                                             temp_strain_material_factor);
//  this->elem_property_card->getPropertyValueFromMaterialCard(TEMP_REF::num(),
//                                                             ref_temp);
//  this->elem_property_card->getPropertyValueFromMaterialCard(THICKNESS_2D_ELEM::num(),
//                                                             thickness);
//  
//  this->elem_property_card->getFactor(material_factor, STRESS_STRAIN_FACTOR::num());
//  
//  double temp = 0.0, temp_sens = 0.0;
//  
//  // iterate over each load case, and ask solver to solve for it
//  for (; load_case_it != load_case_end; load_case_it++)
//    {
//    
//    // get the DOF vector for this elem
//    this->getElemDofValuesForLoadCase(dof_values, *load_case_it);
//    
//    // transform these dofs to the local coordinate system
//    local_dof.zero();
//    transform_mat.right_multiply_vector(dof_values, local_dof);
//    
//    // get the loads for this element
//    this->getNodalTemperatureVector(nodal_temp, *load_case_it);
//    
//    // calculate the factor and the strain and stress
//    B_strain.right_multiply_vector(local_dof, strain);
//    
//    // add the thermal strain component and then calculate the stress
//    for (unsigned int i=0; i<phi.size(); i++)
//      temp += phi[i][0] * nodal_temp(i);
//    
//    strain_plus_thermal.zero();
//    strain_plus_thermal(0) = strain(0) - temp_strain_material_factor * (temp - ref_temp);
//    strain_plus_thermal(1) = strain(1) - temp_strain_material_factor * (temp - ref_temp);
//    strain_plus_thermal(2) = strain(2);
//    
//    // now multiply to calculate the stresses
//    material_factor.right_multiply_vector(strain_plus_thermal, stress);
//    
//    
//    // add the two tensors
//    tensor.zero();
//    tensor(0,0) = strain(0);  // xx
//    tensor(1,1) = strain(1);  // yy
//    tensor(0,1) = strain(2);  // xy
//    tensor(1,0) = strain(2);  // yx
//    return_qty->addStrainTensor(tensor, *load_case_it);
//    
//    tensor.zero();
//    tensor(0,0) = stress(0);  // xx
//    tensor(1,1) = stress(1);  // yy
//    tensor(0,1) = stress(2);  // xy
//    tensor(1,0) = stress(2);  // yx
//    return_qty->addStressTensor(tensor, *load_case_it);
//    
//    tensor.zero();
//    tensor(0,0) = strain_plus_thermal(0);  // xx
//    tensor(1,1) = strain_plus_thermal(1);  // yy
//    tensor(0,1) = strain_plus_thermal(2);  // xy
//    tensor(1,0) = strain_plus_thermal(2);  // yx
//    return_qty->addMechanicalStrainTensor(tensor, *load_case_it);
//    
//    unsigned int dv_ID = 0;
//
//    // now calculate the sensitivities
//    dv_it = dv_begin;
//    for (; dv_it != dv_end; dv_it++)
//      {
//      dv_ID = (*dv_it)->getID();
//      this->clearSensitivityInitialization();
//      
//      // get the dof value sensitivity for this case
//      this->getElemDofValueSensitivityForLoadCaseAndDV(dof_value_sens, *load_case_it, dv_ID);
//      
//      // get the load sensitivity vectors
//      this->getNodalTemperatureVector(nodal_temp_sens, *load_case_it, true, dv_ID);
//      
//      material_factor_sens.zero();
//      temp_strain_material_factor_sens = 0.0;
//      
//      switch ((*dv_it)->getParameterTypeEnumID())
//        {
//        case PROPERTY_PARAMETER_ENUM_ID:
//          {
//            this->reinitForPropertySensitivity(dv_ID);
//            
//            // for a property DV, the strain operator sensitivity will be zero
//            B_strain_sens.zero();
//            
//            // if the property ID is the same as this elements property, then set the 
//            // value of the property sensitivity. Else the value be zero
//            if ( this->elem_property_card->checkGlobalParameterDependence(dv_ID))
//              {
//              this->elem_property_card->getFactorSensitivityForGlobalParameter
//              (material_factor_sens, STRESS_STRAIN_FACTOR::num(), dv_ID);
//              this->elem_property_card->
//                getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//                (ALPHA_EXPANSION::num(), dv_ID, temp_strain_material_factor_sens);
//              this->elem_property_card->
//                getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//                (THICKNESS_2D_ELEM::num(), dv_ID, thickness_sens);
//              }
//            else 
//              {
//              material_factor_sens.zero();
//              
//              temp_strain_material_factor_sens = 0.0;
//              }
//            
//            // also, calculate the sensitivity of the dof vector. The transformation matrix sensitivity 
//            // will be zero for this case, since shape sensitivity will be zero
//            transform_mat.right_multiply_vector(dof_value_sens, local_dof_sens);
//          }
//          break;
//          
//        case SHAPE_PARAMETER_ENUM_ID:
//          {
//            Elem* pert_elem = NULL;
//            pert_elem = this->analysis_discipline.getPerturbedElemForShapeParameter(this->elem_ID,
//                                                                                    *dv_it);
//            this->reinitForShapeSensitivity(dv_ID, pert_elem, 
//                                            (*dv_it)->getPerturbationStepSize());
//            
//            // set the factor_sens to zero
//            material_factor_sens.zero();
//            temp_strain_material_factor_sens = 0.0;
//            
//            // get the transformation matrix sensitivity and calculate the local_dof_sens
//            this->getStructuralT_matrix(&transform_mat_sens, true);
//            
//            local_dof_sens.zero(); scratch_vec.zero();
//            transform_mat_sens.right_multiply_vector(dof_values, local_dof_sens);
//            
//            scratch_vec.zero();
//            transform_mat.right_multiply_vector(dof_value_sens, scratch_vec);
//            
//            local_dof_sens.add(1.0, scratch_vec);
//            
//            // now init the fe_sens to the perturbed DV and calculate the shape sensitivity of the
//            // B matrix
//            // then, reinit it at the base elem, at the element centroid, where the 
//            // element strains and stresses will be evaluated
//            this->attachNewElem(this->perturbed_elem);
//            this->initialize_element(&point_vec);
//
//            // get the necessary data out of this element 
//            const std::vector< std::vector<Real> >& dphi_dx_perturbed = 
//              this->fe_base->get_dphidx();
//            const std::vector< std::vector<Real> >& dphi_dy_perturbed = 
//              this->fe_base->get_dphidy();
//            
//            dphi_dx_sens.zero();
//            dphi_dy_sens.zero();
//            
//            for (unsigned int i=0; i<n_nodes; i++)
//              {
//              dphi_dx_sens(i) = dphi_dx_perturbed[i][0] - dphi_dx[i][0];
//              dphi_dy_sens(i) = dphi_dy_perturbed[i][0] - dphi_dy[i][0];
//              }
//            
//            dphi_dx_sens.scale(1.0/this->perturbation);
//            dphi_dy_sens.scale(1.0/this->perturbation);
//            
//            // calculate B_strain matrix. Since the calculation is being performed at a single point, 
//            // the vector of shape factor derivatives should have size = 1
//            
//            for (unsigned int i=0; i < n_nodes; i++)
//              {
//              B_strain_sens(0,i) = dphi_dx_sens(i);  // epsilon_xx
//              B_strain_sens(0,4*n_nodes+i) = 0.5 * (thickness * dphi_dx_sens(i) + 
//                                                    thickness_sens * dphi_dx[i][0]); // bending contribution
//              B_strain_sens(1,n_nodes+i) = dphi_dy_sens(i);  // epsilon_yy
//              B_strain_sens(1,3*n_nodes+i) = - 0.5 * (thickness * dphi_dy_sens(i) + 
//                                                      thickness_sens * dphi_dy[i][0]); // bending contribution
//              B_strain_sens(2,i) = dphi_dy_sens(i);  // epsilon  xy 
//              B_strain_sens(2,3*n_nodes+i) = -0.5 * (thickness * dphi_dx_sens(i) + 
//                                                     thickness_sens * dphi_dx[i][0]); // bending contribution
//              B_strain_sens(2,n_nodes+i) = dphi_dx_sens(i);  // epsilon xy
//              B_strain_sens(2,4*n_nodes+i) = 0.5 * (thickness * dphi_dy_sens(i) + 
//                                                    thickness_sens * dphi_dy[i][0]); // bending contribution
//              }
//            
//          }
//          break;
//          
//        default:
//          abort();
//          break;
//        }
//      
//      // calculate the strain and stress
//      // calculate the factor and the strain and stress
//      strain_scratch.zero();
//      B_strain.right_multiply_vector(local_dof_sens, strain_sens);
//      B_strain_sens.right_multiply_vector(local_dof, strain_scratch);
//      strain_sens.add(1.0, strain_scratch);
//      
//      // add the thermal strain component and then calculate the stress
//      temp_sens = 0.0;
//      for (unsigned int i=0; i<phi.size(); i++)
//        temp_sens += phi[i][0] * nodal_temp_sens(i);
//      
//      // add the temperture sensitivity effects to the strain sensitivity
//      strain_plus_thermal_sens.zero();
//      strain_plus_thermal_sens(0) = strain_sens(0) - temp_strain_material_factor * temp_sens -
//        temp_strain_material_factor_sens * (temp - ref_temp);
//      strain_plus_thermal_sens(1) = strain_sens(1) - temp_strain_material_factor * temp_sens -
//        temp_strain_material_factor_sens * (temp - ref_temp);
//      strain_plus_thermal_sens(2) = strain_sens(2);
//      
//      // now multiply to calculate the stresses
//      strain_scratch.zero();
//      material_factor_sens.right_multiply_vector(strain_plus_thermal, stress_sens);
//      material_factor.right_multiply_vector(strain_plus_thermal_sens, strain_scratch);
//      stress_sens.add(1.0, strain_scratch);
//      
//      // add the two tensors
//      tensor.zero();
//      tensor(0,0) = strain_sens(0);  // xx
//      tensor(1,1) = strain_sens(1);  // yy
//      tensor(0,1) = strain_sens(2);  // xy
//      tensor(1,0) = strain_sens(2);  // yx
//      return_qty->addStrainTensor(tensor, *load_case_it,  dv_ID);
//      
//      tensor.zero();
//      tensor(0,0) = stress_sens(0);  // xx
//      tensor(1,1) = stress_sens(1);  // yy
//      tensor(0,1) = stress_sens(2);  // xy
//      tensor(1,0) = stress_sens(2);  // yx
//      return_qty->addStressTensor(tensor, *load_case_it, dv_ID);
//      
//      tensor.zero();
//      tensor(0,0) = strain_plus_thermal_sens(0);  // xx
//      tensor(1,1) = strain_plus_thermal_sens(1);  // yy
//      tensor(0,1) = strain_plus_thermal_sens(2);  // xy
//      tensor(1,0) = strain_plus_thermal_sens(2);  // yx
//      return_qty->addMechanicalStrainTensor(tensor, *load_case_it,  dv_ID);
//      
//      }
//    }  
//  
//  return return_qty;
}




