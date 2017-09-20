// $Id: ForcedConvection1D.C,v 1.3 2007-01-15 19:00:04 manav Exp $

// FESystem includes
#include "ThermalElems/ForcedConvection1D.h"
#include "ThermalElems/thermal_elem.h"
#include "Discipline/ThermalAnalysis.h"
#include "Properties/ForcedConvection1DElemDataCard.h"


// libMesh includes
#include "numerics/dense_matrix.h"



FESystemElem::ForcedConvection_1D::ForcedConvection_1D(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::ThermalElem(2, FESystemElem::THERMAL_FORCED_CONVECTION_EDGE2::num(), discipline),
upwind_factor(0.0)
{
  
}




FESystemElem::ForcedConvection_1D::~ForcedConvection_1D()
{
  
}




void
FESystemElem::ForcedConvection_1D::calculate_C(DenseMatrix<double>* matrix, 
                                               const unsigned int design_point_enum_ID,
                                               bool sensitivity_calculation)
{
  static double fluid_factor, fluid_factor_sens, wall_factor, wall_factor_sens, length; 
  fluid_factor = 0.0, fluid_factor_sens = 0.0; wall_factor = 0.0, wall_factor_sens = 0.0; length = 0.0;

  static Elem* elem;
  static Point tmp_vector1; 
  elem = this->geometric_elems_for_DV_map[design_point_enum_ID];
  
  {  // calculate the length of the element
    const Point& point0 = elem->point(0);
    const Point& point1 = elem->point(1);
    tmp_vector1.assign(point1 - point0);
    length = tmp_vector1.size();
  }    

  
  this->elem_property_card->getFactor(fluid_factor, FLUID_CAPACITANCE_FACTOR::num());
  this->elem_property_card->getFactor(wall_factor, WALL_CAPACITANCE_FACTOR::num());

  if (sensitivity_calculation && 
      this->sensitivity_parameter == DesignData::PROPERTY_PARAMETER::num())
    {
    if (!this->elem_property_card->checkElemAndMaterialCardGlobalParameterDependence
        (this->sensitivity_parameter_ID))
      return;
    
    this->elem_property_card->getFactorSensitivityForGlobalParameter
      (fluid_factor_sens, FLUID_CAPACITANCE_FACTOR::num(), this->sensitivity_parameter_ID);
    this->elem_property_card->getFactorSensitivityForGlobalParameter
      (wall_factor_sens, WALL_CAPACITANCE_FACTOR::num(), this->sensitivity_parameter_ID);
    }
  
  static double val1, val2, val3, val4;
  val1=0.0; val2=0.0; val3=0.0; val4=0.0;
    
  if(sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          val1 = fluid_factor_sens * length;
          val2 = this->upwind_factor * fluid_factor_sens * length * .25;
          val3 = wall_factor_sens * length;
          val4 = this->upwind_factor * wall_factor_sens * length * .25;
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          val1 = fluid_factor;
          val2 = this->upwind_factor * fluid_factor * .25;
          val3 = wall_factor;
          val4 = this->upwind_factor * wall_factor * .25;
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else
    {
    val1 = fluid_factor * length;
    val2 = this->upwind_factor * fluid_factor * length * .25;
    val3 = wall_factor * length;
    val4 = this->upwind_factor * wall_factor * length * .25;
    }
  
  // fluid factors
  (*matrix)(0,0) = val1 / 3.0 + (-val2);
  (*matrix)(0,1) = val1 / 6.0 + (-val2);
  (*matrix)(1,0) = val1 / 6.0 + val2;
  (*matrix)(1,1) = val1 / 3.0 + val2;
  
  // wall terms will be different to take care of the numbering
  (*matrix)(2,2) = val3 / 3.0 + val4;
  (*matrix)(2,3) = val3 / 6.0 + val4;
  (*matrix)(3,2) = val3 / 6.0 + (-val4);
  (*matrix)(3,3) = val3 / 3.0 + (-val4);
}







void
FESystemElem::ForcedConvection_1D::calculate_K_c(DenseMatrix<double>* matrix, 
                                         const unsigned int design_point_enum_ID,
                                         bool sensitivity_calculation)
{
  static double fluid_conduction_factor, fluid_conduction_factor_sens,
  fluid_convective_transport_factor, fluid_convective_transport_factor_sens,
  convection_factor, convection_factor_sens,
  wall_conduction_factor, wall_conduction_factor_sens, length, length_sens;

  fluid_conduction_factor=0.0; fluid_conduction_factor_sens=0.0; 
  fluid_convective_transport_factor=0.0; fluid_convective_transport_factor_sens=0.0; 
  convection_factor=0.0; convection_factor_sens=0.0; 
  wall_conduction_factor=0.0; wall_conduction_factor_sens=0.0; length=0.0; length_sens=0.0;
  
  static Elem* elem;
  static Point tmp_vector1; 
  elem = this->geometric_elems_for_DV_map[design_point_enum_ID];
  
  {  // calculate the length of the element
    const Point& point0 = elem->point(0);
    const Point& point1 = elem->point(1);
    tmp_vector1.assign(point1 - point0);
    length = tmp_vector1.size();
  }    
  
  
  this->elem_property_card->getFactor(fluid_conduction_factor, FLUID_CONDUCTIVITY_FACTOR::num());
  this->elem_property_card->getFactor(fluid_convective_transport_factor, FLUID_CONVECTIVE_FACTOR::num());
  this->elem_property_card->getFactor(convection_factor, CONVECTIVE_EXCANGE_FACTOR::num());
  this->elem_property_card->getFactor(wall_conduction_factor, WALL_CONDUCTIVITY_FACTOR::num());
  
  if (sensitivity_calculation && 
      this->sensitivity_parameter == DesignData::PROPERTY_PARAMETER::num())
    {
    if (!this->elem_property_card->checkElemAndMaterialCardGlobalParameterDependence
        (this->sensitivity_parameter_ID))
      return;
    
    this->elem_property_card->getFactorSensitivityForGlobalParameter(fluid_conduction_factor_sens,
                                                                     FLUID_CONDUCTIVITY_FACTOR::num(),
                                                                     this->sensitivity_parameter_ID);
    this->elem_property_card->getFactorSensitivityForGlobalParameter(fluid_convective_transport_factor_sens,
                                                                     FLUID_CONVECTIVE_FACTOR::num(),
                                                                     this->sensitivity_parameter_ID);
    this->elem_property_card->getFactorSensitivityForGlobalParameter(convection_factor_sens,
                                                                     CONVECTIVE_EXCANGE_FACTOR::num(),
                                                                     this->sensitivity_parameter_ID);
    this->elem_property_card->getFactorSensitivityForGlobalParameter(wall_conduction_factor_sens,
                                                                     WALL_CONDUCTIVITY_FACTOR::num(),
                                                                     this->sensitivity_parameter_ID);
    }
  
  static double kf, kw, kh, kh_alpha, kv, kv_alpha;
  kf=0.0;  kw=0.0;  kh=0.0;  kh_alpha=0.0;  kv=0.0;  kv_alpha=0.0;
  
  if(sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          kf = fluid_conduction_factor_sens / length;
          kw = wall_conduction_factor_sens / length;
          kh = convection_factor_sens * length;
          kh_alpha = this->upwind_factor * convection_factor_sens * length * 0.25;
          kv = fluid_convective_transport_factor_sens * 0.5;
          kv_alpha = this->upwind_factor * fluid_convective_transport_factor_sens * 0.5;
        }
        break;
        
      case SHAPE_PARAMETER_ENUM_ID:
        {
          elem = this->geometric_elems_for_DV_map[BASE_PLUS_DELTA_ELEM::num()];
          {  // calculate the length of the element
            const Point& point0 = elem->point(0);
            const Point& point1 = elem->point(1);
            tmp_vector1.assign(point1 - point0);
            length_sens = (tmp_vector1.size() -length)/this->perturbation;
          }
          
          kf = fluid_conduction_factor * (-1.0/ length / length) * length_sens;
          kw = wall_conduction_factor * (-1.0/ length / length) * length_sens;
          kh = convection_factor;
          kh_alpha = this->upwind_factor * convection_factor * 0.25;
          kv = 0.0;
          kv_alpha = 0.0;
        }
        break;
        
      default:
        abort();
        break;
      }
    }
  else
    {
    kf = fluid_conduction_factor / length;
    kw = wall_conduction_factor / length;
    kh = convection_factor * length;
    kh_alpha = this->upwind_factor * convection_factor * length * 0.25;
    kv = fluid_convective_transport_factor * 0.5;
    kv_alpha = this->upwind_factor * fluid_convective_transport_factor * 0.5;
    }
  
  // fluid factors
  (*matrix)(0,0) = (-kv) + (kh / 3.0) + kf + kv_alpha + (-kh_alpha);
  (*matrix)(0,1) = kv + (kh / 6.0) - kf + (-kv_alpha) + (-kh_alpha);
  (*matrix)(1,0) = (-kv) + (kh / 6.0) - kf + (-kv_alpha) + kh_alpha;
  (*matrix)(1,1) = kv + (kh / 3.0) + kf + kv_alpha + kh_alpha;
  
  (*matrix)(0,2) = (-kh / 6.0) + kh_alpha;
  (*matrix)(0,3) = (-kh / 3.0) + kh_alpha;
  (*matrix)(1,2) = (-kh / 3.0) - kh_alpha;
  (*matrix)(1,3) = (-kh / 6.0) - kh_alpha;
  
  // the terms here have been rearranged since the node numbering for the wall nodes on
  // the quad element are in reverse order than the formulation
  (*matrix)(2,0) = (-kh / 6.0) -kh_alpha;
  (*matrix)(2,1) = (-kh / 3.0) -kh_alpha;
  (*matrix)(3,0) = (-kh / 3.0) + kh_alpha;
  (*matrix)(3,1) = (-kh / 6.0) + kh_alpha;
  
  
  // wall terms will be different to take care of the numbering
  (*matrix)(2,2) = (kh / 3.0) + kw + kh_alpha;
  (*matrix)(2,3) = (kh / 6.0) - kw + kh_alpha;
  (*matrix)(3,2) = (kh / 6.0) - kw - kh_alpha;
  (*matrix)(3,3) = (kh / 3.0) + kw - kh_alpha;
}




void
FESystemElem::ForcedConvection_1D::calculate_K_c_Jac(DenseMatrix<double>* matrix,
                                             const unsigned int design_point_enum_ID)
{
  if (!this->analysis_discipline.checkPropertyDependenceOnTemperature())
    return;
  // this needs to be implemented
  matrix->zero();
}




void
FESystemElem::ForcedConvection_1D::calculate_K_h(DenseMatrix<double>* matrix, 
                                                 const unsigned int design_point_enum_ID,
                                                 bool sensitivity_calculation)
{

}





void
FESystemElem::ForcedConvection_1D::calculate_F_qsurf(DenseVector<double>* vector, 
                                             const unsigned int design_point_enum_ID,
                                             bool sensitivity_calculation)
{

}





void
FESystemElem::ForcedConvection_1D::calculate_F_Qvol(DenseVector<double>* vector, 
                                            const unsigned int design_point_enum_ID,
                                            bool sensitivity_calculation)
{

}





void
FESystemElem::ForcedConvection_1D::calculate_F_h(DenseVector<double>* vector, 
                                         const unsigned int design_point_enum_ID,
                                         bool sensitivity_calculation)
{

}





void
FESystemElem::ForcedConvection_1D::calculate_F_sigma(DenseVector<double>* vector, 
                                             const unsigned int design_point_enum_ID,
                                             bool sensitivity_calculation)
{

}





void
FESystemElem::ForcedConvection_1D::calculate_F_sigma_Jac(DenseMatrix<double>* matrix,
                                                 const unsigned int design_point_enum_ID)
{

}

