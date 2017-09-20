// $Id: linear_spring.C,v 1.21.6.1 2007-03-14 22:05:03 manav Exp $

// C++ includes


// FESystem includes
#include "StructuralElems/linear_spring.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic1DElemDataCard.h"

// libMesh includes
#include "geom/edge_edge2.h"


FESystemElem::LinearSpring::LinearSpring(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::StructuralElem(1, FESystemElem::STRUCTURAL_LINEAR_SPRING_EDGE2::num(), discipline)
{
  
}



FESystemElem::LinearSpring::~LinearSpring()
{
	
}


void
FESystemElem::LinearSpring::calculate_K_G(DenseMatrix<double>* matrix, 
                                          const unsigned int design_point_enum_ID,
                                          bool sensitivity_calculation)
{
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;

  abort();
}



void
FESystemElem::LinearSpring::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
}





void
FESystemElem::LinearSpring::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
}




void
FESystemElem::LinearSpring::calculate_M(DenseMatrix<double>* matrix, 
                                        const unsigned int design_point_enum_ID,
                                        bool sensitivity_calculation)
{
  // params not used here
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;
  
}



void
FESystemElem::LinearSpring::calculate_K(DenseMatrix<double>* matrix, 
                                        const unsigned int design_point_enum_ID,
                                        bool sensitivity_calculation)
{
  // param not used here
  (void) design_point_enum_ID;
		
  static double factor, factor_sens, factor_temp_sens; 
  factor = 0.0; factor_sens = 0.0; factor_temp_sens = 0.0; 
  
  this->elem_property_card->getFactor(factor, SPRING_STIFFNESS_FACTOR::num());

  if (sensitivity_calculation)
    {
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (factor_sens, SPRING_STIFFNESS_FACTOR::num(), this->sensitivity_parameter_ID);
        
      default:
        if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
          this->elem_property_card->getFactorSensitivityForLocalParameter
            (factor_temp_sens, SPRING_STIFFNESS_FACTOR::num(), 
             Property::TEMPERATURE::num());
        break;
      }
    }

  // set the factor values in the matrix
  (*matrix)(0,0) = 1.0; // u1,u1
  (*matrix)(0,1) = -1.0; // u1,u2
  (*matrix)(1,0) = -1.0; // u2,u1
  (*matrix)(1,1) = 1.0; // u2,u2
	
	
  // if the property is temperature dependent, then the temperature sensitivity 
  // will be needed
  static DenseVector<double> nodal_temp_sens(2);
  nodal_temp_sens.zero();
  
  if (sensitivity_calculation == true)
    {
    static double avg_temp_sens, avg_factor_temp_sens;
    avg_temp_sens = 0.0; avg_factor_temp_sens = 0.0;
    
    if (this->analysis_discipline.checkPropertyDependenceOnTemperature())
      {
      this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
      for (unsigned int i=0; i < 2; i++)
        avg_temp_sens += nodal_temp_sens(i);
      avg_temp_sens *= 0.5;
      
      avg_factor_temp_sens = factor_temp_sens * avg_temp_sens;
      }
    
    switch (this->sensitivity_parameter)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          matrix->scale(factor_sens + avg_factor_temp_sens);
        }
        break;
				
      case SHAPE_PARAMETER_ENUM_ID:
        {
          // nothing changes for change of shape of the element in its local axis.
          // hence, just zero the matrix. But for temperature dependent property, the 
	  // sensitivity is to be calculated
          matrix->scale(avg_factor_temp_sens);
        }
        break;
				
      default:
        abort();
        break;
      }
    }
  else 
    {
    matrix->scale(factor);
    }
}




void
FESystemElem::LinearSpring::calculate_F_T(DenseVector<double>* vector, 
                                          const unsigned int design_point_enum_ID,
                                          bool sensitivity_calculation)
{
  // params not used here
  (void) vector;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;
}




//void
//FESystemElem::LinearSpring::calculate_F_Pressure(DenseVector<double>* vector, 
//                                                 const unsigned int design_point_enum_ID,
//                                                 bool sensitivity_calculation)
//{
//  vector->zero();
//}




std::auto_ptr<ElemPostProcessQty> 
FESystemElem::LinearSpring::getElementPostProcessQty(std::vector<unsigned int> load_cases, 
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
//    return_qty->addMechanicalStrainTensor(tensor, *load_case_it);
//    
//    dv_it = dv_begin;
//    for (; dv_it != dv_end; dv_it++)
//      {
//      return_qty->addStrainTensor(tensor, *load_case_it, (*dv_it)->getID());
//      return_qty->addStressTensor(tensor, *load_case_it, (*dv_it)->getID());
//      return_qty->addMechanicalStrainTensor(tensor, *load_case_it, (*dv_it)->getID());
//      }
//    }
//  
//  return return_qty;  
}


// void  LinearSpring::calculateStrainOperator(std::vector<DenseMatrix<double> >* operator_vec, 
// 					    const std::vector<Point>* points)
// {
//   // this element has no strain. Hence, this does ont implement anything. For now,
//   // this method will abort the program. But will be changed later

//   abort();
// }
