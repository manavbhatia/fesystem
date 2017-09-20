// $Id: Isotropic2DElemDataCard.h,v 1.10 2006-12-30 09:53:55 manav Exp $

#ifndef __fesystem_isotropic_2d_elem_data_card_h__
#define __fesystem_isotropic_2d_elem_data_card_h__

// FESystem includes
#include "Properties/IsotropicElemDataCard.h"
#include "Properties/MaterialPropertyNameEnums.h"


#ifndef ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_ID
#define ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_ID 4
#else
#error
#endif


#ifndef ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_NAME
#define ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_NAME "ISOTROPIC_2D_ELEM_DATA_CARD"
#else
#error
#endif


#ifndef THICKNESS_2D_ELEM_ENUM_ID
#define THICKNESS_2D_ELEM_ENUM_ID 2
#else
#error
#endif

#ifndef THICKNESS_2D_ELEM_ENUM_NAME
#define THICKNESS_2D_ELEM_ENUM_NAME  "THICKNESS"
#else
#error
#endif



DeclareEnumName(ISOTROPIC_2D_ELEM_DATA_CARD,
                ElemDataCardEnum,
                ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_ID,
                ISOTROPIC_2D_ELEM_DATA_CARD_ENUM_NAME);


DeclareEnumClass(Isotropic_2D_ElemDataEnum);

DeclareEnumName(THICKNESS_2D_ELEM, Isotropic_2D_ElemDataEnum, 
                THICKNESS_2D_ELEM_ENUM_ID , THICKNESS_2D_ELEM_ENUM_NAME);




class Isotropic2D_ElemDataCard: public IsotropicElemDataCard
{
public:
  Isotropic2D_ElemDataCard();
  
  virtual  ~Isotropic2D_ElemDataCard();
  
  /// @returns the property card kind enumeration ID
  virtual unsigned int getPropertyCardKindEnumID() const;
  
  /// @returns the property card kind enumeration name
  virtual const std::string getPropertyCardKindEnumName() const;
  
  /// returns the value of the factor
  virtual void getFactor(double& factor, const unsigned int factor_enum_ID) ;
  
  virtual inline void getFactorSensitivityForGlobalParameter(double& factor,
                                                      const unsigned int factor_enum_ID,
                                                      const unsigned int global_param_ID);
  
  virtual inline void getFactorSensitivityForLocalParameter(double& factor,
                                                     const unsigned int factor_enum_ID,
                                                     const unsigned int local_param_enum_ID);
  
  /// returns the value of the factor
  virtual inline void getFactor(DenseMatrix<double>& factor, const unsigned int factor_enum_ID) ;
  
  virtual inline void getFactorSensitivityForGlobalParameter(DenseMatrix<double>& factor,
                                                      const unsigned int factor_enum_ID,
                                                      const unsigned int global_param_ID);
  
  virtual inline void getFactorSensitivityForLocalParameter(DenseMatrix<double>& factor,
                                                     const unsigned int factor_enum_ID,
                                                     const unsigned int local_param_enum_ID);
  
    
  friend std::istream& operator>> (std::istream& input, Isotropic2D_ElemDataCard& plane_elem);
protected:
    
  /// @returns the property name from enumeration ID of the property
  virtual std::string getPropertyEnumName(const unsigned int enum_ID) const;
  
  /// @returns the property ID from enumeration name of the property
  virtual const unsigned int getPropertyEnumID(const std::string& enum_name) const; 
};




inline
unsigned int
Isotropic2D_ElemDataCard::getPropertyCardKindEnumID() const
{
  return ISOTROPIC_2D_ELEM_DATA_CARD::num();
}


inline
const std::string
Isotropic2D_ElemDataCard::getPropertyCardKindEnumName() const
{
  return ISOTROPIC_2D_ELEM_DATA_CARD::name();
}






inline 
void
Isotropic2D_ElemDataCard::getFactor(double& factor,
                                    const unsigned int factor_enum_ID)
{
  Property::PropertyCard<double>& property = this->getMaterialCard();
  
  factor = 0.0;

  static double val1, val2, val3, val4;
  val1 = 0.0; val2 = 0.0; val3 = 0.0; val4 = 0.0;
  
  switch (factor_enum_ID)
    {
    
    case MEMBRANE_MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        factor = val1*val2;
      }
      break;
      
    case PLATE_MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        factor = val1*pow(val2,3)/12.0;
      }
      break;

    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(ALPHA_EXPANSION::num(), val3);
        property.getPropertyValue(POISSONS_RATIO::num(), val4);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        factor = val1 * val2 * val3 / (1.0 - val4);
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val3);
        
        factor = val1 * val2 * val3;
      }
      break;
      
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        factor = val1 * val2;
      }
      break;
      

    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        factor = val1 * val2;
      }
      break;

      
    case RADIATION_EPSILON_FACTOR_1_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        
        factor = 1.0 / val1;
      }
      break;
      
    case RADIATION_EPSILON_FACTOR_2_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        
        factor = (1.0 - val1) / val1;
      }
      break;

    
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}





inline 
void
Isotropic2D_ElemDataCard::getFactor(DenseMatrix<double>& matrix,
                                    const unsigned int factor_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  if (matrix.m() != 3 || matrix.n() != 3)
    matrix.resize(3,3);
  
  matrix.zero();

  static double val1, val2, val4;
  val1 = 0.0; val2 = 0.0; val4 = 0.0;
  
  switch (factor_enum_ID)
    {
    case STIFFNESS_A_MATRIX_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val4);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        matrix(0,0) = val1 / (1.0 - val4 * val4); // E / (1-nu * nu)
        matrix(1,1) = val1 / (1.0 - val4 * val4); // E / (1-nu * nu)
        matrix(0,1) = val1 * val4 / (1.0 - val4 * val4); // E * nu / (1-nu * nu)
        matrix(1,0) = val1 * val4 / (1.0 - val4 * val4); // E * nu / (1-nu * nu)
        matrix(2,2) = val1 / 2.0 / (1.0 + val4); // E / (2 *(1+nu))
        matrix.scale(val2); // multiply by thickness
      };
      break;
      
    case STIFFNESS_B_MATRIX_FACTOR_ENUM_ID:
      {
        // nothing to be done here, since for isotropic materials, this is zero
      };
      break;

    case STIFFNESS_D_MATRIX_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val4);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        matrix(0,0) = val1 / (1.0 - val4 * val4); // E / (1-nu * nu)
        matrix(1,1) = val1 / (1.0 - val4 * val4); // E / (1-nu * nu)
        matrix(0,1) = val1 * val4 / (1.0 - val4 * val4); // E * nu / (1-nu * nu)
        matrix(1,0) = val1 * val4 / (1.0 - val4 * val4); // E * nu / (1-nu * nu)
        matrix(2,2) = val1 / 2.0 / (1.0 + val4); // E / (2 *(1+nu))
        matrix.scale(val2 * val2 * val2 / 12.0); // multiply by h^3 / 12
      };
      break;
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}




inline 
void
Isotropic2D_ElemDataCard::
getFactorSensitivityForGlobalParameter(double& factor, 
                                       const unsigned int factor_enum_ID,
                                       const unsigned int global_param_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  factor = 0.0;
  
  static double val1, val2, val3, val4;
  static double val1_sens, val2_sens, val3_sens, val4_sens;
  
  val1 = 0.0; val2 = 0.0; val3 = 0.0; val4 = 0.0;
  val1_sens = 0.0; val2_sens = 0.0; val3_sens = 0.0; val4_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    
    case MEMBRANE_MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(),
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(THICKNESS_2D_ELEM::num(),
                                                      global_param_ID,
                                                      val2_sens);
        factor = val1*val2_sens + val1_sens * val2;
      }
      break;
      
    case PLATE_MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(),
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(THICKNESS_2D_ELEM::num(),
                                                      global_param_ID,
                                                      val2_sens);
        factor = val1*(0.25 * pow(val2, 2) * val2_sens) + val1_sens * pow(val2,3) / 12.0;
      }
      break;

    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(ALPHA_EXPANSION::num(), val3);
        property.getPropertyValue(POISSONS_RATIO::num(), val4);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(),
                                                               global_param_ID,
                                                               val1_sens);
        property.getPropertyDerivativeForGlobalParameter(ALPHA_EXPANSION::num(),
                                                    global_param_ID, val3_sens);
        property.getPropertyDerivativeForGlobalParameter(POISSONS_RATIO::num(),
                                                    global_param_ID, val4_sens);
        this->getPropertyDerivativeForGlobalParameter(THICKNESS_2D_ELEM::num(),
                                                 global_param_ID, val2_sens);

        factor = (val1_sens * val2 * val3 + val1 * val2_sens * val3 + 
                  val1 * val2 * val3_sens) / (1.0 - val4) - 
          (val1 * val2 * val3) * val4_sens / (1.0 - val4) / (1.0 - val4) ;
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val3);

        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(),
                                                               global_param_ID, val1_sens);
        property.getPropertyDerivativeForGlobalParameter(SPECIFIC_HEAT::num(),
                                                               global_param_ID, val2_sens);
        this->getPropertyDerivativeForGlobalParameter(THICKNESS_2D_ELEM::num(),
                                                            global_param_ID, val3_sens);
        
        factor = val1_sens * val2 * val3 + val1 * val2_sens * val3 + 
          val1 * val2 * val3_sens;
      }
      break;
      
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter
          (THERMAL_CONDUCTIVITY::num(),
           global_param_ID, val1_sens);
        this->getPropertyDerivativeForGlobalParameter
          (THICKNESS_2D_ELEM::num(),
           global_param_ID, val2_sens);

        factor = val1_sens * val2 + val1 * val2_sens;
      }
      break;
      
      
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter
          (EMISSIVITY::num(), global_param_ID, val1_sens);
        this->getPropertyDerivativeForGlobalParameter
          (THICKNESS_2D_ELEM::num(), global_param_ID, val2_sens);
        
        factor = val1_sens * val2 + val1 * val2_sens;
      }
      break;


    case RADIATION_EPSILON_FACTOR_1_ENUM_ID:
      case RADIATION_EPSILON_FACTOR_2_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        property.getPropertyDerivativeForGlobalParameter
          (EMISSIVITY::num(), global_param_ID, val1_sens);
        
        factor = (-1.0 / val1 / val1) * val1_sens;
      }
      break;
      
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}


inline 
void
Isotropic2D_ElemDataCard::
getFactorSensitivityForGlobalParameter(DenseMatrix<double>& matrix, 
                                       const unsigned int factor_enum_ID,
                                       const unsigned int global_param_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  if (matrix.m() != 3 || matrix.n() != 3)
    matrix.resize(3,3);
  
  static DenseMatrix<double> base_qty(3,3);
  
  base_qty.zero();
  matrix.zero();
  
  static double val1, val2, val3;
  static double val1_sens, val2_sens, val3_sens;
  
  val1 = 0.0; val2 = 0.0; val3 = 0.0; 
  val1_sens = 0.0; val2_sens = 0.0; val3_sens = 0.0;

  switch (factor_enum_ID)
    {
    
    case STIFFNESS_A_MATRIX_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(),
                                                               global_param_ID,
                                                               val1_sens);
        property.getPropertyDerivativeForGlobalParameter(POISSONS_RATIO::num(),
                                                               global_param_ID, 
                                                               val3_sens);
        this->getPropertyDerivativeForGlobalParameter(THICKNESS_2D_ELEM::num(),
                                                            global_param_ID, 
                                                            val2_sens);

        
        matrix(0,0) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(1,1) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(0,1) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(1,0) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(2,2) = val1_sens / 2.0 / (1.0 + val3) + 
          val1 * val3_sens / 2.0 / (1.0 + val3); // E / (2 *(1+nu))

        matrix.scale(val2); // multiply by thickness
        this->getFactor(base_qty, factor_enum_ID);
        matrix.add(val2_sens, base_qty);
      };
      break;
      
    case STIFFNESS_B_MATRIX_FACTOR_ENUM_ID:
      {
        // nothing to be done here, since for isotropic materials, this is zero
      };
      break;
      
    case STIFFNESS_D_MATRIX_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(),
                                                               global_param_ID,
                                                               val1_sens);
        property.getPropertyDerivativeForGlobalParameter(POISSONS_RATIO::num(),
                                                               global_param_ID,
                                                               val3_sens);
        this->getPropertyDerivativeForGlobalParameter(THICKNESS_2D_ELEM::num(),
                                                            global_param_ID,
                                                            val2_sens);

        
        matrix(0,0) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(1,1) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(0,1) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(1,0) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(2,2) = val1_sens / 2.0 / (1.0 + val3) + 
          val1 * val3_sens / 2.0 / (1.0 + val3); // E / (2 *(1+nu))
        matrix.scale(val2 * val2 * val2 / 12.0); // multiply by h^3 / 12
        this->getFactor(base_qty, factor_enum_ID);
        matrix.add(val2 * val2 * val2_sens / 4.0, base_qty);
      };
      break;
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}




inline 
void
Isotropic2D_ElemDataCard::
getFactorSensitivityForLocalParameter(double& factor,
                                      const unsigned int factor_enum_ID,
                                      const unsigned int param_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  factor = 0.0;
  
  static double val1, val2, val3, val4;
  static double val1_sens, val2_sens, val3_sens, val4_sens;
  
  val1 = 0.0; val2 = 0.0; val3 = 0.0; val4 = 0.0;
  val1_sens = 0.0; val2_sens = 0.0; val3_sens = 0.0; val4_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    
    case MEMBRANE_MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(),
							param_enum_ID,
							val1_sens);
        this->getPropertyDerivativeForLocalParameter(THICKNESS_2D_ELEM::num(),
						     param_enum_ID,
						     val2_sens);
        factor = val1*val2_sens + val1_sens * val2;
      }
      break;
      
    case PLATE_MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(),
							param_enum_ID,
							val1_sens);
        this->getPropertyDerivativeForLocalParameter(THICKNESS_2D_ELEM::num(),
						     param_enum_ID,
						     val2_sens);
        factor = val1*(0.25 * pow(val2, 2) * val2_sens) + val1_sens * pow(val2,3) / 12.0;
      }
      break;

      
    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(ALPHA_EXPANSION::num(), val3);
        property.getPropertyValue(POISSONS_RATIO::num(), val4);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(),
                                                               param_enum_ID,
                                                               val1_sens);
        property.getPropertyDerivativeForLocalParameter(ALPHA_EXPANSION::num(),
                                                    param_enum_ID, val3_sens);
        property.getPropertyDerivativeForLocalParameter(POISSONS_RATIO::num(),
                                                    param_enum_ID, val4_sens);
        this->getPropertyDerivativeForLocalParameter(THICKNESS_2D_ELEM::num(),
                                                 param_enum_ID, val2_sens);
        
        factor = (val1_sens * val2 * val3 + val1 * val2_sens * val3 + 
                  val1 * val2 * val3_sens) / (1.0 - val4) - 
          (val1 * val2 * val3) * val4_sens / (1.0 - val4) / (1.0 - val4) ;
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val3);
        
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(),
                                                               param_enum_ID, val1_sens);
        property.getPropertyDerivativeForLocalParameter(SPECIFIC_HEAT::num(),
                                                               param_enum_ID, val2_sens);
        this->getPropertyDerivativeForLocalParameter(THICKNESS_2D_ELEM::num(),
                                                            param_enum_ID, val3_sens);
        
        factor = val1_sens * val2 * val3 + val1 * val2_sens * val3 + 
          val1 * val2 * val3_sens;
      }
      break;
      
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter
          (THERMAL_CONDUCTIVITY::num(),
           param_enum_ID, val1_sens);
        this->getPropertyDerivativeForLocalParameter
          (THICKNESS_2D_ELEM::num(),
           param_enum_ID, val2_sens);
        
        factor = val1_sens * val2 + val1 * val2_sens;
      }
      break;
      

    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter
          (EMISSIVITY::num(),
           param_enum_ID, val1_sens);
        this->getPropertyDerivativeForLocalParameter
          (THICKNESS_2D_ELEM::num(),
           param_enum_ID, val2_sens);
        
        factor = val1_sens * val2 + val1 * val2_sens;
      }
      break;

      
      case RADIATION_EPSILON_FACTOR_1_ENUM_ID:
      case RADIATION_EPSILON_FACTOR_2_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        property.getPropertyDerivativeForLocalParameter
          (EMISSIVITY::num(), param_enum_ID, val1_sens);
        
        factor = (-1.0 / val1 / val1) * val1_sens;
      }
      break;
      
    
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}



inline 
void
Isotropic2D_ElemDataCard::
getFactorSensitivityForLocalParameter(DenseMatrix<double>& matrix,
                                      const unsigned int factor_enum_ID,
                                      const unsigned int param_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  if (matrix.m() != 3 || matrix.n() != 3)
    matrix.resize(3,3);
  
  static DenseMatrix<double> base_qty(3,3);
  
  base_qty.zero();
  matrix.zero();
  
  static double val1, val2, val3;
  static double val1_sens, val2_sens, val3_sens;
  
  val1 = 0.0; val2 = 0.0; val3 = 0.0; 
  val1_sens = 0.0; val2_sens = 0.0; val3_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    
    case STIFFNESS_A_MATRIX_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(),
                                                               param_enum_ID,
                                                               val1_sens);
        property.getPropertyDerivativeForLocalParameter(POISSONS_RATIO::num(),
                                                               param_enum_ID, 
                                                               val3_sens);
        this->getPropertyDerivativeForLocalParameter(THICKNESS_2D_ELEM::num(),
                                                            param_enum_ID, 
                                                            val2_sens);
        
        
        matrix(0,0) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(1,1) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(0,1) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(1,0) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(2,2) = val1_sens / 2.0 / (1.0 + val3) + 
          val1 * val3_sens / 2.0 / (1.0 + val3); // E / (2 *(1+nu))
        
        matrix.scale(val2); // multiply by thickness
        this->getFactor(base_qty, factor_enum_ID);
        matrix.add(val2_sens, base_qty);
      };
      break;
      
    case STIFFNESS_B_MATRIX_FACTOR_ENUM_ID:
      {
        // nothing to be done here, since for isotropic materials, this is zero
      };
      break;
      
    case STIFFNESS_D_MATRIX_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        this->getPropertyValue(THICKNESS_2D_ELEM::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(),
                                                               param_enum_ID,
                                                               val1_sens);
        property.getPropertyDerivativeForLocalParameter(POISSONS_RATIO::num(),
                                                               param_enum_ID,
                                                               val3_sens);
        this->getPropertyDerivativeForLocalParameter(THICKNESS_2D_ELEM::num(),
                                                            param_enum_ID,
                                                            val2_sens);
        
        
        matrix(0,0) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(1,1) = val1_sens / (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E / (1-nu * nu)
        matrix(0,1) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(1,0) = (val1_sens * val3 + val1 * val3_sens )/ (1.0 - val3 * val3) - 
          val1 * 2.0 * val3 * val3_sens / (1.0 - val3 * val3) / 
          (1.0 - val3 * val3); // E * nu / (1-nu * nu)
        matrix(2,2) = val1_sens / 2.0 / (1.0 + val3) + 
          val1 * val3_sens / 2.0 / (1.0 + val3); // E / (2 *(1+nu))
        matrix.scale(val2 * val2 * val2 / 12.0); // multiply by h^3 / 12
        this->getFactor(base_qty, factor_enum_ID);
        matrix.add(val2 * val2 * val2_sens / 4.0, base_qty);
      };
      break;
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}



inline
std::string
Isotropic2D_ElemDataCard::getPropertyEnumName(const unsigned int enum_ID) const
{
  return Isotropic_2D_ElemDataEnum::enumName(enum_ID);
}



inline
const unsigned int
Isotropic2D_ElemDataCard::getPropertyEnumID(const std::string& enum_name) const
{
  return Isotropic_2D_ElemDataEnum::enumID(enum_name);
}


#endif // __fesystem_isotropic_2d_elem_data_card_h__



//
//class Lamina: public PropertyCardBase<double, LaminaDataPropertyEnumType>
//{
//public:
//  Lamina();
//  
//  ~Lamina();
//  
//  unsigned int getMaterialID() const;
//  
//  /// @returns the property card kind enumeration ID
//  virtual const unsigned int getPropertyCardKindEnumID() const;
//  
//  /// @returns the property card kind enumeration name
//  virtual const std::string getPropertyCardKindEnumName() const;
//  
//protected:
//    
//    unsigned int material_id;
//  
//  std::istream& readFromInputStream(std::istream& input);
//};
//
//
//
//
//class LaminatedPlate: public ElemDataCard
//{
// public:
//  LaminatedPlate();
//  
//  virtual ~LaminatedPlate();
//
//  /// @returns the property card kind enumeration ID
//  virtual const unsigned int getPropertyCardKindEnumID() const;
//  
//  /// @returns the property card kind enumeration name
//  virtual const std::string getPropertyCardKindEnumName() const;
//  
//  const DenseMatrix<double>& getMatrix(unsigned int enum_id) const;
//  
//  const DenseMatrix<double>& getMatrix(const std::string& enum_name) const;
//
//  const DenseMatrix<double>& getMatrixSensitivity(unsigned int enum_id, 
//						  const std::string& param_name) const;
//  
//  const DenseMatrix<double>& getMatrixSensitivity(const std::string& enum_name,
//						  const std::string& param_name) const;
//
// protected:
//
//  std::istream& readFromInputStream(std::string& input);
//
//  typedef std::map<unsigned int, Lamina*> LaminaMap;
//
//  Laminate::LaminaMap id_to_lamina_map;
//};
//
  
  
