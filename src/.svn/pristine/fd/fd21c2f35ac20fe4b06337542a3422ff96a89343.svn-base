// $Id: Isotropic1DElemDataCard.h,v 1.6 2006-11-13 00:04:44 manav Exp $

#ifndef __fesystem_isotropic_1d_elem_data_card_h__
#define __fesystem_isotropic_1d_elem_data_card_h__

// FESystem includes
#include "Properties/IsotropicElemDataCard.h"
#include "Properties/MaterialPropertyNameEnums.h"


#ifndef ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_ID
#define ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_ID 3
#else
#error
#endif

#ifndef ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_NAME
#define ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_NAME "ISOTROPIC_1D_ELEM_DATA_CARD"
#else
#error
#endif

#ifndef AREA_ENUM_ID
#define AREA_ENUM_ID  1
#else
#error
#endif

#ifndef AREA_ENUM_NAME
#define AREA_ENUM_NAME "AREA"
#else
#error
#endif

#ifndef THICKNESS_1D_ELEM_ENUM_ID
#define THICKNESS_1D_ELEM_ENUM_ID 2
#else
#error
#endif

#ifndef THICKNESS_1D_ELEM_ENUM_NAME
#define THICKNESS_1D_ELEM_ENUM_NAME  "THICKNESS"
#else
#error
#endif

#ifndef IZZ_ENUM_ID
#define IZZ_ENUM_ID 3
#else
#error
#endif

#ifndef IZZ_ENUM_NAME
#define IZZ_ENUM_NAME "IZZ"
#else
#error
#endif

#ifndef IYY_ENUM_ID
#define IYY_ENUM_ID 4
#else
#error
#endif

#ifndef IYY_ENUM_NAME
#define IYY_ENUM_NAME "IYY"
#else
#error
#endif

#ifndef ROTATION_AXIS_X_ENUM_ID
#define ROTATION_AXIS_X_ENUM_ID 5
#else
#error
#endif

#ifndef ROTATION_AXIS_X_ENUM_NAME
#define ROTATION_AXIS_X_ENUM_NAME "ROTATION_AXIS_X"
#else
#error
#endif

#ifndef ROTATION_AXIS_Y_ENUM_ID
#define ROTATION_AXIS_Y_ENUM_ID 6
#else
#error
#endif

#ifndef ROTATION_AXIS_Y_ENUM_NAME
#define ROTATION_AXIS_Y_ENUM_NAME "ROTATION_AXIS_Y"
#else
#error
#endif


#ifndef ROTATION_AXIS_Z_ENUM_ID
#define ROTATION_AXIS_Z_ENUM_ID 7
#else
#error
#endif

#ifndef ROTATION_AXIS_Z_ENUM_NAME
#define ROTATION_AXIS_Z_ENUM_NAME "ROTATION_AXIS_Z"
#else
#error
#endif


// factor = E * A
#ifndef EA_FACTOR_ENUM_ID
#define EA_FACTOR_ENUM_ID 1
#else
#error
#endif

#ifndef EA_FACTOR_ENUM_NAME
#define EA_FACTOR_ENUM_NAME "EA"
#else
#error
#endif

// factor = G * A
#ifndef GA_FACTOR_ENUM_ID
#define GA_FACTOR_ENUM_ID 2
#else
#error
#endif

#ifndef GA_FACTOR_ENUM_NAME
#define GA_FACTOR_ENUM_NAME "GA"
#else
#error
#endif

// factor = E * IZZ
#ifndef EIZZ_FACTOR_ENUM_ID
#define EIZZ_FACTOR_ENUM_ID 3
#else
#error
#endif

#ifndef EIZZ_FACTOR_ENUM_NAME
#define EIZZ_FACTOR_ENUM_NAME "EIZZ"
#else
#error
#endif

// factor = E * IYY
#ifndef EIYY_FACTOR_ENUM_ID
#define EIYY_FACTOR_ENUM_ID 4
#else
#error
#endif

#ifndef EIYY_FACTOR_ENUM_NAME
#define EIYY_FACTOR_ENUM_NAME "EIYY"
#else
#error
#endif

#ifndef SPRING_STIFFNESS_FACTOR_ENUM_ID
#define SPRING_STIFFNESS_FACTOR_ENUM_ID 5
#else
#error
#endif

#ifndef SPRING_STIFFNESS_FACTOR_ENUM_NAME
#define SPRING_STIFFNESS_FACTOR_ENUM_NAME "STIFFNESS_FACTOR"
#else
#error
#endif


DeclareEnumName(ISOTROPIC_1D_ELEM_DATA_CARD,
                ElemDataCardEnum,
                ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_ID,
                ISOTROPIC_1D_ELEM_DATA_CARD_ENUM_NAME);


DeclareEnumClass(Isotropic_1D_ElemDataEnum);



DeclareEnumName(AREA, Isotropic_1D_ElemDataEnum, AREA_ENUM_ID, AREA_ENUM_NAME);



DeclareEnumName(THICKNESS_1D_ELEM, Isotropic_1D_ElemDataEnum, 
                THICKNESS_1D_ELEM_ENUM_ID , THICKNESS_1D_ELEM_ENUM_NAME);


DeclareEnumName(IZZ, Isotropic_1D_ElemDataEnum,
                IZZ_ENUM_ID, IZZ_ENUM_NAME);

DeclareEnumName(IYY, Isotropic_1D_ElemDataEnum,
                IYY_ENUM_ID, IYY_ENUM_NAME);


DeclareEnumName(ROTATION_AXIS_X, Isotropic_1D_ElemDataEnum,
                ROTATION_AXIS_X_ENUM_ID, 
                ROTATION_AXIS_X_ENUM_NAME);

DeclareEnumName(ROTATION_AXIS_Y, Isotropic_1D_ElemDataEnum,
                ROTATION_AXIS_Y_ENUM_ID, 
                ROTATION_AXIS_Y_ENUM_NAME);

DeclareEnumName(ROTATION_AXIS_Z, Isotropic_1D_ElemDataEnum,
                ROTATION_AXIS_Z_ENUM_ID, 
                ROTATION_AXIS_Z_ENUM_NAME);


DeclareEnumName(EA_FACTOR, ElemDataCardFactorsEnum,
                EA_FACTOR_ENUM_ID, EA_FACTOR_ENUM_NAME);

DeclareEnumName(GA_FACTOR, ElemDataCardFactorsEnum,
                GA_FACTOR_ENUM_ID, GA_FACTOR_ENUM_NAME);

DeclareEnumName(EIZZ_FACTOR, ElemDataCardFactorsEnum,
                EIZZ_FACTOR_ENUM_ID, EIZZ_FACTOR_ENUM_NAME);

DeclareEnumName(EIYY_FACTOR, ElemDataCardFactorsEnum,
                EIYY_FACTOR_ENUM_ID, EIYY_FACTOR_ENUM_NAME);

DeclareEnumName(SPRING_STIFFNESS_FACTOR, ElemDataCardFactorsEnum,
                SPRING_STIFFNESS_FACTOR_ENUM_ID, SPRING_STIFFNESS_FACTOR_ENUM_NAME);




class Isotropic1D_ElemDataCard: public IsotropicElemDataCard
{
public:
  Isotropic1D_ElemDataCard();
  
  virtual  ~Isotropic1D_ElemDataCard();
  
  
  /// @returns the property card kind enumeration ID
  virtual unsigned int getPropertyCardKindEnumID() const;
  
  /// @returns the property card kind enumeration name
  virtual const std::string getPropertyCardKindEnumName() const;
  
  /// returns the value of the factor
  virtual void getFactor(double& factor, const unsigned int factor_enum_ID) ;
  
  virtual void getFactorSensitivityForGlobalParameter(double& factor,
                                                      const unsigned int factor_enum_ID,
                                                      const unsigned int global_param_ID);
  
  virtual void getFactorSensitivityForLocalParameter(double& factor,
                                                     const unsigned int factor_enum_ID,
                                                     const unsigned int local_param_enum_ID);
  
  /// returns the value of the factor
  virtual void getFactor(DenseMatrix<double>& factor, const unsigned int factor_enum_ID) ;
  
  virtual void getFactorSensitivityForGlobalParameter(DenseMatrix<double>& factor,
                                                      const unsigned int factor_enum_ID,
                                                      const unsigned int global_param_ID);
  
  virtual void getFactorSensitivityForLocalParameter(DenseMatrix<double>& factor,
                                                     const unsigned int factor_enum_ID,
                                                     const unsigned int local_param_enum_ID);
  
  friend std::istream& operator>> (std::istream& input, Isotropic1D_ElemDataCard& line_elem);


protected:
    
  /// @returns the property name from enumeration ID of the property
  virtual std::string getPropertyEnumName(const unsigned int enum_ID) const;
  
  /// @returns the property ID from enumeration name of the property
  virtual const unsigned int getPropertyEnumID(const std::string& enum_name) const; 
};


inline
unsigned int
Isotropic1D_ElemDataCard::getPropertyCardKindEnumID() const
{
  return ISOTROPIC_1D_ELEM_DATA_CARD::num();
}





inline
const std::string
Isotropic1D_ElemDataCard::getPropertyCardKindEnumName() const
{
  return ISOTROPIC_1D_ELEM_DATA_CARD::name();
}




inline 
void
Isotropic1D_ElemDataCard::getFactor(double& factor,
                                    const unsigned int factor_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  factor = 0.0;
  static double val1, val2, val3;
  val1 = 0.0; val2 = 0.0; val3 = 0.0;
  
  switch (factor_enum_ID)
    {
    case EA_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        factor = val1*val2;
      }
      break;
      
    case MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        factor = val1*val2;
      }
      break;

    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(ALPHA_EXPANSION::num(), val3);
        this->getPropertyValue(AREA::num(), val2);
        factor = val1*val2*val3;
      }
      break;

    case GA_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        
        this->getPropertyValue(AREA::num(), val2);
        factor = (val1/ 2.0 / (1.0 + val3) ) * val2;
      }
      break;

    case EIZZ_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(IZZ::num(), val2);
        factor = val1*val2;
      }
      break;

    case EIYY_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(IYY::num(), val2);
        factor = val1*val2;
      }
      break;

    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        this->getPropertyValue(AREA::num(), val3);
        
        factor = val1 * val2 * val3;
      }
      break;

    
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        factor = val1 * val2;
      }
      break;
      
      
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        factor = val1 * val2;
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
Isotropic1D_ElemDataCard::getFactorSensitivityForGlobalParameter
(double& factor,
 const unsigned int factor_enum_ID,
 const unsigned int global_param_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  factor = 0.0;
  static double val1, val2, val3, val1_sens, val2_sens, val3_sens;
  val1 = 0.0; val2 = 0.0; val3 = 0.0; val1_sens = 0.0; val2_sens = 0.0; val3_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    case EA_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(AREA::num(), val2);

        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);

        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;

    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(ALPHA_EXPANSION::num(), val3);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        property.getPropertyDerivativeForGlobalParameter(ALPHA_EXPANSION::num(), 
                                                         global_param_ID,
                                                         val3_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        
        factor = val1_sens * val2 * val3 + 
          val1 * val2_sens * val3 + val1 * val2 * val3_sens ;
      }
      break;

    case GA_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        property.getPropertyDerivativeForGlobalParameter(POISSONS_RATIO::num(), 
                                                         global_param_ID,
                                                         val3_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        
        static double G, G_sens;
        G = val1 / 2.0 / (1.0 + val3);
        G_sens = (val1_sens * (1.0 + val3) - val1 * val3_sens)/ 
          (1.0 + val3) / (1.0 + val3);
        factor = G * val2_sens + G_sens * val2;
      }
      break;
      
    case EIZZ_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(IZZ::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(IZZ::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case EIYY_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(IYY::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(IYY::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(SPECIFIC_HEAT::num(), val1);
        property.getPropertyValue(DENSITY::num(), val2);
        this->getPropertyValue(AREA::num(), val3);
        
        property.getPropertyDerivativeForGlobalParameter(SPECIFIC_HEAT::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(), 
                                                         global_param_ID,
                                                         val2_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val3_sens);
        
        factor = val1_sens * val2 * val3 + 
          val1 * val2_sens * val3 + val1 * val2 * val3_sens;
      }
      break;

    
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(THERMAL_CONDUCTIVITY::num(), 
                                                        global_param_ID,
                                                        val1_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(EMISSIVITY::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForGlobalParameter(AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
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
Isotropic1D_ElemDataCard::getFactorSensitivityForLocalParameter
(double& factor,
 const unsigned int factor_enum_ID,
 const unsigned int param_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  factor = 0.0;
  static double val1, val2, val3, val1_sens, val2_sens, val3_sens;
  val1 = 0.0; val2 = 0.0; val3 = 0.0; val1_sens = 0.0; val2_sens = 0.0; val3_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    case EA_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(), 
                                                         param_enum_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                      param_enum_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(), 
                                                        param_enum_ID,
                                                        val1_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                     param_enum_ID,
                                                     val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;

    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(ALPHA_EXPANSION::num(), val3);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(), 
                                                         param_enum_ID,
                                                         val1_sens);
        property.getPropertyDerivativeForLocalParameter(ALPHA_EXPANSION::num(), 
                                                         param_enum_ID,
                                                         val3_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                      param_enum_ID,
                                                      val2_sens);
        
        factor = val1_sens * val2 * val3 + 
          val1 * val2_sens * val3 + val1 * val2 * val3_sens ;
      }
      break;

      
    case GA_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        property.getPropertyValue(POISSONS_RATIO::num(), val3);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(), 
                                                         param_enum_ID,
                                                         val1_sens);
        property.getPropertyDerivativeForLocalParameter(POISSONS_RATIO::num(), 
                                                         param_enum_ID,
                                                         val3_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                      param_enum_ID,
                                                      val2_sens);
        
        static double G, G_sens;
        G = val1 / 2.0 / (1.0 + val3);
        G_sens = (val1_sens * (1.0 + val3) - val1 * val3_sens)/ 
          (1.0 + val3) / (1.0 + val3);
        factor = G * val2_sens + G_sens * val2;
      }
      break;      
      
    case EIZZ_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(IZZ::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(), 
                                                         param_enum_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForLocalParameter(IZZ::num(), 
                                                      param_enum_ID,
                                                      val2_sens);
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case EIYY_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(YOUNGS_MODULUS::num(), val1);
        this->getPropertyValue(IYY::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(), 
                                                         param_enum_ID,
                                                         val1_sens);
        this->getPropertyDerivativeForLocalParameter(IYY::num(), 
                                                      param_enum_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;

    
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(SPECIFIC_HEAT::num(), val1);
        property.getPropertyValue(DENSITY::num(), val2);
        this->getPropertyValue(AREA::num(), val3);
        
        property.getPropertyDerivativeForLocalParameter(SPECIFIC_HEAT::num(), 
                                                         param_enum_ID,
                                                         val1_sens);
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(), 
                                                         param_enum_ID,
                                                         val2_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                      param_enum_ID,
                                                      val3_sens);
        
        factor = val1_sens * val2 * val3 + 
          val1 * val2_sens * val3 + val1 * val2 * val3_sens;
      }
      break;
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(THERMAL_CONDUCTIVITY::num(), 
                                                        param_enum_ID,
                                                        val1_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                     param_enum_ID,
                                                     val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;

    
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        this->getPropertyValue(AREA::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(EMISSIVITY::num(), 
                                                        param_enum_ID,
                                                        val1_sens);
        this->getPropertyDerivativeForLocalParameter(AREA::num(), 
                                                     param_enum_ID,
                                                     val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}





inline
std::string
Isotropic1D_ElemDataCard::getPropertyEnumName(const unsigned int enum_ID) const
{
  return Isotropic_1D_ElemDataEnum::enumName(enum_ID);
}



inline
const unsigned int
Isotropic1D_ElemDataCard::getPropertyEnumID(const std::string& enum_name) const
{
  return Isotropic_1D_ElemDataEnum::enumID(enum_name);
}



#endif // __fesystem_isotropic_1d_elem_data_card_h__
