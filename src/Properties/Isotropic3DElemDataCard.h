// $Id: Isotropic3DElemDataCard.h,v 1.8.6.1 2008-06-03 05:19:52 manav Exp $

#ifndef __fesystem_isotropic_3d_elem_data_card_h__
#define __fesystem_isotropic_3d_elem_data_card_h__

// FESystem includes
#include "Properties/IsotropicElemDataCard.h"
#include "Properties/MaterialPropertyNameEnums.h"


#ifndef ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_ID
#define ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_ID 5
#else
#error
#endif


#ifndef ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_NAME
#define ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_NAME "ISOTROPIC_3D_ELEM_DATA_CARD"
#else
#error
#endif


#ifndef SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_ID
#define SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_ID 19
#else
#error
#endif


#ifndef SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_NAME
#define SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_NAME "SOLID_STIFFNESS_MATRIX_FACTOR"
#else
#error
#endif


DeclareEnumName(ISOTROPIC_3D_ELEM_DATA_CARD,
                ElemDataCardEnum,
                ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_ID,
                ISOTROPIC_3D_ELEM_DATA_CARD_ENUM_NAME);

DeclareEnumClass(Isotropic_3D_ElemDataEnum);


DeclareEnumName(SOLID_STIFFNESS_MATRIX_FACTOR, ElemDataCardFactorsEnum,
                SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_ID, SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_NAME);


class Isotropic3D_ElemDataCard: public IsotropicElemDataCard
{
public:
  
  Isotropic3D_ElemDataCard();
  
  virtual  ~Isotropic3D_ElemDataCard();
  
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
  
  
  friend std::istream& operator>> (std::istream& input, Isotropic3D_ElemDataCard& plane_elem);
protected:
    
    /// @returns the property name from enumeration ID of the property
    virtual std::string getPropertyEnumName(const unsigned int enum_ID) const;
  
  /// @returns the property ID from enumeration name of the property
  virtual const unsigned int getPropertyEnumID(const std::string& enum_name) const; 
};




inline
unsigned int
Isotropic3D_ElemDataCard::getPropertyCardKindEnumID() const
{
  return ISOTROPIC_3D_ELEM_DATA_CARD::num();
}


inline
const std::string
Isotropic3D_ElemDataCard::getPropertyCardKindEnumName() const
{
  return ISOTROPIC_3D_ELEM_DATA_CARD::name();
}






inline 
void
Isotropic3D_ElemDataCard::getFactor(double& factor,
                                    const unsigned int factor_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  factor = 0.0;
  
  static double val1, val2;
  val1 = 0.0; val2 = 0.0;

  switch (factor_enum_ID)
    {
    
    case MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        factor = val1;
      }
      break;
      
    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        static double E, alpha, nu;
        property.getPropertyValue(YOUNGS_MODULUS::num(), E);
        property.getPropertyValue(ALPHA_EXPANSION::num(), alpha);
        property.getPropertyValue(POISSONS_RATIO::num(), nu);
        
        factor =  E*alpha*(1.0+nu)/(1.0-nu-2.0*nu*nu);
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        
        factor = val1 * val2;
      }
      break;
      
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        
        factor = val1;
      }
      break;
      
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(EMISSIVITY::num(), val1);
        
        factor = val1;
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
Isotropic3D_ElemDataCard::getFactor(DenseMatrix<double>& matrix,
                                    const unsigned int factor_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  if (matrix.m() != 6 || matrix.n() != 6)
    matrix.resize(6,6);
  
  matrix.zero();
  
  static double val1, val2, val3, val4;
  val1 = 0.0; val2 = 0.0; val3 = 0.0; val4 = 0.0;
  
  switch (factor_enum_ID)
    {
      
    case SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_ID:
      {
        static double E, nu;

        property.getPropertyValue(YOUNGS_MODULUS::num(), E);
        property.getPropertyValue(POISSONS_RATIO::num(), nu);
        
        matrix(0,0) = E * (1.0-nu) / (1-nu-2.0*nu*nu);
        matrix(1,1) = matrix(0,0);
        matrix(2,2) = matrix(0,0);
        matrix(0,1) = E * nu / (1-nu-2.0*nu*nu);
        matrix(0,2) = matrix(0,1);
        matrix(1,0) = matrix(0,1);
        matrix(1,2) = matrix(0,1);
        matrix(2,0) = matrix(0,1);
        matrix(2,1) = matrix(0,1);
        matrix(3,3) = E / 2.0 / (1.0 + nu); // E / (2 *(1+nu))
        matrix(4,4) = matrix(3,3);
        matrix(5,5) = matrix(3,3);
      };
      break;
            
      case MASS_FACTOR_ENUM_ID:
      default:
        Assert(false,
               FESystemExceptions::ExcEnumCaseNotImplemented
               (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}




inline 
void
Isotropic3D_ElemDataCard::
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
    
    case MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(), 
                                                         global_param_ID,
                                                         val1_sens);
        factor = val1_sens;
      }
      break;
      
    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        static double E, alpha, nu, E_sens, alpha_sens, nu_sens;
        
        property.getPropertyValue(YOUNGS_MODULUS::num(), E);
        property.getPropertyValue(ALPHA_EXPANSION::num(), alpha);
        property.getPropertyValue(POISSONS_RATIO::num(), nu);
        
        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(),
                                                         global_param_ID,
                                                         E_sens);
        property.getPropertyDerivativeForGlobalParameter(ALPHA_EXPANSION::num(),
                                                         global_param_ID, alpha_sens);
        property.getPropertyDerivativeForGlobalParameter(POISSONS_RATIO::num(),
                                                         global_param_ID, nu_sens);
        
        factor = (E_sens * alpha + E * alpha_sens)*(1.0+nu)/(1.0-nu-2.0*nu*nu) + 
          nu_sens *  E*alpha*(2.0+4.0*nu+2.0*nu*nu)/pow((1.0-nu-2.0*nu*nu),2);
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        
        property.getPropertyDerivativeForGlobalParameter(DENSITY::num(),
                                                         global_param_ID, val1_sens);
        property.getPropertyDerivativeForGlobalParameter(SPECIFIC_HEAT::num(),
                                                         global_param_ID, val2_sens);
        
        factor = val1_sens * val2 * + val1 * val2_sens;
      }
      break;
      
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyDerivativeForGlobalParameter
          (THERMAL_CONDUCTIVITY::num(),
           global_param_ID, val1_sens);

        factor = val1_sens;
      }
      break;
      
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyDerivativeForGlobalParameter
        (EMISSIVITY::num(),
         global_param_ID, val1_sens);
        
        factor = val1_sens;
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
Isotropic3D_ElemDataCard::
getFactorSensitivityForGlobalParameter(DenseMatrix<double>& matrix, 
                                       const unsigned int factor_enum_ID,
                                       const unsigned int global_param_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  if (matrix.m() != 6 || matrix.n() != 6)
    matrix.resize(6,6);
  
  matrix.zero();
  
  switch (factor_enum_ID)
    {
    
    case SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_ID:
      {
        static double E, nu, E_sens, nu_sens;
        
        property.getPropertyValue(YOUNGS_MODULUS::num(), E);
        property.getPropertyValue(POISSONS_RATIO::num(), nu);

        property.getPropertyDerivativeForGlobalParameter(YOUNGS_MODULUS::num(), 
                                                         global_param_ID,
                                                         E_sens);
        property.getPropertyDerivativeForGlobalParameter(POISSONS_RATIO::num(),
                                                         global_param_ID,
                                                         nu_sens);

        
        matrix(0,0) = E_sens * (1.0-nu) / (1-nu-2.0*nu*nu) + 
          nu_sens * E * (-2.0*nu + 4.0) * nu / pow((1-nu-2.0*nu*nu),2);
        matrix(1,1) = matrix(0,0);
        matrix(2,2) = matrix(0,0);
        matrix(0,1) = E_sens * nu / (1-nu-2.0*nu*nu) + 
          nu_sens * E * (1.0+2.0*nu*nu) / pow((1-nu-2.0*nu*nu),2);
        matrix(0,2) = matrix(0,1);
        matrix(1,0) = matrix(0,1);
        matrix(1,2) = matrix(0,1);
        matrix(2,0) = matrix(0,1);
        matrix(2,1) = matrix(0,1);
        matrix(3,3) = E_sens / 2.0 / (1.0 + nu) - 
          nu_sens * -0.5 * E / pow((1.0 + nu),2);
        matrix(4,4) = matrix(3,3);
        matrix(5,5) = matrix(3,3);
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
Isotropic3D_ElemDataCard::
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
    
    case MASS_FACTOR_ENUM_ID:
      {
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(), 
                                                        param_enum_ID,
                                                        val1_sens);
        factor = val1_sens;
      }
      break;
      
    case THERMAL_EXPANSION_FACTOR_ENUM_ID:
      {
        static double E, alpha, nu, E_sens, alpha_sens, nu_sens;
        
        property.getPropertyValue(YOUNGS_MODULUS::num(), E);
        property.getPropertyValue(ALPHA_EXPANSION::num(), alpha);
        property.getPropertyValue(POISSONS_RATIO::num(), nu);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(),
                                                         param_enum_ID,
                                                         E_sens);
        property.getPropertyDerivativeForLocalParameter(ALPHA_EXPANSION::num(),
                                                         param_enum_ID, alpha_sens);
        property.getPropertyDerivativeForLocalParameter(POISSONS_RATIO::num(),
                                                         param_enum_ID, nu_sens);
        
        factor = (E_sens * alpha + E * alpha_sens)*(1.0+nu)/(1.0-nu-2.0*nu*nu) + 
          nu_sens *  E*alpha*(2.0+4.0*nu+2.0*nu*nu)/pow((1.0-nu-2.0*nu*nu),2);
      }
      break;
      
    case THERMAL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyValue(DENSITY::num(), val1);
        property.getPropertyValue(SPECIFIC_HEAT::num(), val2);
        
        property.getPropertyDerivativeForLocalParameter(DENSITY::num(),
                                                        param_enum_ID, val1_sens);
        property.getPropertyDerivativeForLocalParameter(SPECIFIC_HEAT::num(),
                                                        param_enum_ID, val2_sens);

        factor = val1_sens * val2 + val1 * val2_sens;
      }
      break;
      
      
    case THERMAL_CONDUCTANCE_FACTOR_ENUM_ID:
      {
        property.getPropertyDerivativeForLocalParameter
          (THERMAL_CONDUCTIVITY::num(),
           param_enum_ID, val1_sens);

        factor = val1_sens;
      }
      break;
      
      
    case THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID:
      {
        property.getPropertyDerivativeForLocalParameter
        (EMISSIVITY::num(), param_enum_ID, val1_sens);
        
        factor = val1_sens;
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
Isotropic3D_ElemDataCard::
getFactorSensitivityForLocalParameter(DenseMatrix<double>& matrix,
                                      const unsigned int factor_enum_ID,
                                      const unsigned int param_enum_ID)
{
  Property::PropertyCard<double>& property =this->getMaterialCard();
  
  if (matrix.m() != 6 || matrix.n() != 6)
    matrix.resize(6,6);
  
  matrix.zero();
    
  switch (factor_enum_ID)
    {
          
    case SOLID_STIFFNESS_MATRIX_FACTOR_ENUM_ID:
      {
        static double E, nu, E_sens, nu_sens;
        
        property.getPropertyValue(YOUNGS_MODULUS::num(), E);
        property.getPropertyValue(POISSONS_RATIO::num(), nu);
        
        property.getPropertyDerivativeForLocalParameter(YOUNGS_MODULUS::num(), 
                                                         param_enum_ID,
                                                         E_sens);
        property.getPropertyDerivativeForLocalParameter(POISSONS_RATIO::num(),
                                                         param_enum_ID,
                                                         nu_sens);
        
        
        matrix(0,0) = E_sens * (1.0-nu) / (1-nu-2.0*nu*nu) + 
          nu_sens * E * (-2.0*nu + 4.0) * nu / pow((1-nu-2.0*nu*nu),2);
        matrix(1,1) = matrix(0,0);
        matrix(2,2) = matrix(0,0);
        matrix(0,1) = E_sens * nu / (1-nu-2.0*nu*nu) + 
          nu_sens * E * (1.0+2.0*nu*nu) / pow((1-nu-2.0*nu*nu),2);
        matrix(0,2) = matrix(0,1);
        matrix(1,0) = matrix(0,1);
        matrix(1,2) = matrix(0,1);
        matrix(2,0) = matrix(0,1);
        matrix(2,1) = matrix(0,1);
        matrix(3,3) = E_sens / 2.0 / (1.0 + nu) - 
          nu_sens * -0.5 * E / pow((1.0 + nu),2);
        matrix(4,4) = matrix(3,3);
        matrix(5,5) = matrix(3,3);
      };
      break;
      
      case MASS_FACTOR_ENUM_ID:
      default:
        Assert(false,
               FESystemExceptions::ExcEnumCaseNotImplemented
               (ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}



inline
std::string
Isotropic3D_ElemDataCard::getPropertyEnumName(const unsigned int enum_ID) const
{
  return Isotropic_3D_ElemDataEnum::enumName(enum_ID);
}



inline
const unsigned int
Isotropic3D_ElemDataCard::getPropertyEnumID(const std::string& enum_name) const
{
  return Isotropic_3D_ElemDataEnum::enumID(enum_name);
}



#endif // __fesystem_isotropic_3d_elem_data_card_h__



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
