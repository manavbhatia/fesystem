// $Id: ElemDataCard.h,v 1.9 2006-12-30 02:11:03 manav Exp $

#ifndef __fesystem_elem_data_card_h__
#define __fesystem_elem_data_card_h__

// C++ includes
#include <iostream>


// FESystem includes
#include "Properties/PropertyCard.h"
#include "Utilities/NameEnumHandler.h"
#include "Properties/PropertyDatabase.h"


// libmesh includes
#include "numerics/dense_matrix.h"


// Forward declerations
namespace Property
{
  class PropertyDatabase;
}



#ifndef MASS_FACTOR_ENUM_ID
#define MASS_FACTOR_ENUM_ID 6
#else
#error
#endif


#ifndef MASS_FACTOR_ENUM_NAME
#define MASS_FACTOR_ENUM_NAME "MASS"
#else
#error
#endif


#ifndef MEMBRANE_MASS_FACTOR_ENUM_ID
#define MEMBRANE_MASS_FACTOR_ENUM_ID 17
#else
#error
#endif


#ifndef MEMBRANE_MASS_FACTOR_ENUM_NAME
#define MEMBRANE_MASS_FACTOR_ENUM_NAME "MEMBRANE_MASS"
#else
#error
#endif


#ifndef PLATE_MASS_FACTOR_ENUM_ID
#define PLATE_MASS_FACTOR_ENUM_ID 18
#else
#error
#endif


#ifndef PLATE_MASS_FACTOR_ENUM_NAME
#define PLATE_MASS_FACTOR_ENUM_NAME "PLATE_MASS"
#else
#error
#endif


#ifndef THERMAL_EXPANSION_FACTOR_ENUM_ID
#define THERMAL_EXPANSION_FACTOR_ENUM_ID 7
#else
#error
#endif


#ifndef THERMAL_EXPANSION_FACTOR_ENUM_NAME
#define THERMAL_EXPANSION_FACTOR_ENUM_NAME "THERMAL_EXPANSION"
#else
#error
#endif


#ifndef THERMAL_CAPACITANCE_FACTOR_ENUM_ID
#define THERMAL_CAPACITANCE_FACTOR_ENUM_ID 8
#else
#error
#endif


#ifndef THERMAL_CAPACITANCE_FACTOR_ENUM_NAME
#define THERMAL_CAPACITANCE_FACTOR_ENUM_NAME "THERMAL_CAPACITANCE"
#else
#error
#endif



#ifndef THERMAL_CONDUCTANCE_FACTOR_ENUM_ID
#define THERMAL_CONDUCTANCE_FACTOR_ENUM_ID 9
#else
#error
#endif


#ifndef THERMAL_CONDUCTANCE_FACTOR_ENUM_NAME
#define THERMAL_CONDUCTANCE_FACTOR_ENUM_NAME "THERMAL_CONDUCTANCE"
#else
#error
#endif


#ifndef THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID
#define THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID 10
#else
#error
#endif


#ifndef THERMAL_EMITTED_LOAD_FACTOR_ENUM_NAME
#define THERMAL_EMITTED_LOAD_FACTOR_ENUM_NAME "THERMAL_EMITTED_LOAD"
#else
#error
#endif



#ifndef STIFFNESS_A_MATRIX_FACTOR_ENUM_ID
#define STIFFNESS_A_MATRIX_FACTOR_ENUM_ID 11
#else
#error
#endif


#ifndef STIFFNESS_A_MATRIX_FACTOR_ENUM_NAME
#define STIFFNESS_A_MATRIX_FACTOR_ENUM_NAME "STIFFNESS_A_MATRIX"
#else
#error
#endif

#ifndef STIFFNESS_B_MATRIX_FACTOR_ENUM_ID
#define STIFFNESS_B_MATRIX_FACTOR_ENUM_ID 12
#else
#error
#endif


#ifndef STIFFNESS_B_MATRIX_FACTOR_ENUM_NAME
#define STIFFNESS_B_MATRIX_FACTOR_ENUM_NAME "STIFFNESS_B_MATRIX"
#else
#error
#endif


#ifndef STIFFNESS_D_MATRIX_FACTOR_ENUM_ID
#define STIFFNESS_D_MATRIX_FACTOR_ENUM_ID 13
#else
#error
#endif


#ifndef STIFFNESS_D_MATRIX_FACTOR_ENUM_NAME
#define STIFFNESS_D_MATRIX_FACTOR_ENUM_NAME "STIFFNESS_D_MATRIX"
#else
#error
#endif

#ifndef RADIATION_EPSILON_FACTOR_1_ENUM_ID
#define RADIATION_EPSILON_FACTOR_1_ENUM_ID 14
#else
#error
#endif


#ifndef RADIATION_EPSILON_FACTOR_1_ENUM_NAME
#define RADIATION_EPSILON_FACTOR_1_ENUM_NAME "RADIATION_EPSILON_FACTOR_1"
#else
#error
#endif


#ifndef RADIATION_EPSILON_FACTOR_2_ENUM_ID
#define RADIATION_EPSILON_FACTOR_2_ENUM_ID 15
#else
#error
#endif


#ifndef RADIATION_EPSILON_FACTOR_2_ENUM_NAME
#define RADIATION_EPSILON_FACTOR_2_ENUM_NAME "RADIATION_EPSILON_FACTOR_2"
#else
#error
#endif


#ifndef STRESS_STRAIN_FACTOR_ENUM_ID
#define STRESS_STRAIN_FACTOR_ENUM_ID 16
#else
#error
#endif


#ifndef STRESS_STRAIN_FACTOR_ENUM_NAME
#define STRESS_STRAIN_FACTOR_ENUM_NAME "STRESS_STRAIN"
#else
#error
#endif



// declare a class for element data types
DeclareEnumClass(ElemDataCardEnum);



DeclareEnumClass(ElemDataCardFactorsEnum);

DeclareEnumName(MASS_FACTOR, ElemDataCardFactorsEnum,
                MASS_FACTOR_ENUM_ID, MASS_FACTOR_ENUM_NAME);

DeclareEnumName(MEMBRANE_MASS_FACTOR, ElemDataCardFactorsEnum,
                MEMBRANE_MASS_FACTOR_ENUM_ID, MEMBRANE_MASS_FACTOR_ENUM_NAME);

DeclareEnumName(PLATE_MASS_FACTOR, ElemDataCardFactorsEnum,
                PLATE_MASS_FACTOR_ENUM_ID, PLATE_MASS_FACTOR_ENUM_NAME);

DeclareEnumName(THERMAL_EXPANSION_FACTOR, ElemDataCardFactorsEnum,
                THERMAL_EXPANSION_FACTOR_ENUM_ID, THERMAL_EXPANSION_FACTOR_ENUM_NAME);

DeclareEnumName(THERMAL_CAPACITANCE_FACTOR, ElemDataCardFactorsEnum,
                THERMAL_CAPACITANCE_FACTOR_ENUM_ID, THERMAL_CAPACITANCE_FACTOR_ENUM_NAME);


DeclareEnumName(THERMAL_CONDUCTANCE_FACTOR, ElemDataCardFactorsEnum,
                THERMAL_CONDUCTANCE_FACTOR_ENUM_ID, THERMAL_CONDUCTANCE_FACTOR_ENUM_NAME);

DeclareEnumName(THERMAL_EMITTED_LOAD_FACTOR, ElemDataCardFactorsEnum,
                THERMAL_EMITTED_LOAD_FACTOR_ENUM_ID, THERMAL_EMITTED_LOAD_FACTOR_ENUM_NAME);

DeclareEnumName(STIFFNESS_A_MATRIX_FACTOR, ElemDataCardFactorsEnum,
                STIFFNESS_A_MATRIX_FACTOR_ENUM_ID, STIFFNESS_A_MATRIX_FACTOR_ENUM_NAME);

DeclareEnumName(STIFFNESS_B_MATRIX_FACTOR, ElemDataCardFactorsEnum,
                STIFFNESS_B_MATRIX_FACTOR_ENUM_ID, STIFFNESS_B_MATRIX_FACTOR_ENUM_NAME);

DeclareEnumName(STIFFNESS_D_MATRIX_FACTOR, ElemDataCardFactorsEnum,
                STIFFNESS_D_MATRIX_FACTOR_ENUM_ID, STIFFNESS_D_MATRIX_FACTOR_ENUM_NAME);

DeclareEnumName(RADIATION_EPSILON_FACTOR_1, ElemDataCardFactorsEnum,
                RADIATION_EPSILON_FACTOR_1_ENUM_ID, RADIATION_EPSILON_FACTOR_1_ENUM_NAME);

DeclareEnumName(RADIATION_EPSILON_FACTOR_2, ElemDataCardFactorsEnum,
                RADIATION_EPSILON_FACTOR_2_ENUM_ID, RADIATION_EPSILON_FACTOR_2_ENUM_NAME);

DeclareEnumName(STRESS_STRAIN_FACTOR, ElemDataCardFactorsEnum,
                STRESS_STRAIN_FACTOR_ENUM_ID, STRESS_STRAIN_FACTOR_ENUM_NAME);




class ElemDataCard: public Property::PropertyCard<double>
{
public: 
  ElemDataCard();
  
  ~ElemDataCard();
  
  /// attached the property database to this card.
  inline void attachPropertyDatabase(Property::PropertyDatabase& database);
  
  /// initializes the card and the material card for the given local parameters
  virtual void reinitElemAndMaterialCardForLocalParameters
    (std::map<unsigned int, double> *value_map = NULL) = 0;
  
  /// initializes the card and the material card for the given global parameters
  virtual void reinitElemAndMaterialCardForGlobalParameters
    (std::map<unsigned int, double> &value_map) = 0;

  /// clears initialization of the card and the material card for local parameters
  virtual void clearElemAndMaterialCardLocalParameterInitialization() = 0;
  
  /// clears initialization of the card and the material card only for local parameters
  /// specified in the input vector
  virtual void 
    partialClearElemAndMaterialCardLocalParameterInitialization
    (std::vector<unsigned int>& param_ids) = 0;

  /// clears initialization of the card and the material card for global parameters
  virtual void clearElemAndMaterialCardGlobalParameterInitialization() = 0;
  
  /// @returns true if either the element data card, or the material property card are
  /// dependent on the specified global parameter
  virtual bool checkElemAndMaterialCardGlobalParameterDependence
    (const unsigned int global_param_ID) = 0;
  
  /// @returns true if either the element data card, or the material property card are
  /// dependent on the specified local parameter
  virtual bool checkElemAndMaterialCardLocalParameterDependence
    (const unsigned int local_param_ID) = 0;
  
  
  /// returns the value of the factor
  virtual void getFactor(double& factor, const unsigned int factor_enum_ID) =0;
  
  virtual void getFactorSensitivityForGlobalParameter(double& factor,
                                                      const unsigned int factor_enum_ID,
                                                      const unsigned int global_param_ID)=0;
  
  virtual void getFactorSensitivityForLocalParameter(double& factor,
                                                     const unsigned int factor_enum_ID,
                                                     const unsigned int local_param_enum_ID)=0;

  /// returns the value of the factor
  virtual void getFactor(DenseMatrix<double>& factor, const unsigned int factor_enum_ID) =0;
  
  virtual void getFactorSensitivityForGlobalParameter(DenseMatrix<double>& factor,
                                                      const unsigned int factor_enum_ID,
                                                      const unsigned int global_param_ID)=0;
  
  virtual void getFactorSensitivityForLocalParameter(DenseMatrix<double>& factor,
                                                     const unsigned int factor_enum_ID,
                                                     const unsigned int local_param_enum_ID)=0;
  
  
  virtual  std::istream& readFromInputStream(std::istream& input) = 0;
  
  /// returns value from material card
  virtual void 
    getPropertyValueFromMaterialCard(const unsigned int prop_enum_ID, double& value,
                                     const unsigned int layer = 0) = 0;
  
  /// returns value of the derivative of a property from the material card
  virtual void
    getPropertyValueDerivativeForGlobalParameterFromMaterialCard
    (const unsigned int prop_enum_ID, 
     const unsigned int global_param_ID,
     double& value,
     const unsigned int layer = 0) = 0;

  /// returns value of the derivative of a property from the material card
  virtual void
    getPropertyValueDerivativeForLocalParameterFromMaterialCard
    (const unsigned int prop_enum_ID, 
     const unsigned int local_param_ID,
     double& value,
     const unsigned int layer = 0) = 0;
  
protected:

  Property::PropertyDatabase *property_database;
};





inline 
void
ElemDataCard::attachPropertyDatabase(Property::PropertyDatabase &database)
{
  this->property_database = &database;
}





#endif // __fesystem_elem_data_card_h__

