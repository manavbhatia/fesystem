// $Id: ForcedConvection1DElemDataCard.h,v 1.3.6.2 2007-07-19 06:30:38 manav Exp $

#ifndef __fesystem_forced_convection_1d_elem_data_card_h__
#define __fesystem_forced_convection_1d_elem_data_card_h__

// FESystem includes
#include "Properties/ElemDataCard.h"
#include "Properties/MaterialPropertyNameEnums.h"


#ifndef FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_ID
#define FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_ID 6
#else
#error
#endif

#ifndef FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_NAME
#define FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_NAME "FORCED_CONVECTION_1D_ELEM_DATA_CARD"
#else
#error
#endif

#ifndef WALL_AREA_ENUM_ID
#define WALL_AREA_ENUM_ID  1
#else
#error
#endif

#ifndef WALL_AREA_ENUM_NAME
#define WALL_AREA_ENUM_NAME "WALL_AREA"
#else
#error
#endif

#ifndef FLUID_AREA_ENUM_ID
#define FLUID_AREA_ENUM_ID  2
#else
#error
#endif

#ifndef FLUID_AREA_ENUM_NAME
#define FLUID_AREA_ENUM_NAME "FLUID_AREA"
#else
#error
#endif


#ifndef CONTACT_PERIMETER_ENUM_ID
#define CONTACT_PERIMETER_ENUM_ID 3
#else
#error
#endif

#ifndef CONTACT_PERIMETER_ENUM_NAME
#define CONTACT_PERIMETER_ENUM_NAME  "CONTACT_PERIMETER"
#else
#error
#endif

#ifndef MASS_FLOW_RATE_ENUM_ID
#define MASS_FLOW_RATE_ENUM_ID 4
#else
#error
#endif

#ifndef MASS_FLOW_RATE_ENUM_NAME
#define MASS_FLOW_RATE_ENUM_NAME "MASS_FLOW_RATE"
#else
#error
#endif

#ifndef CONVECTION_COEFFICIENT_ENUM_ID
#define CONVECTION_COEFFICIENT_ENUM_ID 5
#else
#error
#endif

#ifndef CONVECTION_COEFFICIENT_ENUM_NAME
#define CONVECTION_COEFFICIENT_ENUM_NAME "CONVECTION_COEFFICIENT"
#else
#error
#endif


#ifndef FLUID_CONVECTIVE_FACTOR_ENUM_ID
#define FLUID_CONVECTIVE_FACTOR_ENUM_ID 1
#else
#error
#endif

#ifndef FLUID_CONVECTIVE_FACTOR_ENUM_NAME
#define FLUID_CONVECTIVE_FACTOR_ENUM_NAME "FLUID_CONVECTIVE_FACTOR"
#else
#error
#endif


#ifndef CONVECTIVE_EXCANGE_FACTOR_ENUM_ID
#define CONVECTIVE_EXCANGE_FACTOR_ENUM_ID 2
#else
#error
#endif

#ifndef CONVECTIVE_EXCANGE_FACTOR_ENUM_NAME
#define CONVECTIVE_EXCANGE_FACTOR_ENUM_NAME "CONVECTIVE_EXCANGE_FACTOR"
#else
#error
#endif


#ifndef WALL_CONDUCTIVITY_FACTOR_ENUM_ID
#define WALL_CONDUCTIVITY_FACTOR_ENUM_ID 3
#else
#error
#endif

#ifndef WALL_CONDUCTIVITY_FACTOR_ENUM_NAME
#define WALL_CONDUCTIVITY_FACTOR_ENUM_NAME "WALL_CONDUCTIVITY"
#else
#error
#endif


#ifndef FLUID_CONDUCTIVITY_FACTOR_ENUM_ID
#define FLUID_CONDUCTIVITY_FACTOR_ENUM_ID 4
#else
#error
#endif

#ifndef FLUID_CONDUCTIVITY_FACTOR_ENUM_NAME
#define FLUID_CONDUCTIVITY_FACTOR_ENUM_NAME "FLUID_CONDUCTIVITY"
#else
#error
#endif

#ifndef FLUID_CAPACITANCE_FACTOR_ENUM_ID
#define FLUID_CAPACITANCE_FACTOR_ENUM_ID 5
#else
#error
#endif

#ifndef FLUID_CAPACITANCE_FACTOR_ENUM_NAME
#define FLUID_CAPACITANCE_FACTOR_ENUM_NAME "FLUID_CAPACITANCE"
#else
#error
#endif


#ifndef WALL_CAPACITANCE_FACTOR_ENUM_ID
#define WALL_CAPACITANCE_FACTOR_ENUM_ID 6
#else
#error
#endif

#ifndef WALL_CAPACITANCE_FACTOR_ENUM_NAME
#define WALL_CAPACITANCE_FACTOR_ENUM_NAME "WALL_CAPACITANCE"
#else
#error
#endif

DeclareEnumName(FORCED_CONVECTION_1D_ELEM_DATA_CARD,
                ElemDataCardEnum,
                FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_ID,
                FORCED_CONVECTION_1D_ELEM_DATA_CARD_ENUM_NAME);


DeclareEnumClass(ForcedConvection_1D_ElemDataEnum);

DeclareEnumClass(ForcedConvection_1D_ElemDataCardFactorsEnum);



DeclareEnumName(WALL_AREA, ForcedConvection_1D_ElemDataEnum, WALL_AREA_ENUM_ID, WALL_AREA_ENUM_NAME);


DeclareEnumName(FLUID_AREA, ForcedConvection_1D_ElemDataEnum, 
                FLUID_AREA_ENUM_ID , FLUID_AREA_ENUM_NAME);


DeclareEnumName(CONTACT_PERIMETER, ForcedConvection_1D_ElemDataEnum,
                CONTACT_PERIMETER_ENUM_ID, CONTACT_PERIMETER_ENUM_NAME);


DeclareEnumName(MASS_FLOW_RATE, ForcedConvection_1D_ElemDataEnum,
                MASS_FLOW_RATE_ENUM_ID, MASS_FLOW_RATE_ENUM_NAME);


DeclareEnumName(CONVECTION_COEFFICIENT, ForcedConvection_1D_ElemDataEnum,
                CONVECTION_COEFFICIENT_ENUM_ID, CONVECTION_COEFFICIENT_ENUM_NAME);


DeclareEnumName(FLUID_CAPACITANCE_FACTOR, ForcedConvection_1D_ElemDataCardFactorsEnum,
                FLUID_CAPACITANCE_FACTOR_ENUM_ID, FLUID_CAPACITANCE_FACTOR_ENUM_NAME);


DeclareEnumName(WALL_CAPACITANCE_FACTOR, ForcedConvection_1D_ElemDataCardFactorsEnum,
                WALL_CAPACITANCE_FACTOR_ENUM_ID, WALL_CAPACITANCE_FACTOR_ENUM_NAME);


DeclareEnumName(FLUID_CONVECTIVE_FACTOR, ForcedConvection_1D_ElemDataCardFactorsEnum,
                FLUID_CONVECTIVE_FACTOR_ENUM_ID, FLUID_CONVECTIVE_FACTOR_ENUM_NAME);


DeclareEnumName(CONVECTIVE_EXCANGE_FACTOR, ForcedConvection_1D_ElemDataCardFactorsEnum,
                CONVECTIVE_EXCANGE_FACTOR_ENUM_ID, CONVECTIVE_EXCANGE_FACTOR_ENUM_NAME);


DeclareEnumName(WALL_CONDUCTIVITY_FACTOR, ForcedConvection_1D_ElemDataCardFactorsEnum,
                WALL_CONDUCTIVITY_FACTOR_ENUM_ID, WALL_CONDUCTIVITY_FACTOR_ENUM_NAME);


DeclareEnumName(FLUID_CONDUCTIVITY_FACTOR, ForcedConvection_1D_ElemDataCardFactorsEnum,
                FLUID_CONDUCTIVITY_FACTOR_ENUM_ID, FLUID_CONDUCTIVITY_FACTOR_ENUM_NAME);




class ForcedConvection_1D_ElemDataCard: public ElemDataCard
{
public:
  ForcedConvection_1D_ElemDataCard();
  
  virtual  ~ForcedConvection_1D_ElemDataCard();
    
  /// @returns the property card kind enumeration ID
  virtual unsigned int getPropertyCardKindEnumID() const;
  
  /// @returns the property card kind enumeration name
  virtual const std::string getPropertyCardKindEnumName() const;

  /// returns the wall material card ID for this card
  inline unsigned int getWallMaterialCardID() const;
  
  /// returns the fluid material card ID for this card
  inline unsigned int getFluidMaterialCardID() const;

  /// returns the wall material card for this card
  inline Property::PropertyCard<double>& getWallMaterialCard();
  
  /// returns the fluid material card for this card
  inline Property::PropertyCard<double>& getFluidMaterialCard();

  /// @returns true if either the element data card, or the material property card are
  /// dependent on the specified global parameter
  virtual inline bool checkElemAndMaterialCardGlobalParameterDependence
    (const unsigned int global_param_ID);
  
  /// @returns true if either the element data card, or the material property card are
  /// dependent on the specified local parameter
  virtual inline bool checkElemAndMaterialCardLocalParameterDependence
    (const unsigned int local_param_ID);
  
  
  /// initializes the card and the material card for the given local parameters
  virtual inline void reinitElemAndMaterialCardForLocalParameters
    (std::map<unsigned int, double> *value_map = NULL);
  
  /// initializes the card and the material card for the given global parameters
  virtual inline void reinitElemAndMaterialCardForGlobalParameters
    (std::map<unsigned int, double> &value_map);
  
  /// clears initialization of the card and the material card for local parameters
  virtual inline void clearElemAndMaterialCardLocalParameterInitialization();
  
  /// clears initialization of the card and the material card for local parameters
  virtual inline void 
    partialClearElemAndMaterialCardLocalParameterInitialization
    (std::vector<unsigned int >& param_ids);
  
  /// clears initialization of the card and the material card for global parameters
  virtual inline void clearElemAndMaterialCardGlobalParameterInitialization();

  virtual inline std::istream& readFromInputStream(std::istream& input);

  /// returns the value of the factor
  virtual void getFactor(double& factor, const unsigned int factor_enum_ID);
  
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
  
  virtual void 
    getPropertyValueFromMaterialCard(const unsigned int prop_enum_ID, double& value,
                                     const unsigned int layer = 0)
  {
    // params not used here
    (void) prop_enum_ID;
    (void) value;
    (void) layer;

    Assert(false, ExcInternalError());
  }
  
  /// returns value of the derivative of a property from the material card
  virtual void
    getPropertyValueDerivativeForGlobalParameterFromMaterialCard
    (const unsigned int prop_enum_ID, 
     const unsigned int global_param_ID,
     double& value,
     const unsigned int layer = 0)
  {
    // params not used here
    (void) prop_enum_ID;
    (void) global_param_ID;
    (void) value;
    (void) layer;
    
    Assert(false, ExcInternalError());
  }
  
  /// returns value of the derivative of a property from the material card
  virtual void
    getPropertyValueDerivativeForLocalParameterFromMaterialCard
    (const unsigned int prop_enum_ID, 
     const unsigned int local_param_ID,
     double& value,
     const unsigned int layer = 0)
  {
    // params not used here
    (void) prop_enum_ID;
    (void) local_param_ID;
    (void) value;
    (void) layer;
    
    Assert(false, ExcInternalError());
  }
  
  friend inline std::istream& operator>> (std::istream& input, ForcedConvection_1D_ElemDataCard& card);

protected:
    
  /// @returns the property name from enumeration ID of the property
  virtual std::string getPropertyEnumName(const unsigned int enum_ID) const;
  
  /// @returns the property ID from enumeration name of the property
  virtual const unsigned int getPropertyEnumID(const std::string& enum_name) const; 
  
  /// material card ID for this card
  unsigned int fluid_material_ID;
  
  /// material card ID for this card
  unsigned int wall_material_ID;

  /// this is a pointer to the material property card for this element, and is stored for 
  /// convenience
  Property::PropertyCard<double> *fluid_material_card;
  
  /// this is a pointer to the material property card for this element, and is stored for 
  /// convenience
  Property::PropertyCard<double> *wall_material_card;
};



inline
unsigned int
ForcedConvection_1D_ElemDataCard::getPropertyCardKindEnumID() const
{
  return FORCED_CONVECTION_1D_ELEM_DATA_CARD::num();
}





inline
const std::string
ForcedConvection_1D_ElemDataCard::getPropertyCardKindEnumName() const
{
  return FORCED_CONVECTION_1D_ELEM_DATA_CARD::name();
}




inline
std::string
ForcedConvection_1D_ElemDataCard::getPropertyEnumName(const unsigned int enum_ID) const
{
  return ForcedConvection_1D_ElemDataEnum::enumName(enum_ID);
}



inline
const unsigned int
ForcedConvection_1D_ElemDataCard::getPropertyEnumID(const std::string& enum_name) const
{
  return ForcedConvection_1D_ElemDataEnum::enumID(enum_name);
}



inline
bool
ForcedConvection_1D_ElemDataCard::checkElemAndMaterialCardGlobalParameterDependence
(const unsigned int global_param_ID)
{
  if (this->checkGlobalParameterDependence(global_param_ID) &&
      this->getFluidMaterialCard().checkGlobalParameterDependence(global_param_ID) &&
      this->getWallMaterialCard().checkGlobalParameterDependence(global_param_ID))
    return true;
  else
    return false;
}




inline
bool
ForcedConvection_1D_ElemDataCard::checkElemAndMaterialCardLocalParameterDependence
(const unsigned int local_param_ID)
{
  if (this->checkLocalParameterDependence(local_param_ID) &&
      this->getFluidMaterialCard().checkLocalParameterDependence(local_param_ID) &&
      this->getWallMaterialCard().checkLocalParameterDependence(local_param_ID))
    return true;
  else
    return false;
}



inline 
void
ForcedConvection_1D_ElemDataCard::reinitElemAndMaterialCardForLocalParameters
(std::map<unsigned int, double> *value_map)
{
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  
  this->reinitLocalParameters(value_map);
  fluid_property.reinitLocalParameters(value_map);
  wall_property.reinitLocalParameters(value_map);
}




inline
void
ForcedConvection_1D_ElemDataCard::reinitElemAndMaterialCardForGlobalParameters
(std::map<unsigned int, double> &value_map)
{
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  
  this->reinitGlobalParameters(value_map);
  fluid_property.reinitGlobalParameters(value_map);
  wall_property.reinitGlobalParameters(value_map);
}






inline
void
ForcedConvection_1D_ElemDataCard::clearElemAndMaterialCardLocalParameterInitialization()
{
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  
  this->clearLocalParameterInitialization();
  fluid_property.clearLocalParameterInitialization();
  wall_property.clearLocalParameterInitialization();
}



inline
void
ForcedConvection_1D_ElemDataCard::
partialClearElemAndMaterialCardLocalParameterInitialization(std::vector<unsigned int>& param_ids)
{
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  
  this->partialClearLocalParameterInitialization(param_ids);
  fluid_property.partialClearLocalParameterInitialization(param_ids);
  wall_property.partialClearLocalParameterInitialization(param_ids);
}



inline
void
ForcedConvection_1D_ElemDataCard::clearElemAndMaterialCardGlobalParameterInitialization()
{
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  
  this->clearGlobalParameterInitialization();
  fluid_property.clearGlobalParameterInitialization();
  wall_property.clearGlobalParameterInitialization();
}




inline unsigned int 
ForcedConvection_1D_ElemDataCard::getFluidMaterialCardID() const
{
  Assert(this->fluid_material_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(this->fluid_material_ID));
  
  return this->fluid_material_ID;
}


inline unsigned int 
ForcedConvection_1D_ElemDataCard::getWallMaterialCardID() const
{
  Assert(this->wall_material_ID != FESystemNumbers::InvalidID,
         FESystemExceptions::ExcInvalidID(this->wall_material_ID));
  
  return this->wall_material_ID;
}



inline
Property::PropertyCard<double>&
ForcedConvection_1D_ElemDataCard::getFluidMaterialCard()
{
  if (this->fluid_material_card == NULL)
    this->fluid_material_card =
      &(this->property_database->getMaterialPropertyCardFromID<double>(this->fluid_material_ID));
  
  
  return *(this->fluid_material_card);
}



inline
Property::PropertyCard<double>&
ForcedConvection_1D_ElemDataCard::getWallMaterialCard()
{
  if (this->wall_material_card == NULL)
    this->wall_material_card =
      &(this->property_database->getMaterialPropertyCardFromID<double>(this->wall_material_ID));
  
  
  return *(this->wall_material_card);
}



inline 
void
ForcedConvection_1D_ElemDataCard::getFactor(double& factor,
                                            const unsigned int factor_enum_ID)
{
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  
  factor = 0.0;
  static double val1, val2, val3;
  val1 = 0.0; val2 = 0.0, val3 = 0.0;
  
  switch (factor_enum_ID)
    {
    case FLUID_CONVECTIVE_FACTOR_ENUM_ID:
      {
        // mass_flow * cf
        fluid_property.getPropertyValue(SPECIFIC_HEAT::num(), val1);
        this->getPropertyValue(MASS_FLOW_RATE::num(), val2);
        factor = val1*val2;
      }
      break;
      
    case CONVECTIVE_EXCANGE_FACTOR_ENUM_ID:
      {
        // h p
        this->getPropertyValue(CONTACT_PERIMETER::num(), val1);
        this->getPropertyValue(CONVECTION_COEFFICIENT::num(), val2);
        factor = val1*val2;
      }
      break;

    case WALL_CONDUCTIVITY_FACTOR_ENUM_ID:
      {
        // kw Aw
        wall_property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(WALL_AREA::num(), val2);
        factor = val1*val2;
      }
      break;

    case FLUID_CONDUCTIVITY_FACTOR_ENUM_ID:
      {
        // kf Af
        fluid_property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        this->getPropertyValue(FLUID_AREA::num(), val2);
        factor = val1*val2;
      }
      break;

    case FLUID_CAPACITANCE_FACTOR_ENUM_ID:
      {
        // rhof Af Cf
        fluid_property.getPropertyValue(DENSITY::num(), val1);
        fluid_property.getPropertyValue(SPECIFIC_HEAT::num(), val3);
        this->getPropertyValue(FLUID_AREA::num(), val2);
        factor = val1 * val2 * val3;
      }
      break;
      
    case WALL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        // rhow Aw Cw
        wall_property.getPropertyValue(DENSITY::num(), val1);
        wall_property.getPropertyValue(SPECIFIC_HEAT::num(), val3);
        this->getPropertyValue(WALL_AREA::num(), val2);
        factor = val1 * val2 * val3;
      }
      break;
      
default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ForcedConvection_1D_ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}




inline 
void
ForcedConvection_1D_ElemDataCard::getFactorSensitivityForGlobalParameter
(double& factor,
 const unsigned int factor_enum_ID,
 const unsigned int global_param_ID)
{
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  
  factor = 0.0;
  static double val1, val2, val1_sens, val2_sens, val3, val3_sens;
  val1 = 0.0; val2 = 0.0; val1_sens = 0.0; val2_sens = 0.0; val3 = 0.0; val3_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    case FLUID_CONVECTIVE_FACTOR_ENUM_ID:
      {
        // mass_flow * cf
        fluid_property.getPropertyValue(SPECIFIC_HEAT::num(), val1);
        fluid_property.getPropertyDerivativeForGlobalParameter(SPECIFIC_HEAT::num(),
                                                               global_param_ID,
                                                               val1_sens);

        this->getPropertyValue(MASS_FLOW_RATE::num(), val2);
        this->getPropertyDerivativeForGlobalParameter(MASS_FLOW_RATE::num(),
                                                      global_param_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case CONVECTIVE_EXCANGE_FACTOR_ENUM_ID:
      {
        // h p
        this->getPropertyValue(CONTACT_PERIMETER::num(), val1);
        this->getPropertyDerivativeForGlobalParameter(CONTACT_PERIMETER::num(), 
                                                      global_param_ID,
                                                      val1_sens);
        
        this->getPropertyValue(CONVECTION_COEFFICIENT::num(), val2);
        this->getPropertyDerivativeForGlobalParameter(CONVECTION_COEFFICIENT::num(),
                                                      global_param_ID,
                                                      val2_sens);

        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;

    case WALL_CONDUCTIVITY_FACTOR_ENUM_ID:
      {
        // kw Aw
        wall_property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        wall_property.getPropertyDerivativeForGlobalParameter(THERMAL_CONDUCTIVITY::num(),
                                                              global_param_ID,
                                                              val1_sens);

        this->getPropertyValue(WALL_AREA::num(), val2);
        this->getPropertyDerivativeForGlobalParameter(WALL_AREA::num(), 
                                                      global_param_ID,
                                                      val2_sens);

        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;

    case FLUID_CONDUCTIVITY_FACTOR_ENUM_ID:
      {
        // kf Af
        fluid_property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        fluid_property.getPropertyDerivativeForGlobalParameter(THERMAL_CONDUCTIVITY::num(),
                                                               global_param_ID,
                                                               val1_sens);

        this->getPropertyValue(FLUID_AREA::num(), val2);
        this->getPropertyDerivativeForGlobalParameter(FLUID_AREA::num(),
                                                      global_param_ID,
                                                      val2_sens);

        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
            
    case FLUID_CAPACITANCE_FACTOR_ENUM_ID:
      {
        // rhof Af Cf
        fluid_property.getPropertyValue(DENSITY::num(), val1);
        fluid_property.getPropertyDerivativeForGlobalParameter(DENSITY::num(),
                                                               global_param_ID,
                                                               val1_sens);
        fluid_property.getPropertyValue(SPECIFIC_HEAT::num(), val3);
        fluid_property.getPropertyDerivativeForGlobalParameter(SPECIFIC_HEAT::num(),
                                                               global_param_ID,
                                                               val3_sens);
        this->getPropertyValue(FLUID_AREA::num(), val2);
        this->getPropertyDerivativeForGlobalParameter(FLUID_AREA::num(),
                                                      global_param_ID,
                                                      val2_sens);
        factor = val1_sens * val2 * val3 + val1 * val2_sens * val3 + val1 * val2 * val3_sens;
      }
      break;
      
    case WALL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        // rhow Aw Cw
        wall_property.getPropertyValue(DENSITY::num(), val1);
        wall_property.getPropertyDerivativeForGlobalParameter(DENSITY::num(),
                                                               global_param_ID,
                                                               val1_sens);
        wall_property.getPropertyValue(SPECIFIC_HEAT::num(), val3);
        wall_property.getPropertyDerivativeForGlobalParameter(SPECIFIC_HEAT::num(),
                                                               global_param_ID,
                                                               val3_sens);
        this->getPropertyValue(WALL_AREA::num(), val2);
        this->getPropertyDerivativeForGlobalParameter(WALL_AREA::num(),
                                                      global_param_ID,
                                                      val2_sens);
        factor = val1_sens * val2 * val3 + val1 * val2_sens * val3 + val1 * val2 * val3_sens;
      }
      break;

    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ForcedConvection_1D_ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}





inline 
void
ForcedConvection_1D_ElemDataCard::getFactorSensitivityForLocalParameter
(double& factor,
 const unsigned int factor_enum_ID,
 const unsigned int param_enum_ID)
{
  Property::PropertyCard<double>& wall_property = this->getWallMaterialCard();
  Property::PropertyCard<double>& fluid_property = this->getFluidMaterialCard();
  
  factor = 0.0;
  static double val1, val2, val1_sens, val2_sens, val3, val3_sens;
  val1 = 0.0; val2 = 0.0; val1_sens = 0.0; val2_sens = 0.0; val3 = 0.0; val3_sens = 0.0;
  
  switch (factor_enum_ID)
    {
    case FLUID_CONVECTIVE_FACTOR_ENUM_ID:
      {
        // mass_flow * cf
        fluid_property.getPropertyValue(SPECIFIC_HEAT::num(), val1);
        fluid_property.getPropertyDerivativeForLocalParameter(SPECIFIC_HEAT::num(),
                                                               param_enum_ID,
                                                               val1_sens);
        
        this->getPropertyValue(MASS_FLOW_RATE::num(), val2);
        this->getPropertyDerivativeForLocalParameter(MASS_FLOW_RATE::num(),
                                                      param_enum_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case CONVECTIVE_EXCANGE_FACTOR_ENUM_ID:
      {
        // h p
        this->getPropertyValue(CONTACT_PERIMETER::num(), val1);
        this->getPropertyDerivativeForLocalParameter(CONTACT_PERIMETER::num(), 
                                                      param_enum_ID,
                                                      val1_sens);
        
        this->getPropertyValue(CONVECTION_COEFFICIENT::num(), val2);
        this->getPropertyDerivativeForLocalParameter(CONVECTION_COEFFICIENT::num(),
                                                      param_enum_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case WALL_CONDUCTIVITY_FACTOR_ENUM_ID:
      {
        // kw Aw
        wall_property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        wall_property.getPropertyDerivativeForLocalParameter(THERMAL_CONDUCTIVITY::num(),
                                                              param_enum_ID,
                                                              val1_sens);
        
        this->getPropertyValue(WALL_AREA::num(), val2);
        this->getPropertyDerivativeForLocalParameter(WALL_AREA::num(), 
                                                      param_enum_ID,
                                                      val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case FLUID_CONDUCTIVITY_FACTOR_ENUM_ID:
      {
        // kf Af
        fluid_property.getPropertyValue(THERMAL_CONDUCTIVITY::num(), val1);
        fluid_property.getPropertyDerivativeForLocalParameter(THERMAL_CONDUCTIVITY::num(),
                                                               param_enum_ID,
                                                               val1_sens);
        
        this->getPropertyValue(FLUID_AREA::num(), val2);
        this->getPropertyDerivativeForLocalParameter(FLUID_AREA::num(),
                                                     param_enum_ID,
                                                     val2_sens);
        
        factor = val1 * val2_sens + val1_sens * val2;
      }
      break;
      
    case FLUID_CAPACITANCE_FACTOR_ENUM_ID:
      {
        // rhof Af Cf
        fluid_property.getPropertyValue(DENSITY::num(), val1);
        fluid_property.getPropertyDerivativeForLocalParameter(DENSITY::num(),
                                                               param_enum_ID,
                                                               val1_sens);
        fluid_property.getPropertyValue(SPECIFIC_HEAT::num(), val3);
        fluid_property.getPropertyDerivativeForLocalParameter(SPECIFIC_HEAT::num(),
                                                               param_enum_ID,
                                                               val3_sens);
        this->getPropertyValue(FLUID_AREA::num(), val2);
        this->getPropertyDerivativeForLocalParameter(FLUID_AREA::num(),
                                                      param_enum_ID,
                                                      val2_sens);
        factor = val1_sens * val2 * val3 + val1 * val2_sens * val3 + val1 * val2 * val3_sens;
      }
      break;
      
    case WALL_CAPACITANCE_FACTOR_ENUM_ID:
      {
        // rhow Aw Cw
        wall_property.getPropertyValue(DENSITY::num(), val1);
        wall_property.getPropertyDerivativeForLocalParameter(DENSITY::num(),
                                                              param_enum_ID,
                                                              val1_sens);
        wall_property.getPropertyValue(SPECIFIC_HEAT::num(), val3);
        wall_property.getPropertyDerivativeForLocalParameter(SPECIFIC_HEAT::num(),
                                                              param_enum_ID,
                                                              val3_sens);
        this->getPropertyValue(WALL_AREA::num(), val2);
        this->getPropertyDerivativeForLocalParameter(WALL_AREA::num(),
                                                      param_enum_ID,
                                                      val2_sens);
        factor = val1_sens * val2 * val3 + val1 * val2_sens * val3 + val1 * val2 * val3_sens;
      }
      break;

    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (ForcedConvection_1D_ElemDataCardFactorsEnum::enumName(factor_enum_ID)));
    }
}



inline
std::istream& 
ForcedConvection_1D_ElemDataCard::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  // read the beginning data of the card
  tag = this->getPropertyCardKindEnumName();
  FESystemIOUtility::readFromInput(input, tag);
  
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->cardID);
  
  FESystemIOUtility::readFromInput(input, "FLUID_MATERIAL_PROPERTY_CARD_ID", this->fluid_material_ID);
  FESystemIOUtility::readFromInput(input, "WALL_MATERIAL_PROPERTY_CARD_ID", this->wall_material_ID);
  
  // now read in the property base values
  {
    unsigned int n_prop = 0, prop_ID = 0;
    double value = 0.0;
    bool insert_success = false;
    FESystemIOUtility::readFromInput(input, "N_PROPERTIES", n_prop);
    
    Property::PropertyCard<double>::PropertyMap::iterator prop_it;
    
    for (unsigned int i=0; i<n_prop; i++)
      {
      tag.clear(); value = 0.0;
      input >> tag;
      prop_ID = this->getPropertyEnumID(tag);
      input >> value;
      // now insert the value in the map
      insert_success = this->property_reference_values.insert
        (Property::PropertyCard<double>::PropertyMap::value_type
         (prop_ID, value)).second;
      
      Assert(insert_success, 
             Property::PropertyCardBase::ExcDuplicatePropertySpecified(tag));
      
      insert_success = this->property_enumID_to_internal_ID_map.insert
        (Property::PropertyCardBase::IDMap::value_type(prop_ID, i)).second;
      
      // if the previous assert was true, this one should also be true.
      Assert(insert_success, 
             Property::PropertyCardBase::ExcDuplicatePropertySpecified(tag));
      
      }
  }
  
  // now read the parameters
  this->readParametersFromInputStream(input);
  
  // now read in the property-parameter function ID table
  {
    // first initialize the table
    unsigned int n_prop = this->getNProperties(),
    n_param = this->getNParameters();
    
    this->parameter_function_ID_table->reinit(n_prop, n_param);
    this->parameter_function_ID_table->reset_values(FESystemNumbers::InvalidID);
    
    TableIndices<3> indices(n_prop, n_param, 2);
    this->parameter_function_value_table->reinit(indices);
    this->parameter_function_value_table->reset_values(0.0);
    
    unsigned int n_func_IDs, prop_ID = 0, param_ID = 0, func_ID = 0,
      prop_internal_ID = 0, param_internal_ID = 0;
    
    FESystemIOUtility::readFromInput(input, "N_FUNCTION_IDS", n_func_IDs);
    
    for (unsigned int i=0; i<n_func_IDs; i++)
      {
      tag.clear(); 
      input >> tag;
      prop_ID = this->getPropertyEnumID(tag);
      input >> param_ID;
      input >> func_ID;
      
      // get the internal IDS
      prop_internal_ID = this->getPropertyInternalID(prop_ID);
      param_internal_ID = this->getParameterInternalID(param_ID);
      
      // and, insert the value in the table
      // first make sure that this property hasn't already been set. 
      Assert((*this->parameter_function_ID_table)(prop_internal_ID, param_internal_ID)  ==
             FESystemNumbers::InvalidID,
             FESystemExceptions::ExcDuplicateID("Property Function",prop_internal_ID));
      
      // and, insert the value in the table
      (*this->parameter_function_ID_table)(prop_internal_ID, param_internal_ID) = func_ID;
      }
  }
  
  // the card should end with the END tag
  tag = this->getPropertyCardKindEnumName();
  FESystemIOUtility::readFromInput(input, tag);
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
  
}




inline
ForcedConvection_1D_ElemDataCard::ForcedConvection_1D_ElemDataCard():
ElemDataCard()
{
  
}

inline
ForcedConvection_1D_ElemDataCard::~ForcedConvection_1D_ElemDataCard()
{
  
}


inline
std::istream& operator>> (std::istream& input, ForcedConvection_1D_ElemDataCard& card)
{
  return card.readFromInputStream(input);
}



inline
void
ForcedConvection_1D_ElemDataCard::getFactor(DenseMatrix<double>& factor,
                                    const unsigned int factor_enum_ID)
{
  // params not used for now
  (void) factor;
  (void) factor_enum_ID;
  
  Assert(false, ExcPureFunctionCalled());
}



inline
void
ForcedConvection_1D_ElemDataCard::getFactorSensitivityForGlobalParameter
(DenseMatrix<double>& factor,
 const unsigned int factor_enum_ID,
 const unsigned int global_param_ID)
{
  // params not used here
  (void) factor;
  (void) factor_enum_ID;
  (void) global_param_ID;
  
  Assert(false, ExcPureFunctionCalled());
}


inline
void
ForcedConvection_1D_ElemDataCard::getFactorSensitivityForLocalParameter
(DenseMatrix<double>& factor,
 const unsigned int factor_enum_ID,
 const unsigned int param_enum_ID)
{
  // params not used here
  (void) factor;
  (void) factor_enum_ID;
  (void) param_enum_ID;
   
  Assert(false, ExcPureFunctionCalled());
}


#endif // __fesystem_forced_convection_1d_elem_data_card_h__
