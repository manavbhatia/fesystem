// $Id: PropertyCardParameter.h,v 1.3 2006-09-05 20:41:51 manav Exp $

#ifndef __fesystem_property_card_parameter_h__
#define __fesystem_property_card_parameter_h__

//C++ include
#include <iostream>
#include <string>

// FESystem includes
#include "Utilities/NameEnumHandler.h"
#include "Utilities/InputOutputUtility.h"
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"


namespace Property
{
  
  /// the enumeration type class for the parameters
  DeclareEnumClass(PropertyCardParameterType);
  
  /// now define the parameter enum types
#define LOCAL_PROPERTY_PARAMETER_ENUM_ID 1
#define LOCAL_PROPERTY_PARAMETER_ENUM_NAME "LOCAL_PROPERTY_PARAMETER"
  DeclareEnumName(LOCAL_PROPERTY_PARAMETER, 
                  PropertyCardParameterType,
                  LOCAL_PROPERTY_PARAMETER_ENUM_ID,
                  LOCAL_PROPERTY_PARAMETER_ENUM_NAME);
  
#define GLOBAL_PROPERTY_PARAMETER_ENUM_ID 2
#define GLOBAL_PROPERTY_PARAMETER_ENUM_NAME "GLOBAL_PROPERTY_PARAMETER"
  DeclareEnumName(GLOBAL_PROPERTY_PARAMETER, 
                  PropertyCardParameterType,
                  GLOBAL_PROPERTY_PARAMETER_ENUM_ID,
                  GLOBAL_PROPERTY_PARAMETER_ENUM_NAME);  
  
  /// declare the local parameter types
  DeclareEnumClass(LocalParameterType);
  
#define TEMPERATURE_ENUM_ID 1
#define TEMPERATURE_ENUM_NAME "TEMPERATURE"
  DeclareEnumName(TEMPERATURE, 
                  LocalParameterType,
                  TEMPERATURE_ENUM_ID,
                  TEMPERATURE_ENUM_NAME);  
  
  
  /// a parameter type definition for storing the parameters defined in the input
  class ParameterBase
    {
public:
      /// constructor
      ParameterBase();
      
      /// destructor
      virtual  ~ParameterBase();
      
      /// @returns the ID of the parameter
      inline unsigned int getID() const;
      
      /// @returns the name of the parameter
      inline const std::string& getName() const;
      
      /// @returns the reference value of the parameter
      inline double getReferenceValue() const;
      
      /// @returs the current value of the parameter
      inline double getCurrentValue() const;
      
      /// @returns a writable reference to the current value of the parameter
      inline void setCurrentValue(const double value);
      
      /// @returns enum ID of the type of parameter
      virtual unsigned int getParameterTypeEnumID() const = 0;
      
      /// @returns enum name of the type of parameter
      virtual const std::string getParameterTypeEnumName() const = 0;
      
protected:
        
        /// reads from the input stream
        virtual std::istream& readFromInputStream(std::istream& input) = 0;
      
      /// the internal ID of this object
      unsigned int ID;
      
      /// name of the parameter
      std::string name;
      
      /// the reference value of the parameter at which properties are defined 
      /// in the class
      double ref_value;
      
      /// current value of the parameter
      double current_value;
    };
  
  
  /// a data structure for local parameters
  class LocalParameter: public Property::ParameterBase
    {
public:
      /// constructor
      LocalParameter();
      
      /// destructor
      virtual ~LocalParameter();
      
      /// @returns enum ID of the type of parameter
      inline unsigned int getLocalParameterEnumID() const;
      
      /// @returns enum name of the type of parameter
      inline const std::string getLocalParameterEnumName() const;
      
      /// @returns enum ID of the type of parameter
      inline unsigned int getParameterTypeEnumID() const;
      
      /// @returns enum name of the type of parameter
      inline const std::string getParameterTypeEnumName() const;
      
      /// overloaded input operator
      friend std::istream& operator>>(std::istream& input, LocalParameter& param);
      
protected:
        
        /// reads from the input stream
        std::istream& readFromInputStream(std::istream& input);
      
      /// enumeration ID of the local parameter name
      unsigned int param_enum_ID;
    };
  
  
  /// a data structure for global parameters
  class GlobalParameter: public Property::ParameterBase
    {
public:
      /// constructor
      GlobalParameter();
      
      /// destructor
      virtual  ~GlobalParameter();
      
      /// @returns the ID of the parameter in the global database
      inline unsigned int getGlobalParameterID() const;
      
      /// @returns enum ID of the type of parameter
      inline unsigned int getParameterTypeEnumID() const;
      
      /// @returns enum name of the type of parameter
      inline const std::string getParameterTypeEnumName() const;
      
      /// overloaded input operator
      friend std::istream& operator>>(std::istream& input, GlobalParameter& param);
      
protected:
        
        /// reads from the input stream
        std::istream& readFromInputStream(std::istream& input);
      
      /// ID of the global parameter in the parameter database
      unsigned int global_parameter_ID;
    };
  
  
  
  inline 
    unsigned int  
    ParameterBase::getID() const
    {
      Assert(this->ID != FESystemNumbers::InvalidID,
             FESystemExceptions::ExcInvalidID(this->ID));
      
      return this->ID;
    }
  
  
  inline 
    const std::string& 
    ParameterBase::getName() const
    {
      return this->name;
    }
  
  
  inline 
    double 
    ParameterBase::getReferenceValue() const
    {
      return this->ref_value;
    }
  
  inline 
    double 
    ParameterBase::getCurrentValue() const
    {
      return this->current_value;
    }
  
  inline 
    void
    ParameterBase::setCurrentValue(const double value)
    {
      this->current_value = value;
    }
  
  
  
  inline 
    unsigned int
    LocalParameter::getLocalParameterEnumID() const
    {
      Assert(this->param_enum_ID != FESystemNumbers::InvalidID,
             FESystemExceptions::ExcInvalidID(this->param_enum_ID));
      
      return this->param_enum_ID;
    }
  
  
  inline 
    const std::string 
    LocalParameter::getLocalParameterEnumName() const
    {
      return LocalParameterType::enumName(this->param_enum_ID);
    }
  
  
  inline 
    unsigned int
    LocalParameter::getParameterTypeEnumID() const
    {
      return Property::LOCAL_PROPERTY_PARAMETER::num();
    }
  
  
  inline
    const std::string 
    LocalParameter::getParameterTypeEnumName() const
    {
      return Property::LOCAL_PROPERTY_PARAMETER::name();
    }
  
  
  
  inline 
    unsigned int
    GlobalParameter::getGlobalParameterID() const
    {
      Assert(this->global_parameter_ID != FESystemNumbers::InvalidID,
             FESystemExceptions::ExcInvalidID(this->global_parameter_ID));
      
      return this->global_parameter_ID;
    }
  
  
  inline
    unsigned int
    GlobalParameter::getParameterTypeEnumID() const
    {
      return Property::GLOBAL_PROPERTY_PARAMETER::num();
    }
  
  
  inline
    const std::string 
    GlobalParameter::getParameterTypeEnumName() const
    {
      return Property::GLOBAL_PROPERTY_PARAMETER::name();
    }
  
}

#endif // __fesystem_property_card_parameter_h__



