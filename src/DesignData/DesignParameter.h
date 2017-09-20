// $Id: DesignParameter.h,v 1.4 2007-01-15 18:59:43 manav Exp $

#ifndef __fesystem_design_parameter_h__
#define __fesystem_design_parameter_h__

// C++ includes
#include <iostream>
#include <string>

// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "Utilities/NameEnumHandler.h"
#include "FESystem/FESystemExceptions.h"



#ifndef ANALYTIC_SENSITIVITY_ENUM_ID
#define ANALYTIC_SENSITIVITY_ENUM_ID 1
#else
#error
#endif

#ifndef ANALYTIC_SENSITIVITY_ENUM_NAME
#define ANALYTIC_SENSITIVITY_ENUM_NAME "ANALYTIC_SENSITIVITY"
#else
#error
#endif


#ifndef EULER_FD_SENSITIVITY_ENUM_ID
#define EULER_FD_SENSITIVITY_ENUM_ID 2
#else
#error
#endif

#ifndef EULER_FD_SENSITIVITY_ENUM_NAME
#define EULER_FD_SENSITIVITY_ENUM_NAME "EULER_FD_SENSITIVITY"
#else
#error
#endif


namespace DesignData
{
  
  DeclareEnumClass(DesignParameterTypeEnum);
  
  
  DeclareEnumClass(SensitivityMethodEnum);
  
  
  DeclareEnumName(ANALYTIC_SENSITIVITY, SensitivityMethodEnum, 
                  ANALYTIC_SENSITIVITY_ENUM_ID,
                  ANALYTIC_SENSITIVITY_ENUM_NAME);
  
  
  DeclareEnumName(EULER_FD_SENSITIVITY, SensitivityMethodEnum, 
                  EULER_FD_SENSITIVITY_ENUM_ID,
                  EULER_FD_SENSITIVITY_ENUM_NAME);
  
  
  /// this class forms the base class for all design parameters.
  class DesignParameter
    {
public:
      DesignParameter();
      
      virtual ~DesignParameter();
      
      inline const unsigned int getID() const;
      
      inline const std::string& getName() const;
      
      virtual const unsigned int getParameterTypeEnumID() const = 0;
      
      virtual const std::string getParameterTypeEnumName() const = 0;
      
      inline const unsigned int getSensitivityMethodEnumID() const;
      
      inline const std::string getSensitivityMethodEnumName() const;
      
      inline const double getValue() const;
      
      inline const double getPerturbationStepSize() const;
      
protected:
        
        virtual std::istream& readFromInputStream(std::istream& input) = 0;  
      
      unsigned int ID;
      
      std::string name;
      
      double value;
      
      unsigned int sensitivity_method;
      
      double finite_difference_step_size;
      
      DeclException0(ExcFESystemNonFDSensitivity);
    };
  
}



inline
const unsigned int 
DesignData::DesignParameter::getID() const
{
  // make sure that the ID is a valid number before returning it
  Assert(this->ID != FESystemNumbers::InvalidID, 
         FESystemExceptions::ExcInvalidID(this->ID));
  
  return this->ID;
}


inline
const std::string& 
DesignData::DesignParameter::getName() const
{return this->name;}


inline
const unsigned int
DesignData::DesignParameter::getSensitivityMethodEnumID() const
{return this->sensitivity_method;}


inline
const std::string 
DesignData::DesignParameter::getSensitivityMethodEnumName() const
{return DesignData::SensitivityMethodEnum::enumName(this->sensitivity_method);}


const double 
DesignData::DesignParameter::getValue() const
{
  return this->value;
}



const double 
DesignData::DesignParameter::getPerturbationStepSize() const
{
  // make sure that the sensitivity method is finite difference before 
  // returning the step size
  Assert(this->sensitivity_method == DesignData::EULER_FD_SENSITIVITY::num(),
         ExcFESystemNonFDSensitivity());
  
  return this->finite_difference_step_size;
}



#endif // __fesystem_design_parameter_h__

