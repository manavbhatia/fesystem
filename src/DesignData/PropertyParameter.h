// $Id: PropertyParameter.h,v 1.3.6.1 2008-08-21 20:42:16 manav Exp $

#ifndef __fesystem_property_parameter_h__
#define __fesystem_property_parameter_h__

// C++ includes
#include <iostream>
#include <string>

// FESystem includes
#include "DesignData/DesignParameter.h"


#ifndef PROPERTY_PARAMETER_ENUM_ID
#define PROPERTY_PARAMETER_ENUM_ID 1
#else
#error
#endif

#ifndef PROPERTY_PARAMETER_ENUM_NAME
#define PROPERTY_PARAMETER_ENUM_NAME "PROPERTY_PARAMETER"
#else
#error
#endif

namespace DesignData
{
  
  DeclareEnumName(PROPERTY_PARAMETER, DesignData::DesignParameterTypeEnum, 
                  PROPERTY_PARAMETER_ENUM_ID,
                  PROPERTY_PARAMETER_ENUM_NAME);
  
  
  
  class PropertyParameter: public DesignParameter
    {
public: 
      PropertyParameter();
      
      virtual ~PropertyParameter();
      
      inline virtual const unsigned int getParameterTypeEnumID() const;
      
      inline virtual const std::string getParameterTypeEnumName() const;  
      
      friend std::istream& operator>>(std::istream& input, 
                                      DesignData::PropertyParameter& param)
      {
        param.readFromInputStream(input);
        return input;
      }
      
      
      std::istream& readFromInputStream(std::istream& input);  

protected: 
        
      
    };
}


inline
const unsigned int  
DesignData::PropertyParameter::getParameterTypeEnumID() const
{
  return DesignData::PROPERTY_PARAMETER::num();
}


inline
const std::string
DesignData::PropertyParameter::getParameterTypeEnumName() const
{
  return DesignData::PROPERTY_PARAMETER::name();
}




#endif // __fesystem_design_parameter_h__

