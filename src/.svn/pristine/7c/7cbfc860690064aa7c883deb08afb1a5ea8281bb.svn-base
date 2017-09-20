// $Id: ShapeParameter.h,v 1.3 2006-09-05 20:41:41 manav Exp $

#ifndef __fesystem_shape_parameter_h__
#define __fesystem_shape_parameter_h__

// C++ includes
#include <iostream>
#include <string>

// FESystem includes
#include "DesignData/DesignParameter.h"

#ifndef SHAPE_PARAMETER_ENUM_ID 
#define SHAPE_PARAMETER_ENUM_ID 2
#else
#error
#endif


#ifndef SHAPE_PARAMETER_ENUM_NAME
#define SHAPE_PARAMETER_ENUM_NAME "SHAPE_PARAMETER"
#else
#error
#endif



namespace DesignData
{

DeclareEnumName(SHAPE_PARAMETER, DesignParameterTypeEnum, 
                SHAPE_PARAMETER_ENUM_ID,
                SHAPE_PARAMETER_ENUM_NAME);




class ShapeParameter: public DesignParameter
{
public: 
  ShapeParameter();
  
  virtual ~ShapeParameter();
  
  inline virtual const unsigned int getParameterTypeEnumID() const;
  
  inline virtual const std::string getParameterTypeEnumName() const;  
  
  inline const unsigned int getPerturbedMeshID(const unsigned int discipline) const;
  
  friend std::istream& operator>>(std::istream& input, 
                                  DesignData::ShapeParameter& param);
  
  std::istream& readFromInputStream(std::istream& input);

  
protected: 
    
    typedef std::map<unsigned int, unsigned int> IDMap;
    
  ShapeParameter::IDMap perturbed_mesh_ID_map;
};
}



inline
const unsigned int  
DesignData::ShapeParameter::getParameterTypeEnumID() const
{
  return DesignData::SHAPE_PARAMETER::num();
}


inline
const std::string
DesignData::ShapeParameter::getParameterTypeEnumName() const
{
  return DesignData::SHAPE_PARAMETER::name();
}



inline 
const unsigned int
DesignData::ShapeParameter::getPerturbedMeshID(const unsigned int discipline) const
{
  DesignData::ShapeParameter::IDMap::const_iterator it, end;
  it = this->perturbed_mesh_ID_map.find(discipline);
  end = this->perturbed_mesh_ID_map.end();
  
  Assert(it != end, 
         FESystemExceptions::ExcIDDoesNotExist("Discipline", discipline));
  
  return it->second;
}


#endif // __fesystem_shape_parameter_h__

