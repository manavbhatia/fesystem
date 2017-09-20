// $Id: ThermalDisciplineInfo.h,v 1.3 2006-09-05 20:09:54 manav Exp $

#ifndef __fesystem_thermal_discipline_info_h__
#define __fesystem_thermal_discipline_info_h__

// FESystem includes
#include "Discipline/DisciplineInfo.h"
#include "Discipline/ThermalAnalysis.h"


#ifndef THERMAL_DISCIPLINE_INFO_ENUM_ID
#define THERMAL_DISCIPLINE_INFO_ENUM_ID 1
#else
#error
#endif

#ifndef THERMAL_DISCIPLINE_INFO_ENUM_NAME
#define THERMAL_DISCIPLINE_INFO_ENUM_NAME "THERMAL_DISCIPLINE_INFO"
#else
#error
#endif



namespace Discipline
{
  
  DeclareEnumName(THERMAL_DISCIPLINE_INFO, Discipline::DisciplineInfoTypeEnum,
                  THERMAL_DISCIPLINE_INFO_ENUM_ID, THERMAL_DISCIPLINE_INFO_ENUM_NAME);
  
  class ThermalDisciplineInfo: public DisciplineInfo
  {
public:
    
    ThermalDisciplineInfo();
    virtual ~ThermalDisciplineInfo();
    
    virtual inline unsigned int getDisciplineEnumID() const;
    virtual inline const std::string getDisciplineEnumName() const;
    inline unsigned int getNRadiationCavities() const;
    inline const std::vector<unsigned int>& getRadiationCavityIDs() const;
    virtual std::istream& readFromInputStream(std::istream& input);
    
    friend std::istream& operator >> (std::istream& input, 
                                      Discipline::ThermalDisciplineInfo& info);
    
protected:
      
      std::vector<unsigned int> radiation_cavity_IDs;
  };
}




inline 
unsigned int
Discipline::ThermalDisciplineInfo::getDisciplineEnumID() const
{
  return Discipline::THERMAL_DISCIPLINE::num();
}





inline 
const std::string
Discipline::ThermalDisciplineInfo::getDisciplineEnumName() const
{
  return Discipline::THERMAL_DISCIPLINE::name();
}





inline 
unsigned int
Discipline::ThermalDisciplineInfo::getNRadiationCavities() const
{
  return this->radiation_cavity_IDs.size();
}




inline 
const std::vector<unsigned int>&
Discipline::ThermalDisciplineInfo::getRadiationCavityIDs() const
{
  return this->radiation_cavity_IDs;
}


#endif //  __fesystem_thermal_discipline_info_h__
