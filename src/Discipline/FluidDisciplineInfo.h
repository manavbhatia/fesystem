// $Id: FluidDisciplineInfo.h,v 1.1.6.1 2008-02-25 04:21:42 manav Exp $

#ifndef __fesystem_structural_discipline_info_h__
#define __fesystem_structural_discipline_info_h__

// FESystem includes
#include "Discipline/DisciplineInfo.h"
#include "Discipline/StructuralAnalysis.h"


#ifndef STRUCTURAL_DISCIPLINE_INFO_ENUM_ID
#define STRUCTURAL_DISCIPLINE_INFO_ENUM_ID 2
#else
#error
#endif

#ifndef STRUCTURAL_DISCIPLINE_INFO_ENUM_NAME
#define STRUCTURAL_DISCIPLINE_INFO_ENUM_NAME "STRUCTURAL_DISCIPLINE_INFO"
#else
#error
#endif



namespace Discipline
{
  DeclareEnumName(STRUCTURAL_DISCIPLINE_INFO, Discipline::DisciplineInfoTypeEnum,
                  STRUCTURAL_DISCIPLINE_INFO_ENUM_ID, STRUCTURAL_DISCIPLINE_INFO_ENUM_NAME);
  
  
  class StructuralDisciplineInfo: public DisciplineInfo
    {
public:
      
      StructuralDisciplineInfo();
      virtual ~StructuralDisciplineInfo();
      
      virtual inline unsigned int getDisciplineEnumID() const;
      virtual inline const std::string getDisciplineEnumName() const;
      virtual std::istream& readFromInputStream(std::istream& input);
      
      friend std::istream& operator >> (std::istream& input, 
                                        Discipline::StructuralDisciplineInfo& info);
      
protected:
        
    };
}




inline 
unsigned int
Discipline::StructuralDisciplineInfo::getDisciplineEnumID() const
{
  return Discipline::STRUCTURAL_DISCIPLINE::num();
}





inline 
const std::string
Discipline::StructuralDisciplineInfo::getDisciplineEnumName() const
{
  return Discipline::STRUCTURAL_DISCIPLINE::name();
}


#ifndef THERMOELASTICITY_DISCIPLINE_ENUM_ID
#define THERMOELASTICITY_DISCIPLINE_ENUM_ID 3
#else
#error
#endif


#ifndef THERMOELASTICITY_DISCIPLINE_ENUM_NAME
#define THERMOELASTICITY_DISCIPLINE_ENUM_NAME "THERMOELASTICITY_DISCIPLINE"
#else
#error
#endif



namespace Discipline
{
  DeclareEnumName(THERMOELASTICITY_DISCIPLINE, AnalysisDisciplineEnum,
                  THERMOELASTICITY_DISCIPLINE_ENUM_ID, 
                  THERMOELASTICITY_DISCIPLINE_ENUM_NAME);
}  


#endif //  __fesystem_structural_discipline_info_h__
