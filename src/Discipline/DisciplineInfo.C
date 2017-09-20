// $Id: DisciplineInfo.C,v 1.3 2006-09-05 20:09:54 manav Exp $


// FESystem includes
#include "Discipline/DisciplineInfo.h"
#include "FESystem/FESystemNumbers.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "Discipline/PistonTheoryInfo.h"


Discipline::DisciplineInfo::DisciplineInfo():
ID(FESystemNumbers::InvalidID),
mesh_ID(FESystemNumbers::InvalidID),
analysis_type(FESystemNumbers::InvalidID)
{
  
}




Discipline::DisciplineInfo::~DisciplineInfo()
{
  
}



std::auto_ptr<Discipline::DisciplineInfo> 
Discipline::createDisciplineInfo(const unsigned int info_enum_id)
{
  
  std::auto_ptr<Discipline::DisciplineInfo> info(NULL);
  
  switch(info_enum_id)
    {
    case THERMAL_DISCIPLINE_INFO_ENUM_ID:
      info.reset(new Discipline::ThermalDisciplineInfo());
      break;
      
    case STRUCTURAL_DISCIPLINE_INFO_ENUM_ID:
      info.reset(new Discipline::StructuralDisciplineInfo());
      break;

      case PISTON_THEORY_INFO_ENUM_ID:
        info.reset(new Discipline::PistonTheoryInfo());
        break;
        
    default:
      Assert(false, ExcInternalError());
    }
  
  return info;
}


std::istream& 
Discipline::operator>>(std::istream& input, Discipline::DisciplineInfo& info)
{
  info.readFromInputStream(input);
  return input;
}

