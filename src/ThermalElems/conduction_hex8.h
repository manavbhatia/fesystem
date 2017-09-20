// $Id: conduction_hex8.h,v 1.10 2006-10-29 05:09:10 manav Exp $

#ifndef __fesystem_conduction_hex8_h__
#define __fesystem_conduction_hex8_h__

// C++ includes


// FESystem includes
#include "ThermalElems/thermal_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef THERMAL_CONDUCTION_HEX8_ENUM_ID
#define THERMAL_CONDUCTION_HEX8_ENUM_ID 2
#else
#error
#endif


#ifndef THERMAL_CONDUCTION_HEX8_ENUM_NAME
#define THERMAL_CONDUCTION_HEX8_ENUM_NAME "THERMAL_CONDUCTION_HEX8"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(THERMAL_CONDUCTION_HEX8,
                                THERMAL_CONDUCTION_HEX8_ENUM_ID,
                                THERMAL_CONDUCTION_HEX8_ENUM_NAME,
                                HEX8)


namespace FESystemElem
{
  
/**
*	a HEX8 conduction element
 */

  class Conduction_Hex8: public FESystemElem::ThermalElem
  {
public:
    Conduction_Hex8(Discipline::AnalysisDisciplineBase& discipline);
    
    ~Conduction_Hex8();
    
    /// @returns the number of dofs in the element
    virtual unsigned int getNDofs();
    
    
private:
      
  };
}


inline
unsigned int FESystemElem::Conduction_Hex8::getNDofs()
{
  return this->getNNodes();
}


#endif //__fesystem_conduction_1d_h__
