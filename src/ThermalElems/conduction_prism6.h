// $Id: conduction_prism6.h,v 1.9 2006-10-29 05:09:10 manav Exp $

#ifndef __fesystem_conduction_prism6_h__
#define __fesystem_conduction_prism6_h__

// C++ includes


// FESystem includes
#include "ThermalElems/thermal_elem.h"

// libMesh includes

namespace Discipline
{
  class AnalysisDisciplineBase;
}



#ifndef THERMAL_CONDUCTION_PRISM6_ENUM_ID
#define THERMAL_CONDUCTION_PRISM6_ENUM_ID 3
#else
#error
#endif


#ifndef THERMAL_CONDUCTION_PRISM6_ENUM_NAME
#define THERMAL_CONDUCTION_PRISM6_ENUM_NAME "THERMAL_CONDUCTION_PRISM6"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(THERMAL_CONDUCTION_PRISM6,
                                THERMAL_CONDUCTION_PRISM6_ENUM_ID,
                                THERMAL_CONDUCTION_PRISM6_ENUM_NAME,
                                PRISM6)


namespace FESystemElem
{
/**
*	a PRISM6 conduction element
 */

  class Conduction_Prism6: public FESystemElem::ThermalElem
  {
public:
    Conduction_Prism6(Discipline::AnalysisDisciplineBase& discipline);
    
    ~Conduction_Prism6();
    
    /// @returns the number of dofs in the element
    virtual  unsigned int getNDofs();
    
private:
		
  };
}



inline
unsigned int FESystemElem::Conduction_Prism6::getNDofs()
{
  return this->getNNodes();
}


#endif //__fesystem_conduction_prism6_h__
