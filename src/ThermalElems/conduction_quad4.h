// $Id: conduction_quad4.h,v 1.10 2006-10-29 05:09:10 manav Exp $

#ifndef __fesystem_conduction_quad4_h__
#define __fesystem_conduction_quad4_h__

// C++ includes


// FESystem includes
#include "ThermalElems/thermal_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}



#ifndef THERMAL_CONDUCTION_QUAD4_ENUM_ID
#define THERMAL_CONDUCTION_QUAD4_ENUM_ID 4
#else
#error
#endif


#ifndef THERMAL_CONDUCTION_QUAD4_ENUM_NAME
#define THERMAL_CONDUCTION_QUAD4_ENUM_NAME "THERMAL_CONDUCTION_QUAD4"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(THERMAL_CONDUCTION_QUAD4,
                                THERMAL_CONDUCTION_QUAD4_ENUM_ID,
                                THERMAL_CONDUCTION_QUAD4_ENUM_NAME,
                                QUAD4)



namespace FESystemElem
{
  /**
  *	a QUAD4 conduction element
   */

  class Conduction_Quad4: public FESystemElem::ThermalElem
  {
public:
    Conduction_Quad4(Discipline::AnalysisDisciplineBase& discipline);
    
    ~Conduction_Quad4();
    
    /// @returns the number of dofs in the element
    virtual unsigned int getNDofs();
    
private:
      
  };
}


inline
unsigned int FESystemElem::Conduction_Quad4::getNDofs()
{
  return this->getNNodes();
}


#endif //__fesystem_conduction_quad4_h__
