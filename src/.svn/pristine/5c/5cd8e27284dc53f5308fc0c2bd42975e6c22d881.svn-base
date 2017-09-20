// $Id: conduction_tri3.h,v 1.10 2006-10-29 05:09:10 manav Exp $

#ifndef __fesystem_conduction_tri3_h__
#define __fesystem_conduction_tri3_h__

// C++ includes


// FESystem includes
#include "ThermalElems/thermal_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}



#ifndef THERMAL_CONDUCTION_TRI3_ENUM_ID
#define THERMAL_CONDUCTION_TRI3_ENUM_ID 5
#else
#error
#endif


#ifndef THERMAL_CONDUCTION_TRI3_ENUM_NAME
#define THERMAL_CONDUCTION_TRI3_ENUM_NAME "THERMAL_CONDUCTION_TRI3"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(THERMAL_CONDUCTION_TRI3,
                                THERMAL_CONDUCTION_TRI3_ENUM_ID,
                                THERMAL_CONDUCTION_TRI3_ENUM_NAME,
                                TRI3)



namespace FESystemElem
{
  /**
  *	a 2D 3 noded triangular conduction element
   */
  
  class Conduction_Tri3: public FESystemElem::ThermalElem
  {
public:
    Conduction_Tri3(Discipline::AnalysisDisciplineBase& discipline);
    
    ~Conduction_Tri3();
    
    /// @returns the number of dofs in the element
    virtual unsigned int getNDofs();
    
private:
      
  };
}



inline
unsigned int FESystemElem::Conduction_Tri3::getNDofs()
{
  return this->getNNodes();
}


#endif //__fesystem_conduction_tri3_h__
