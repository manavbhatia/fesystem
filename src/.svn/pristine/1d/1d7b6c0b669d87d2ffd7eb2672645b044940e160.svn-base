// $Id: conduction_1d.h,v 1.10 2006-10-29 05:09:10 manav Exp $

#ifndef __fesystem_conduction_1d_h__
#define __fesystem_conduction_1d_h__

// C++ includes


// FESystem includes
#include "ThermalElems/thermal_elem.h"

namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef THERMAL_CONDUCTION_EDGE2_ENUM_ID
#define THERMAL_CONDUCTION_EDGE2_ENUM_ID 1
#else
#error
#endif


#ifndef THERMAL_CONDUCTION_EDGE2_ENUM_NAME
#define THERMAL_CONDUCTION_EDGE2_ENUM_NAME "THERMAL_CONDUCTION_EDGE2"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(THERMAL_CONDUCTION_EDGE2,
                                THERMAL_CONDUCTION_EDGE2_ENUM_ID,
                                THERMAL_CONDUCTION_EDGE2_ENUM_NAME,
                                EDGE2)


namespace FESystemElem
{
  
/**
*	a 1-D conduction element
 */

  class Conduction_1D: public FESystemElem::ThermalElem
  {
public:
    Conduction_1D(Discipline::AnalysisDisciplineBase& discipline);
    
    ~Conduction_1D();
    
    /// @returns the number of dofs in the element
    virtual unsigned int getNDofs();
    
private:
	
  };
}


inline
unsigned int FESystemElem::Conduction_1D::getNDofs()
{
  return this->getNNodes();
}


#endif //__fesystem_conduction_1d_h__
