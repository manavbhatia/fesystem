// $Id: membrane_tri3.h,v 1.10 2006-10-29 05:09:01 manav Exp $

#ifndef __fesystem_membrane_tri3_h__
#define __fesystem_membrane_tri3_h__

// C++ includes

// FESystem includes
#include "StructuralElems/membrane.h"

namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef STRUCTURAL_MEMBRANE_TRI3_ENUM_ID
#define STRUCTURAL_MEMBRANE_TRI3_ENUM_ID 11
#else
#error
#endif


#ifndef STRUCTURAL_MEMBRANE_TRI3_ENUM_NAME
#define STRUCTURAL_MEMBRANE_TRI3_ENUM_NAME "STRUCTURAL_MEMBRANE_TRI3"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(STRUCTURAL_MEMBRANE_TRI3,
                                STRUCTURAL_MEMBRANE_TRI3_ENUM_ID,
                                STRUCTURAL_MEMBRANE_TRI3_ENUM_NAME,
                                TRI3)


namespace FESystemElem
{
  /**
  *	a class definig the quad membrane structural element
   */
  
  class MembraneTri3: public FESystemElem::Membrane
  {
public:
    MembraneTri3(Discipline::AnalysisDisciplineBase& discipline);
    
    ~MembraneTri3();
    
  };
}

#endif // __fesystem_membrane_tri3_h__
