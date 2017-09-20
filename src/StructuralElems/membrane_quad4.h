// $Id: membrane_quad4.h,v 1.14 2006-10-29 05:09:01 manav Exp $

#ifndef __fesystem_membrane_quad4_h__
#define __fesystem_membrane_quad4_h__

// C++ includes

// FESystem includes
#include "StructuralElems/membrane.h"

namespace Discipline
{
  class AnalysisDisciplineBase;
}



#ifndef STRUCTURAL_MEMBRANE_QUAD4_ENUM_ID
#define STRUCTURAL_MEMBRANE_QUAD4_ENUM_ID 10
#else
#error
#endif


#ifndef STRUCTURAL_MEMBRANE_QUAD4_ENUM_NAME
#define STRUCTURAL_MEMBRANE_QUAD4_ENUM_NAME "STRUCTURAL_MEMBRANE_QUAD4"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(STRUCTURAL_MEMBRANE_QUAD4,
                                STRUCTURAL_MEMBRANE_QUAD4_ENUM_ID,
                                STRUCTURAL_MEMBRANE_QUAD4_ENUM_NAME,
                                QUAD4)



namespace FESystemElem
{
  /**
  *	a class definig the quad membrane structural element
   */
  
  class MembraneQuad4: public FESystemElem::Membrane
  {
public:
    MembraneQuad4(Discipline::AnalysisDisciplineBase& discipline);
    
    ~MembraneQuad4();
    
  };
}

#endif // __fesystem_membrane_quad4_h__
