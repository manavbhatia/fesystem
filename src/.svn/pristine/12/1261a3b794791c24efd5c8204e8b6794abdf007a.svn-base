// $Id:$
/*
 *  PistonTheoryTri3.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 12/25/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

#ifndef __fesystem_piston_theory_tri3_elem_h__
#define __fesystem_piston_theory_tri3_elem_h__

// FESystem includes
#include "PanelMethods/PistonTheoryElem.h"

namespace Discipline
{
  class AnalysisDisciplineBase;
}



#ifndef PISTON_THEORY_TRI3_ENUM_ID
#define PISTON_THEORY_TRI3_ENUM_ID 18
#else
#error
#endif


#ifndef PISTON_THEORY_TRI3_ENUM_NAME
#define PISTON_THEORY_TRI3_ENUM_NAME "PISTON_THEORY_TRI3"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(PISTON_THEORY_TRI3,
                                PISTON_THEORY_TRI3_ENUM_ID,
                                PISTON_THEORY_TRI3_ENUM_NAME,
                                TRI3)


namespace FESystemElem
{
  /**
   *	a 2D 3 noded triangular conduction element
   */
  
  class PistonTheoryTri3: public FESystemElem::PistonTheoryElem
  {
  public:
    PistonTheoryTri3(Discipline::AnalysisDisciplineBase& discipline);
    
    ~PistonTheoryTri3();
        
  private:
    
  };
}




#endif //__fesystem_piston_theory_tri3_elem_h__

