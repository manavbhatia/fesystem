// $Id: PlateDKT.h,v 1.3 2006-10-29 05:09:01 manav Exp $

#ifndef __fesystem_plate_DKT_h__
#define __fesystem_plate_DKT_h__

// C++ includes

// FESystem includes
#include "StructuralElems/PlateDKBatoz.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef STRUCTURAL_PLATE_DKT_ENUM_ID
#define STRUCTURAL_PLATE_DKT_ENUM_ID 16
#else
#error
#endif


#ifndef STRUCTURAL_PLATE_DKT_ENUM_NAME
#define STRUCTURAL_PLATE_DKT_ENUM_NAME "STRUCTURAL_PLATE_DKT"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(STRUCTURAL_PLATE_DKT,
                                STRUCTURAL_PLATE_DKT_ENUM_ID,
                                STRUCTURAL_PLATE_DKT_ENUM_NAME,
                                TRI3)





namespace FESystemElem
{
  
  
  class PlateDKT: public FESystemElem::PlateDKBatoz
  {
public:
    
    PlateDKT(Discipline::AnalysisDisciplineBase& discipline);
    
    ~PlateDKT();
    
    
    
        
protected:
      
    
  };
}


#endif // __fesystem_plate_DKT_h__
