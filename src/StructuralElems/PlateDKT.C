// $Id: PlateDKT.C,v 1.1 2006-09-10 05:49:03 manav Exp $

// C++ includes

// FESystem includes
#include "StructuralElems/PlateDKT.h"


FESystemElem::PlateDKT::PlateDKT(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::PlateDKBatoz(2, FESystemElem::STRUCTURAL_PLATE_DKT::num(), discipline)
{
  
}





FESystemElem::PlateDKT::~PlateDKT()
{
	
}

