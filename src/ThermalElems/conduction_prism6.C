// $Id: conduction_prism6.C,v 1.5 2006-09-05 20:41:39 manav Exp $

// C++ includes


// FESystem includes
#include "ThermalElems/conduction_prism6.h"

// libMesh includes
#include "geom/cell_prism6.h"


FESystemElem::Conduction_Prism6::Conduction_Prism6(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::ThermalElem(3, FESystemElem::THERMAL_CONDUCTION_PRISM6::num(), discipline)
{
}


FESystemElem::Conduction_Prism6::~Conduction_Prism6()
{
	
}
