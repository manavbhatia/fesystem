// $Id: conduction_1d.C,v 1.7 2006-09-05 20:41:39 manav Exp $

// C++ includes


// FESystem includes
#include "ThermalElems/conduction_1d.h"

// libMesh includes
#include "geom/edge_edge2.h"

FESystemElem::Conduction_1D::Conduction_1D(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::ThermalElem(1, FESystemElem::THERMAL_CONDUCTION_EDGE2::num(), discipline)
{
}



FESystemElem::Conduction_1D::~Conduction_1D()
{
	
}
