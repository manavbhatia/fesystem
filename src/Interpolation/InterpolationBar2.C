// $Id: InterpolationBar2.C,v 1.5 2006-09-05 20:41:48 manav Exp $

// C++ includes


// FESystem includes
#include "Interpolation/InterpolationBar2.h"

// libMesh includes
#include "geom/edge_edge2.h"

InterpolationBar2::InterpolationBar2(Discipline::AnalysisDisciplineBase& discipline):
FEInterpolationElem(1, FESystemElem::FE_INTERPOLATION_EDGE2::num(), discipline)
{
}

InterpolationBar2::~InterpolationBar2()
{
	
}
