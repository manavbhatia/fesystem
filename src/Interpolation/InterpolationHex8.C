// $Id: InterpolationHex8.C,v 1.5 2006-09-05 20:41:48 manav Exp $

// C++ includes


// FESystem includes
#include "Interpolation/InterpolationHex8.h"

// libMesh includes
#include "geom/cell_hex8.h"


InterpolationHex8::InterpolationHex8(Discipline::AnalysisDisciplineBase& discipline):
FEInterpolationElem(3, FESystemElem::FE_INTERPOLATION_HEX8::num(), discipline)
{

}

InterpolationHex8::~InterpolationHex8()
{
	
}
