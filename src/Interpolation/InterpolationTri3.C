// $Id: InterpolationTri3.C,v 1.1 2007-01-16 04:09:29 manav Exp $

// C++ includes


// FESystem includes
#include "Interpolation/InterpolationTri3.h"

// libMesh includes
#include "geom/face_tri3.h"


InterpolationTri3::InterpolationTri3(Discipline::AnalysisDisciplineBase& discipline):
FEInterpolationElem(2, FESystemElem::FE_INTERPOLATION_TRI3::num(), discipline)
{
}

InterpolationTri3::~InterpolationTri3()
{
  
}
