// $Id: membrane_tri3.C,v 1.5 2006-09-05 20:41:58 manav Exp $

// C++ includes


// FESystem includes
#include "StructuralElems/membrane_tri3.h"

// libMesh includes
#include "geom/face_tri3.h"


FESystemElem::MembraneTri3::MembraneTri3(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::Membrane(2, FESystemElem::STRUCTURAL_MEMBRANE_TRI3::num(), discipline)
{

}





FESystemElem::MembraneTri3::~MembraneTri3()
{
	
}


