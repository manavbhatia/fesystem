// $Id:$
/*
 *  PistonTheoryTri3.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 12/25/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

// FESystem includes
#include "PistonTheoryTri3.h"

// libMesh includes
#include "geom/face_tri3.h"

FESystemElem::PistonTheoryTri3::PistonTheoryTri3(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::PistonTheoryElem(2, FESystemElem::PISTON_THEORY_TRI3::num(), discipline)
{
}



FESystemElem::PistonTheoryTri3::~PistonTheoryTri3()
{
	
}


