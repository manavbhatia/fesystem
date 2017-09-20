/*
 *  PistonTheoryInfo.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */


// $Id: PistonTheoryInfo.C,v 1.4 2006/11/13 10:11:57 manav Exp $

// FESystem includes
#include "Discipline/PistonTheoryInfo.h"
#include "Utilities/InputOutputUtility.h"


Discipline::PistonTheoryInfo::PistonTheoryInfo():
DisciplineInfo(),
theory_order(0),
element_set_ID(0)
{
  
}



Discipline::PistonTheoryInfo::~PistonTheoryInfo()
{
}



std::istream& 
Discipline::PistonTheoryInfo::readFromInputStream(std::istream& input)
{
  
  FESystemIOUtility::readFromInput(input, Discipline::PISTON_THEORY_INFO::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  // this is the ID of the structural mesh on which the analysis will be performed
  FESystemIOUtility::readFromInput(input, "MESH_ID", this->mesh_ID);
    
  /// order of the piston theory
  FESystemIOUtility::readFromInput(input, "ORDER", this->theory_order);
    
  FESystemIOUtility::readFromInput(input, Discipline::PISTON_THEORY_INFO::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  
  return input;
}



