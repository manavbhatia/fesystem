// $Id: InterpolationCase.C,v 1.5.4.2 2007-06-13 14:57:20 manav Exp $

// C++ includes
#include <string>
#include <cassert>

// FESystem includes
#include "Interpolation/InterpolationCase.h"
#include "Utilities/InputOutputUtility.h"
#include "Discipline/AnalysisDisciplineBase.h"

// libMesh includes



InterpolationCase::InterpolationCase():
interpolation_case_ID(0),
interpolation_type(INVALID_INTERPOLATION_TYPE),
to_discipline_enum_ID(0),
from_discipline_enum_ID(0)
{

}

 
InterpolationCase::~InterpolationCase()
{
  
}

  
std::istream& InterpolationCase::readFromInputStream(std::istream& input)
{
  std::string tag;
  
  // should begin with a begin InterpolationCase tag
  FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASE");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  // read the interpolation case ID
  FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASE_ID", 
                                   this->interpolation_case_ID);

  // read the interpolation type
  tag.clear();
  FESystemIOUtility::readFromInput(input, "INTERPOLATION_TYPE", tag);
  
  if (tag == "DIRECT")
    this->interpolation_type = InterpolationCase::DIRECT;
  else if (tag == "FE")
    this->interpolation_type = InterpolationCase::FE;
  else
    abort();


  // read the mesh IDs
  tag.clear();
  FESystemIOUtility::readFromInput(input, "FROM_DISCIPLINE_ENUM_ID", tag);
  this->from_discipline_enum_ID = Discipline::AnalysisDisciplineEnum::enumID(tag);

  tag.clear();
  FESystemIOUtility::readFromInput(input, "TO_DISCIPLINE_ENUM_ID", tag);
  this->to_discipline_enum_ID = Discipline::AnalysisDisciplineEnum::enumID(tag);
  AssertThrow(this->from_discipline_enum_ID != this->to_discipline_enum_ID,
              ExcInternalError());

//  // read the to element pairs
//  unsigned int n_elem_sets = 0, elem_set_1 = 0, elem_set_2 = 0;
//
//  FESystemIOUtility::readFromInput(input, "ELEMENT_SET_PAIRS");
//  FESystemIOUtility::readFromInput(input, "BEGIN");
//  FESystemIOUtility::readFromInput(input, "N_ELEM_SETS", n_elem_sets);
//
//  for (unsigned int i=0; i<n_elem_sets; i++ )
//    {
//      input >> elem_set_1;
//      input >> elem_set_2;
//      this->elem_set_pairs.push_back(std::pair<unsigned int, unsigned int>
//				     (elem_set_1, elem_set_2));
//    }
//
//  // read the to mesh ID
//  FESystemIOUtility::readFromInput(input, "ELEMENT_SET_PAIRS");
//  FESystemIOUtility::readFromInput(input, "END");
  

  // end with the end interpolation tag
  FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASE");
  FESystemIOUtility::readFromInput(input, "END");

  return input;
}

