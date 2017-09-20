// $Id: ElemSet.C,v 1.3 2006-10-29 05:09:18 manav Exp $

// C++ includes
#include <string>
#include <cassert>


// FESystem includes
#include "ElemSet.h"

// libMesh includes


ElemSet::ElemSet():
  set_ID(0),
  mesh_ID(0)
{

}



ElemSet::~ElemSet()
{

}


std::istream& ElemSet::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_elems = 0, elem_ID = 0;

  input >> tag;
  assert (tag == "BEGIN_ELEM_SET");
  input >> this->set_ID;
  
  tag.clear();
  input >> tag;
  assert (tag == "MESH_ID");
  input >> this->mesh_ID;

  tag.clear();
  input >> tag;
  assert (tag == "N_ELEMS");
  input >> n_elems;

  bool insert;

  for (unsigned int i=0; i<n_elems; i++)
    {
      input >> elem_ID;
      insert = this->elem_IDs.insert(elem_ID).second;
      assert (insert == true);
    }

  tag.clear();
  input >> tag;
  assert (tag== "END_ELEM_SET");

  return input;
}
