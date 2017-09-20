// $Id: ElemSet.h,v 1.3 2006-09-05 20:41:44 manav Exp $

#ifndef __elem_set_h__
#define __elem_set_h__

// C++ includes
#include <set>
#include <iostream>


// FESystem includes


// libMesh includes


/// this class defines a set of elements on a mesh
class ElemSet
{
 public:

  /// constructor
  ElemSet();

  /// destructor
  ~ElemSet();

  /// returns the ID of this elem set
  unsigned int ID() const;

  
  /// get mesh ID to which these elements belong
  inline unsigned int getMeshID() const;

  /// get the set of elems
  inline const std::set<unsigned int>& getElemIDs() const ;

  /// read from input stream
  /// @param input an istream from which the data for this
  /// element set will be read
  std::istream& readFromInputStream(std::istream& input);

  /// returns true if the given elem ID is contained in the set
  /// @param elem_ID ID of the elem
  /// @param mesh_ID ID of the mesh that contains the elem
  inline bool containsElem(unsigned int elem_ID, 
			   unsigned int mesh_ID) const;
  
 protected:
  /// set ID
  unsigned int set_ID;

  /// Mesh ID to which this belongs
  unsigned int mesh_ID;

  /// set of elem IDs
  std::set<unsigned int> elem_IDs;

};



inline 
unsigned int ElemSet::ID() const
{
  return this->set_ID;
}


inline
unsigned int ElemSet::getMeshID() const
{
  return this->mesh_ID;
}


inline 
const std::set<unsigned int>& ElemSet::getElemIDs() const
{
  return this->elem_IDs;
}

inline 
bool ElemSet::containsElem(unsigned int elem,
			   unsigned int mesh) const
{
  if (mesh != this->mesh_ID)
    return false;

  std::set<unsigned int>::const_iterator it = 
    this->elem_IDs.find(elem);
  
  if (it == this->elem_IDs.end())
    return false;
  else
    return true;
}


#endif // __elem_set_h__
