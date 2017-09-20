// $Id: InterpolationCase.h,v 1.3.4.2 2007-06-13 14:57:20 manav Exp $

#ifndef __fesystem_interpolation_case_h__
#define __fesystem_interpolation_case_h__


//  C++ includes
#include <iostream>
#include <vector>

// FESystem includes


// libMesh Includes



/// this defines an interpolation case that will contain the necessary 
/// data to interpolate the data from a mesh to another 

class InterpolationCase
{
 public:
  enum InterpolationType{
    DIRECT,
    FE,
    INVALID_INTERPOLATION_TYPE 
  };

  /// constructor
  InterpolationCase();

  /// destructor
  ~InterpolationCase();

  /// method to read from an input stream
  /// @param input input stream from which the data is read
  std::istream& readFromInputStream(std::istream& input);
  
  /// returns the type of interpolation used for this case
  InterpolationCase::InterpolationType type() const;

  /// returns case ID
  unsigned int ID() const;
  
//  /// returns a reference ot the element set pairs for this case
//  inline std::vector<std::pair<unsigned int, unsigned int> >& getElemSetPairs();


  /// returns the mesh from which the data is interpolated
  unsigned int fromDiscipline() const;


  /// returns the mesh to which the data is interpolated
  unsigned int toDiscipline() const;
  

 protected:
  
  /// case ID
  unsigned int interpolation_case_ID;

  /// Interpolation type
  InterpolationType interpolation_type;

  /// to mesh ID
  unsigned int to_discipline_enum_ID;
  
  /// from mesh ID
  unsigned int from_discipline_enum_ID;
  
//  /// elem set ID pairs
//  std::vector<std::pair<unsigned int, unsigned int> >
//    elem_set_pairs;

};


inline 
InterpolationCase::InterpolationType InterpolationCase::type() const
{
  return this->interpolation_type;
}



inline
unsigned int InterpolationCase::ID() const
{
  return this->interpolation_case_ID;
}


//inline
//std::vector<std::pair<unsigned int, unsigned int> >& 
//InterpolationCase::getElemSetPairs() 
//{
//  return this->elem_set_pairs;
//}



inline 
unsigned int InterpolationCase::fromDiscipline() const
{
  return this->from_discipline_enum_ID;
}


inline 
unsigned int InterpolationCase::toDiscipline() const
{
  return this->to_discipline_enum_ID;
}


#endif // __fesystem_interpolation_case_h__
