// $Id: GeometricEntity.h,v 1.2 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_entity_h__
#define __geometric_entity_h__ 

// C++ includes
#include <iostream>

// FESystem includes


// libMesh includes


// Forward declerations
class GeometricModel;


/// this is a base class for all geometric entities

class GeometricEntity
{
 public:

  /// type of the geometric entity
  enum GeometricEntityType
    {
      POINT,
      STRAIGHT_LINE,
      CIRCULAR_ARC,
      ELLIPTIC_ARC,
      SURFACE
    };

  /// constructor
  GeometricEntity(GeometricModel& );

  /// destructor 
  virtual ~GeometricEntity();

  /// returns the unique ID for this entity
  inline unsigned int ID() const;

  /// returns if this entity is base or duplicate. true indicates entity is a duplicate
  inline bool isDuplicate() const;

  /// return the ID of the base quantity if this quantity is a duplicate
  inline unsigned int baseEntityID() const;

  /// returns the orientation wrt the base elem, +1 for same, -1 for opposite
  inline int orientationWithBaseEntity() const;
  
  /// returns the type of this quantity
  virtual GeometricEntityType type() const = 0;

  /// this is an abstract method which checks if this object is a duplicate of the given object
  virtual bool checkIfDuplicate(const GeometricEntity& ) = 0;


 protected:

  /// abstract method to read from an input
  virtual void readFromInput(std::istream& ) = 0;

  /// abstract method to write to output
  virtual void writeToOutput(std::ostream& ) = 0;

  /// unique ID
  unsigned int unique_ID;

  /// base elem ID
  unsigned int base_entity_ID;

  /// if this entity has the same orientation as the base entity
  bool same_orientation_as_base_entity;

  /// indicates if duplicate
  bool is_duplicate;

  /// this is a reference to the geometric model to which this entity belongs. This is a very 
  /// expensive way of storing information, since each entity then keeps a reference to the 
  /// model, which is redundant. This will be changed in fututre, where approproate methods will
  /// be added to each entity
  GeometricModel& geometric_model;

  /// friend method to read from input stream
  inline friend std::istream& operator>> (std::istream& , GeometricEntity& );

  /// friend method to write to output stream
  inline friend std::ostream& operator<< (std::ostream&, GeometricEntity& );

};

inline 
unsigned int GeometricEntity::ID() const
{
  return this->unique_ID;
}

inline 
unsigned int GeometricEntity::baseEntityID() const
{
  // if this is a unique object, return the unique ID, else return the base entity ID
  if (this->is_duplicate == true)
    return this->base_entity_ID;
  else 
    return this->unique_ID;
}


inline
int GeometricEntity::orientationWithBaseEntity() const
{
  if (this->same_orientation_as_base_entity == true)
    return +1;
  else
    return -1;
}


inline
bool GeometricEntity::isDuplicate() const
{
  return this->is_duplicate;
}


inline
std::istream& operator>> (std::istream& input, GeometricEntity& entity)
{
  entity.readFromInput(input);
  return input;
}


inline
std::ostream& operator<< (std::ostream& output, GeometricEntity& entity)
{
  entity.writeToOutput(output);
  return output;
}

#endif // __geometric_entity_h__
