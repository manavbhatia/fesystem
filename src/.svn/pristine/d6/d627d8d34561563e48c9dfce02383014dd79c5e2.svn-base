// $Id: GeometricCircularArc.h,v 1.2 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_circular_arc_h__
#define __geometric_circular_arc_h__

// C++ includes
#include <vector>

// FESystem includes
#include "GeometricLine.h"
#include "GeometricPoint.h"

// libMesh includes


/// this class defines a geometric point
class GeometricCircularArc: public GeometricLine
{
 public:

  /// constructor
  GeometricCircularArc(GeometricModel& );

  /// destructor
  virtual ~GeometricCircularArc();

  /// returns the type of this entity
  inline virtual GeometricEntityType type() const;

  /// method to check if this object is a duplicate of the given object
  virtual bool checkIfDuplicate(const GeometricEntity& );

 protected:

  /// method to read from input
  virtual void readFromInput(std::istream& );

  /// method to write to output
  virtual void writeToOutput(std::ostream& );

  /// location of this point
  std::vector<GeometricPoint*> points;
};


inline
GeometricEntity::GeometricEntityType GeometricCircularArc::type() const
{
  return GeometricEntity::CIRCULAR_ARC;
}

#endif // __geometric_circular_arc_h__
