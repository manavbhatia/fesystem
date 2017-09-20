// $Id: GeometricStraightLine.h,v 1.2 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_straight_line_h__
#define __geometric_straight_line_h__

// C++ includes
#include <vector>

// FESystem includes
#include "GeometricPoint.h"
#include "GeometricLine.h"

// libMesh includes


/// this class defines a geometric point
class GeometricStraightLine: public GeometricLine
{
 public:

  /// constructor
  GeometricStraightLine(GeometricModel& );

  /// destructor
  virtual ~GeometricStraightLine();

  /// returns the type of this entity
  inline virtual GeometricEntityType type() const;

  /// method to check if this iobject is a duplicate of the given object
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
GeometricEntity::GeometricEntityType GeometricStraightLine::type() const
{
  return GeometricEntity::STRAIGHT_LINE;
}

#endif // __geometric_straight_line_h__
