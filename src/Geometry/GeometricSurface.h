// $Id: GeometricSurface.h,v 1.2 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_surface_h__
#define __geometric_surface_h__


// C++ includes


// FESystem includes
#include "GeometricEntity.h"
#include "GeometricLine.h"

// libMesh includes


typedef std::pair<GeometricLine*, bool> LineAndOrientationPair;

/// this class defines a geometric surface, that is composed of lines. It could 
/// be 3 or 4 sided surface
class GeometricSurface : public GeometricEntity
{
public:

  /// constructor 
  GeometricSurface(GeometricModel& );

  /// destructor 
  virtual ~GeometricSurface();

  /// returns the type of this entity
  inline virtual GeometricEntityType type() const;

  /// returns the lines of this surface
  inline const std::vector<LineAndOrientationPair>& getLines() const;

  /// method to check if this object is a duplicate of the given object
  virtual bool checkIfDuplicate(const GeometricEntity& );

protected:

  /// method to read from input
  virtual void readFromInput(std::istream& );

  /// method to write to output 
  virtual void writeToOutput(std::ostream& );

  /// location of this point
  std::vector<LineAndOrientationPair> lines;
};


inline
GeometricEntity::GeometricEntityType GeometricSurface::type() const
{
  return GeometricEntity::SURFACE; 
}


inline
const std::vector<LineAndOrientationPair>& GeometricSurface::getLines() const
{
  return this->lines;
}

#endif // __geometric_surface_h__
