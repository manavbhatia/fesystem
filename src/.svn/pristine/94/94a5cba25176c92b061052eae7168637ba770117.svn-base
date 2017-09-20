// $Id: GeometricSet.h,v 1.2 2006-09-05 20:41:56 manav Exp $

#ifndef __geometric_set_h__
#define __geometric_set_h__


// C++ includes
#include <vector>
#include <string>

// FESystem includes



// libMesh includes


// forward declerations
class GeometricPoint;
class GeometricLine;
class GeometricSurface;


/// 
/// this class defines a set of points, lines, surfaces, etc that together form a set 
///

class GeometricSet
{
 public:
  /// constructor, takes the name of the set as an arguement
  GeometricSet(const std::string& );

  /// destructor
  ~GeometricSet();

  /// method to return the name of the set
  inline std::string getName() const;

  /// method to add a point to this set
  inline void addEntity(GeometricPoint* );

  /// method to add a line to this set
  inline void addEntity(GeometricLine* );

  /// method to add a point to this set
  inline void addEntity(GeometricSurface* );

  /// method to return the vector of points
  inline const std::vector<GeometricPoint*>& getPoints() const;

  /// method to return the vector of lines
  inline const std::vector<GeometricLine*>& getLines() const;

  /// method to return the vector of surfaces
  inline const std::vector<GeometricSurface*>& getSurfaces() const;
 
 protected:

  /// name of this set
  const std::string name;
  
  /// vector of points in this set
  std::vector<GeometricPoint*> points;

  /// vector of lines in this set
  std::vector<GeometricLine*> lines;

  /// vector of surfaces in this set
  std::vector<GeometricSurface*> surfaces;
};



inline std::string GeometricSet::getName() const
{return this->name; }


inline void GeometricSet::addEntity(GeometricPoint* point)
{
  this->points.push_back(point);
}


inline void GeometricSet::addEntity(GeometricLine* line)
{
  this->lines.push_back(line);
}


inline void GeometricSet::addEntity(GeometricSurface* surface)
{
  this->surfaces.push_back(surface);
}


inline const std::vector<GeometricPoint*>& GeometricSet::getPoints() const
{ return this->points;}


inline const std::vector<GeometricLine*>& GeometricSet::getLines() const
{ return this->lines;}


inline const std::vector<GeometricSurface*>& GeometricSet::getSurfaces() const
{ return this->surfaces;}




#endif // __geometric_set_h__
