// $Id: GeometricEllipticArc.C,v 1.3 2006-09-05 20:41:56 manav Exp $

// C++ includes


// FESystem includes
#include "Geometry/GeometricEllipticArc.h"
#include "Geometry/GeometricModel.h"

// libMesh includes
#include "geom/point.h"
#include "numerics/type_vector.h"


GeometricEllipticArc::GeometricEllipticArc(GeometricModel& model)
  :GeometricLine(model)
{

}


GeometricEllipticArc::~GeometricEllipticArc()
{

}


bool GeometricEllipticArc::checkIfDuplicate(const GeometricEntity& arc_ref)
{
  // make sure that the two arcs are not the same object
  assert (&arc_ref != this);

  // also, its type should be an elliptic arc
  if (arc_ref.type() != GeometricEntity::ELLIPTIC_ARC)
    return false;

  // convert this line's reference to a geometric arc
  const GeometricEllipticArc& arc = dynamic_cast<const GeometricEllipticArc&>(arc_ref);

  // get the points for this arc
  const std::vector<GeometricPoint*>& arc_points = arc.getPoints();
  
  // compare the point IDs of this line with that of the given line
  // the comparisons should be made only with the base entity IDs
  unsigned int this_point1 = this->points[0]->baseEntityID();
  unsigned int this_point2 = this->points[3]->baseEntityID();
  unsigned int this_center = this->points[1]->baseEntityID();
  unsigned int this_major_axis_point = this->points[2]->baseEntityID();
  // calculate the unit vector along the mahor axis
  Location center, major_axis_point;
  center = this->geometric_model.getPointForID(this_center).getLocation();
  major_axis_point = this->geometric_model.getPointForID(this_major_axis_point).getLocation();

  Point this_major_axis_unit_vec(major_axis_point.coord1 - center.coord1,
				 major_axis_point.coord2 - center.coord2,
				 major_axis_point.coord3 - center.coord3);
  TypeVector<double> this_unit_vec(this_major_axis_unit_vec.unit());

  unsigned int given_point1 = arc_points[0]->baseEntityID();
  unsigned int given_point2 = arc_points[3]->baseEntityID();
  unsigned int given_center = arc_points[1]->baseEntityID();
  unsigned int given_major_axis_point = arc_points[2]->baseEntityID();
  // calculate the unit vector along the major axis
  center = this->geometric_model.getPointForID(given_center).getLocation();
  major_axis_point = this->geometric_model.getPointForID(given_major_axis_point).getLocation();

  Point given_major_axis_unit_vec(major_axis_point.coord1 - center.coord1,
				 major_axis_point.coord2 - center.coord2,
				 major_axis_point.coord3 - center.coord3);
  TypeVector<double> given_unit_vec(given_major_axis_unit_vec.unit());


  // if the two centers are not the same, the the arcs are different
  if (this_center != given_center)
    return false;

  // next compare the unit vectors. They should be the same vectors, either parallel or 
  // anti-parallel
  double dot_prod = this_unit_vec * given_unit_vec;
  double machine_eps = 1.0e-12;

  // if the two vectors do not have a unity dot product, they are not coincidental
  if ( fabs(fabs(dot_prod)-1.0) > machine_eps)
    return false;
  
  bool dupl = false;

  // finally compare the end points
  if (this_point1 == given_point1 && this_point2 == given_point2)
    {
      dupl = true;
      this->same_orientation_as_base_entity = true;
    }
  else if (this_point2 == given_point1 && this_point1 == given_point2)
    {
      dupl = true;
      this->same_orientation_as_base_entity = false;
    }

  if (dupl == true)
    {
      this->is_duplicate = true;
      this->base_entity_ID = arc.baseEntityID();

      return true;
    }
  else 
    return false;
}



void GeometricEllipticArc::readFromInput(std::istream& input)
{
  input >> this->unique_ID;

  unsigned int point_num = 0;

  GeometricPoint* point = NULL;

  for (unsigned int i=0; i<4; i++)
    {
      input >> point_num;
      point = &(this->geometric_model.getPointForID(point_num));
      this->points.push_back(point);
    }
}



void GeometricEllipticArc::writeToOutput(std::ostream& output)
{
  output << this->unique_ID << "   ";
  output << this->points[0]->ID() << "   ";
  output << this->points[1]->ID() << "   ";
  output << this->points[2]->ID() << "   ";
  output << this->points[3]->ID() << std::endl;
}
