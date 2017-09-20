// $Id: GeometricCircularArc.C,v 1.3 2006-09-05 20:41:56 manav Exp $

// C++ includes


// FESystem includes
#include "Geometry/GeometricCircularArc.h"
#include "Geometry/GeometricModel.h"

// libMesh includes


GeometricCircularArc::GeometricCircularArc(GeometricModel& model)
  :GeometricLine(model)
{

}


GeometricCircularArc::~GeometricCircularArc()
{

}


bool GeometricCircularArc::checkIfDuplicate(const GeometricEntity& arc_ref)
{
  // make sure that the two arcs are not the same object
  assert (&arc_ref != this);

  // also, its type should be a circular arc
  if (arc_ref.type() != GeometricEntity::CIRCULAR_ARC)
    return false;

  // convert this line's reference to a geometric arc
  const GeometricCircularArc& arc = dynamic_cast<const GeometricCircularArc&>(arc_ref);

  // get the points for this arc
  const std::vector<GeometricPoint*>& arc_points = arc.getPoints();
  
  // compare the point IDs of this line with that of the given line
  // the comparisons should be made only with the base_entity_IDs
  unsigned int this_point1 = this->points[0]->baseEntityID();
  unsigned int this_point2 = this->points[2]->baseEntityID(); 
  unsigned int this_center = this->points[1]->baseEntityID();

  unsigned int given_point1 = arc_points[0]->baseEntityID();
  unsigned int given_point2 = arc_points[2]->baseEntityID(); 
  unsigned int given_center = arc_points[1]->baseEntityID();



  // if the two centers are not the same, then the arcs are different
  if (this_center != given_center)
    return false;

  bool dupl = false;

  // next compare the two ends, which would be the same as comparing a straight line
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


void GeometricCircularArc::readFromInput(std::istream& input)
{
  input >> this->unique_ID;

  unsigned int point_num = 0;

  GeometricPoint* point = NULL;

  for (unsigned int i=0; i<3; i++)
    {
      input >> point_num;
      point = &(this->geometric_model.getPointForID(point_num));
      this->points.push_back(point);
    }
}



void GeometricCircularArc::writeToOutput(std::ostream& output)
{
  output << this->unique_ID << "   ";
  output << this->points[0]->ID() << "   ";
  output << this->points[1]->ID() << "   ";
  output << this->points[2]->ID() << std::endl;
}
