// $Id: GeometricStraightLine.C,v 1.2 2006-09-05 20:41:56 manav Exp $

// C++ includes


// FESystem includes
#include "GeometricStraightLine.h"
#include "GeometricModel.h"

// libMesh includes





GeometricStraightLine::GeometricStraightLine(GeometricModel& model)
  :GeometricLine(model)
{

}


GeometricStraightLine::~GeometricStraightLine()
{

}


bool GeometricStraightLine::checkIfDuplicate(const GeometricEntity& line_ref)
{
  // make sure that the two lines are not the same object
  assert (&line_ref != this);
  
  // also, this type should be a line
  if (line_ref.type() != GeometricEntity::STRAIGHT_LINE)
    return false;

  // convert this line's reference to a geometric line
  const GeometricStraightLine& line = dynamic_cast<const GeometricStraightLine&>(line_ref);

  // get the points of this line
  const std::vector<GeometricPoint*>& line_points = line.getPoints();

  // compare the points IDs of this line with that of the given line
  // the comparisons should be made only with the base_entity_IDs
  unsigned int this_point1 = this->points[0]->baseEntityID();
  unsigned int this_point2 = this->points[1]->baseEntityID();

  unsigned int given_point1 = line_points[0]->baseEntityID();
  unsigned int given_point2 = line_points[1]->baseEntityID();

  bool dupl = false;


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
      this->base_entity_ID = line.baseEntityID();

      return true;
    }
  else 
    return  false;
}


void GeometricStraightLine::readFromInput(std::istream& input)
{
  input >> this->unique_ID;

  unsigned int point_num = 0;
  // first point
  input >> point_num;
  GeometricPoint  *point = &(this->geometric_model.getPointForID(point_num));
  this->points.push_back(point);
  
  // second point
  input >> point_num;
  point = &(this->geometric_model.getPointForID(point_num));
  this->points.push_back(point);  
}


void GeometricStraightLine::writeToOutput(std::ostream& output)
{
  output << this->unique_ID << "   ";
  output << this->points[0]->ID() << "   ";
  output << this->points[1]->ID() << std::endl;
}
