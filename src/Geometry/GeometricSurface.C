// $Id: GeometricSurface.C,v 1.4 2006-09-05 20:41:56 manav Exp $

// C++ includes
#include <memory>

// FESystem includes
#include "GeometricSurface.h"
#include "GeometricModel.h"


// libMesh includes


GeometricSurface::GeometricSurface(GeometricModel& model)
  :GeometricEntity(model)
{

}


GeometricSurface::~GeometricSurface()
{

}

bool GeometricSurface::checkIfDuplicate(const GeometricEntity& surf_ref)
{
  // make sure that the two surfaces are not the same object
  assert (&surf_ref != this);

  const unsigned int n_lines = this->lines.size();

  // also, its type should be a surface
  if (surf_ref.type() != GeometricEntity::SURFACE)
    return false;

  // convert this entity's reference to a geometric surface
  const GeometricSurface& surf = dynamic_cast<const GeometricSurface&>(surf_ref);

  // get the lines for the given surface
  const std::vector<LineAndOrientationPair>& surf_lines = surf.getLines();
  if (n_lines != surf_lines.size())
    return false;

  // get the line IDs for the two surfaces, and take into account the orientation.
  // the tests for now will be simple enough where it is assumed that for the two surfaces, if the 
  // bounding lines are the same, they are duplicate. Then, if that is the case, then the orientation
  // of one of the lines in bot surfaces will give the orientation of the surface. i.e., if line 1 has ID 45 
  // in surface 1 and -45 in surface 2, then the two orientations are the opposite in direction.
  
  unsigned int id_sum=0;

  for (unsigned int i=0; i < n_lines; i++)
    for (unsigned int j=0; j < n_lines; j++)
      {
	if (this->lines[i].first->baseEntityID() == surf_lines[j].first->baseEntityID())
	  id_sum++;
      }

  if (id_sum > n_lines)
    abort(); // a surface cannot have multiple lines by same ID
  else if (id_sum < n_lines)
    return false; // the two surfaces do not share the same boundary
  else 
    {
      // this implies that id_sum == n_lines
      
      // now set the base entity ID and return
      this->is_duplicate = true;
      this->base_entity_ID = surf.baseEntityID();

      // check the orientation of any one line to establish the orientation of this surface
      // with respect to the given surface
      for (unsigned int j=0; j < n_lines; j++)
	if (this->lines[0].first->baseEntityID() == surf_lines[j].first->baseEntityID())
	  {
	    if (this->lines[0].second == surf_lines[j].second)
	      this->same_orientation_as_base_entity = true;
	    else 
	      this->same_orientation_as_base_entity = false;

	    break;
	  }

      return false;
    }
}


void GeometricSurface::readFromInput(std::istream& input)
{
  input >> this->unique_ID;

  unsigned int n_lines = 0;
  int line_num=0;
  input >> n_lines;

  GeometricLine* line = NULL;

  for (unsigned int i=0; i<n_lines; i++)
    {
      input >> line_num;
      line = &(this->geometric_model.getLineForID(line_num));
      if (line_num > 0) // same orientation
	this->lines.push_back(std::make_pair<GeometricLine*, bool>(line,true));
      else // opposite orientation
	this->lines.push_back(std::make_pair<GeometricLine*, bool>(line,false));
    }
}



void GeometricSurface::writeToOutput(std::ostream& output)
{
  output << this->unique_ID << "   ";
  unsigned int n_lines = this->lines.size();
  int line_num = 0;
  
  for (unsigned int i=0; i < n_lines; i++)
    {
      line_num = this->lines[i].first->ID();
      if (this->lines[i].second == false)
	line_num *= -1;
      output << line_num << "   ";
    }

  output << std::endl;
}
