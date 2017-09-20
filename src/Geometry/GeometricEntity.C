// $Id: GeometricEntity.C,v 1.3 2006-09-05 20:41:56 manav Exp $

// C++ includes


// FESystem includes
#include "Geometry/GeometricEntity.h"

// libMesh includes



GeometricEntity::GeometricEntity(GeometricModel& model):
  unique_ID(0),
  base_entity_ID(0),
  same_orientation_as_base_entity(false),
  is_duplicate(false),
  geometric_model(model)
{

}


GeometricEntity::~GeometricEntity()
{

}
