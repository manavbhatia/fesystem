// $Id: LoadCombination.C,v 1.1.2.1 2007-05-15 20:38:52 manav Exp $

// C++ includes


// FESystem includes
#include "Loads/LoadCombination.h"

std::auto_ptr<Loads::VolumeLoadCombination>
Loads::createVolumeLoadCombination(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::VolumeLoadCombination> return_ptr;
  
  switch (load_type_enum_ID)
    {
    case SCALAR_VOLUME_LOAD_ENUM_ID:
      return_ptr.reset(new Loads::ScalarVolumeLoadCombination);
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return return_ptr;
}


std::auto_ptr<Loads::NodalLoadCombination>
Loads::createNodalLoadCombination(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::NodalLoadCombination> return_ptr;

  switch (load_type_enum_ID)
    {
    case NODAL_POINT_LOAD_ENUM_ID:
      return_ptr.reset(new Loads::NodalPointLoadCombination);
      break;

    default:
      Assert(false, ExcInternalError());
    }
  
  return return_ptr;
}



std::auto_ptr<Loads::DirichletBoundaryConditionCombination>
Loads::createDirichletBoundaryConditionCombination(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::DirichletBoundaryConditionCombination> return_ptr;
  
  switch (load_type_enum_ID)
    {
    case DIRICHLET_BOUNDARY_CONDITION_ENUM_ID:
      return_ptr.reset(new Loads::DirichletBoundaryConditionCombination);
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return return_ptr;
}


std::auto_ptr<Loads::SurfaceLoadCombination>
Loads::createSurfaceLoadCombination(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::SurfaceLoadCombination> return_ptr;

  switch (load_type_enum_ID)
    {
    case SCALAR_SURFACE_LOAD_ENUM_ID:
      return_ptr.reset(new Loads::ScalarSurfaceLoadCombination);
      break;

    case SURFACE_CONVECTION_LOAD_ENUM_ID:
      return_ptr.reset(new Loads::SurfaceConvectionLoadCombination);
      break;

    case SURFACE_RADIATION_LOAD_ENUM_ID:
      return_ptr.reset(new Loads::SurfaceRadiationLoadCombination);
      break;

    default:
      Assert(false, ExcInternalError());
    }
  
  return return_ptr;
}

