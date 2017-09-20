// $Id: load.C,v 1.6.6.2 2007-05-11 05:14:07 manav Exp $


// C++ includes


// FESystem includes
#include "Loads/load.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/PistonTheory.h"


unsigned int 
getLoadClassEnumIDForLoadName(unsigned int load_name_enum_ID)
{
  unsigned int load_class_enum_ID = FESystemNumbers::InvalidID;
  
  switch (load_name_enum_ID)
    {
    case NODAL_FORCE_ENUM_ID:
    case NODAL_TEMPERATURE_ENUM_ID:
    case NODAL_HEAT_LOAD_ENUM_ID:
      load_class_enum_ID = NODAL_POINT_LOAD::num();
      break;
      
    case SURFACE_PRESSURE_ENUM_ID:
      load_class_enum_ID = SCALAR_SURFACE_LOAD::num();
      break;

    case DISPLACEMENT_BOUNDARY_CONDITION_ENUM_ID:      
    case TEMPERATURE_BOUNDARY_CONDITION_ENUM_ID:
      load_class_enum_ID = DIRICHLET_BOUNDARY_CONDITION::num();
      break;

    case VOLUME_HEAT_LOAD_ENUM_ID:
      load_class_enum_ID = SCALAR_VOLUME_LOAD::num();
      break;

    case SURFACE_HEAT_LOAD_ENUM_ID:
      case PISTON_THEORY_SURFACE_ENUM_ID:
      load_class_enum_ID = SCALAR_SURFACE_LOAD::num();
      break;

    case SURFACE_RADIATION_HEAT_LOAD_ENUM_ID:
      load_class_enum_ID = SURFACE_RADIATION_LOAD::num();
      break;

    case SURFACE_CONVECTION_HEAT_LOAD_ENUM_ID:
      load_class_enum_ID = SURFACE_CONVECTION_LOAD::num();
      break;      
    
    default:
      Assert(false, ExcInternalError());
    }
  
  return load_class_enum_ID;
}


std::auto_ptr<SurfaceLoad> 
createSurfaceLoad(const unsigned int load_kind_enum_ID)
{
  std::auto_ptr<SurfaceLoad> load;

  switch (load_kind_enum_ID)
    {
    case SURFACE_RADIATION_LOAD_ENUM_ID:
      load.reset(new SurfaceRadiationLoad);
      break;

    case SCALAR_SURFACE_LOAD_ENUM_ID:
      load.reset(new ScalarSurfaceLoad);
      break;

    case SURFACE_CONVECTION_LOAD_ENUM_ID:
      load.reset(new SurfaceConvectionLoad);
      break;

    default:
      Assert(false, ExcInternalError());
    }

  return load;
}


std::auto_ptr<NodalLoad> 
createNodalLoad(const unsigned int nodal_load_kind_enum_ID,
                const unsigned int dim)
{
  std::auto_ptr<NodalLoad> load;

  switch (nodal_load_kind_enum_ID)
    {
    case NODAL_POINT_LOAD_ENUM_ID:
      {
        NodalPointLoad *ptr = new NodalPointLoad;
        ptr->setNDofs(dim);
        load.reset(ptr);
      }
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return load;
}



std::auto_ptr<VolumeLoad> 
createVolumeLoad(const unsigned int volume_load_kind_enum_ID)
{
  std::auto_ptr<VolumeLoad> load;
  
  switch (volume_load_kind_enum_ID)
    {
    case SCALAR_VOLUME_LOAD_ENUM_ID:
      load.reset(new ScalarVolumeLoad);
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return load;
}




std::auto_ptr<DirichletBoundaryCondition> 
createBoundaryConditionLoad(const unsigned int bc_load_kind_enum_ID)
{
  std::auto_ptr<DirichletBoundaryCondition> load;
  
  switch (bc_load_kind_enum_ID)
    {
    case DIRICHLET_BOUNDARY_CONDITION_ENUM_ID:
      load.reset(new DirichletBoundaryCondition);
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return load;
}

