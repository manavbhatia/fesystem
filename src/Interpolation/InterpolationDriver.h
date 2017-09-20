// $Id: InterpolationDriver.h,v 1.5 2006-09-05 20:41:48 manav Exp $

#ifndef __interpolation_driver_h__
#define __interpolation_driver_h__

// C++ includes
#include <vector>
#include <memory>

// FESystem includes


// libMesh includes


// Forward declerations
namespace FESystem
{
  class FESystemController;
}

namespace MeshDS
{
  class FEMesh;
  class FEMeshData;
}

class InterpolationCase;
class InterpolationBase;
class DofMap;

/// this class will handle all aspects of interpolation. It will take the fesystem controller
/// and the interpolation case as an argurment, and based on the interpolation kind requested
/// it will perform the interpolation 
class InterpolationDriver
{
public:
  /// constructor
  /// @param controller the FESystemController class that owns this object
  /// #param case InterpolatinCase that this driver will solve
  InterpolationDriver(FESystem::FESystemController& controller, 
                      InterpolationCase& interp_case );


  /// destructor
  ~InterpolationDriver();

  
  /// method to interpolate and create the loads for next analysis
  void interpolateAndCreateLoads(const std::vector<unsigned int>& load_cases);

  /// returns the fesystem controller for this class
  inline  FESystem::FESystemController& getFESystemController();

  /// returns the interpolation case for this driver
  inline InterpolationCase& getInterpolationCase();

  /// returns the interpolation mesh
  inline MeshDS::FEMesh* getInterpolationMesh();

  
  /// returns the interpolation mesh data
  inline MeshDS::FEMeshData* getInterpolationMeshData();


  /// returns the interpolation mesh ID
  inline unsigned int getInterpolationMeshID();


protected:

  /// initialization function
  void init();

  /// interpolates the laods for all load cases and DVs and
  /// stores the interpolated solution in the database
  void interpolate(const std::vector<unsigned int>& load_cases);


  /// writes loads from the interpolated values
  void createLoads(const std::vector<unsigned int>& load_cases);


  /// FESystemController object that owns this object
  FESystem::FESystemController& fesystem_controller;

  /// InterpolationCase for this driver
  InterpolationCase& interpolation_case;

  /// the interpolation base that will process the interpolation
  std::auto_ptr<InterpolationBase> interpolation_base;

  /// mesh ID
  unsigned int mesh_ID;

  /// mesh to which data will be interpolated
  MeshDS::FEMesh* mesh;

  /// mesh data to which data will be interpolated
  MeshDS::FEMeshData* mesh_data;

  /// dofmap of the mesh to which the data will be interpolated
  DofMap* dof_map;
};




inline
FESystem::FESystemController& 
InterpolationDriver::getFESystemController()
{
  return this->fesystem_controller;
}


inline
InterpolationCase&
InterpolationDriver::getInterpolationCase()
{
  return this->interpolation_case;
}



inline 
MeshDS::FEMesh* 
InterpolationDriver::getInterpolationMesh()
{
  return this->mesh;
}

  

inline 
MeshDS::FEMeshData*
InterpolationDriver::getInterpolationMeshData()
{
  return this->mesh_data;
}


inline 
unsigned int
InterpolationDriver::getInterpolationMeshID()
{
  return this->mesh_ID;
}


#endif // __interpolation_driver_h__
