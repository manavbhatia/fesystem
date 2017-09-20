// $Id: InterpolationBase.h,v 1.5.4.2 2007-06-13 14:57:20 manav Exp $

#ifndef __fesystem_interpolation_base_h__
#define __fesystem_interpolation_base_h__


// C++ includes


// FESystem includes


// libMesh includes
#include "numerics/numeric_vector.h"

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
class DofMap;

// Forward decleration
class InterpolationDriver;
class Point;
class Elem;

/// this base class defines the interpolation interface for various
/// interpolation methods
class InterpolationBase
{
public:
  /// constructor
  InterpolationBase(FESystem::FESystemController& controller, 
                    const InterpolationCase& interp_case);
  
  /// destructor 
  virtual ~InterpolationBase();
  
  
  /// returns the fesystem controller for this class
  FESystem::FESystemController& getFESystemController();
  
  /// returns the interpolation case for this driver
  const InterpolationCase& getInterpolationCase();
  
//  /// returns the interpolation mesh
//  MeshDS::FEMesh* getInterpolationMesh();
//  
//  
//  /// returns the interpolation mesh data
//  MeshDS::FEMeshData* getInterpolationMeshData();
//  
//  
//  /// returns the interpolation mesh ID
//  unsigned int getInterpolationMeshID();
  
  /// interpolates the laods for all load cases and DVs
  virtual void interpolate() = 0;

protected:
    
//  /// initialization function
//  void init();
    
  
//  /// returns the element which contains this point. This method searches 
//  /// only in the element set that contains this element
//  /// @param point the point which has to be mapped
//  /// @param elem_ID the elem in the interpolation mesh which contanis this point
//  /// @param mesh_ID the mesh which contains elem.
//  Elem* getElemContainingPoint(Point& point,
//                               unsigned int elem_ID,
//                               unsigned int mesh_ID);
  
  
//  /// returns the value of the DOFs of the element. The values are obtained 
//  /// from the vector sol for the interpolation
//  /// @param elem element for which the dofs need to be obtained
//  /// @param vector a vector in which the values will be returned
//  void getDofValuesForSourceElem(Elem* elem, 
//                                 DenseVector<double>& vector);
  
//  /// method to interpolate and create the loads for next analysis
//  void interpolateAndCreateLoads(const std::vector<unsigned int>& load_cases);
  
  
//  /// writes loads from the interpolated values
//    void createLoads(const std::vector<unsigned int>& load_cases);
  
  
  /// FESystemController object that owns this object
  FESystem::FESystemController& fesystem_controller;
  
  /// InterpolationCase for this driver
  const InterpolationCase& interpolation_case;
  
//  /// mesh ID
//  unsigned int mesh_ID;
//  
//  /// mesh to which data will be interpolated
//  MeshDS::FEMesh* mesh;
//  
//  /// mesh data to which data will be interpolated
//  MeshDS::FEMeshData* mesh_data;
//  
//  /// dofmap of the mesh to which the data will be interpolated
//  DofMap* dof_map;
//
//  /// reference to the solution to be interpolated
//  NumericVector<double>* sol_to_interpolate;
//  
};


inline
FESystem::FESystemController& 
InterpolationBase::getFESystemController()
{
  return this->fesystem_controller;
}


inline
const InterpolationCase&
InterpolationBase::getInterpolationCase()
{
  return this->interpolation_case;
}



//inline 
//MeshDS::FEMesh* 
//InterpolationBase::getInterpolationMesh()
//{
//  Assert(this->mesh != NULL, ExcInternalError());
//  
//  return this->mesh;
//}
//
//
//
//inline 
//MeshDS::FEMeshData*
//InterpolationBase::getInterpolationMeshData()
//{
//  Assert(this->mesh_data != NULL, ExcInternalError());
//  
//  return this->mesh_data;
//}
//
//
//inline 
//unsigned int
//InterpolationBase::getInterpolationMeshID()
//{
//  Assert(this->mesh_ID != FESystemNumbers::InvalidID, ExcInternalError());
//  
//  return this->mesh_ID;
//}


#endif //__fesystem_interpolation_base_h__
