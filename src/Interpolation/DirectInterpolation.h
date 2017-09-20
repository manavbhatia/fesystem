// $Id: DirectInterpolation.h,v 1.1.2.2 2007-06-13 14:57:20 manav Exp $

#ifndef __fesystem_direct_interpolation_h__
#define __fesystem_direct_interpolation_h__


// C++ includes


// FESystem includes
#include "Interpolation/InterpolationBase.h"

// libMesh includes


// Forward declerations



/// this transfers a specified variable from one mesh to the other. This assumes that
/// the two mesh have the same topology, hence, no interpolation is needed, as the 
/// data from the nodes of one finite element mesh can be consistently applied to nodes 
/// of the other mesh
class DirectInterpolation: public InterpolationBase
{
public:
  /// constructor
  DirectInterpolation(FESystem::FESystemController& controller,
                      const InterpolationCase& interp_case);
  
  /// destructor 
  virtual ~DirectInterpolation();
  
  /// interpolates the laods for all load cases and DVs
  virtual void interpolate();
  
protected:
    
//  /// initialization function
//  void init();
  
  /// writes loads from the interpolated values
  void createLoads();
};


#endif //__fesystem_interpolation_base_h__
