// $Id:$
/*
 *  NonlinearStructuralAnalysisDriver.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 11/22/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

#ifndef __fesystem_updated_lagrangian_nonlinear_structural_analysis_driver_h__ 
#define __fesystem_updated_lagrangian_nonlinear_structural_analysis_driver_h__ 

// C++ includes


// FESystem includes
#include "AnalysisDriver/TransientAnalysisDriver.h"



namespace Driver {
  
  /// This class derives from the \p NonlinearTransientAnalysisDriver class and (re)implements some 
  /// methods to suit well for use for a nonlinear structural analysis. Nonlinear structural analysis
  /// uses a time parameter to define the load steps for a static solution. The transient analysis with 
  /// inertial forces is also dependent on the same approach. 
  
  class UpdatedLagrangianNonlinearStructuralAnalysisDriver: protected NonlinearTransientAnalysisDriver
    {
    public: 
      /// constructor 
      UpdatedLagrangianNonlinearStructuralAnalysisDriver(const unsigned int ID,
                                        FESystem::FESystemController& controller);
      
      /// destructor
      virtual ~UpdatedLagrangianNonlinearStructuralAnalysisDriver();
      
      
      /// this method reimplements the solution method, since it requires updating the mesh at every 
      /// iteration
      virtual void setNewIteration(const double time,
                                   const unsigned int iter_num,
                                   const std::vector<NumericVector<double>*>& sol_vec);
      
    protected: 
      
      
    };
  
}



#endif // __fesystem_updated_lagrangian_nonlinear_structural_analysis_driver_h__

