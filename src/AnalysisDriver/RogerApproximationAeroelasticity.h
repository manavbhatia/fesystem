/*
 *  RogerApproximationAeroelasticity.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/20/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

#ifndef __fesystem_roger_approximation_aeroelasticity_driver_h__
#define __fesystem_roger_approximation_aeroelasticity_driver_h__

// C++ includes
#include <vector>

// FESystem incudes
#include "AnalysisDriver/AeroelasticityAnalysisDriverBase.h"


// libmesh includes
#include "numerics/petsc_matrix.h"
#include "petscis.h"


namespace FESystem {
  class FESystemController;
}



#ifndef ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_ID
#define ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_ID 6
#else
#error
#endif

#ifndef ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_NAME
#define ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_NAME "ROGER_APPROXIMATION_AEROELASTICITY_DRIVER"
#else
#error
#endif


namespace Driver{
  
  DeclareEnumName(ROGER_APPROXIMATION_AEROELASTICITY_DRIVER, Driver::AnalysisDriverTypeEnum,
                  ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_ID,
                  ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_NAME);

  
  class RogerApproximationAeroelasticityDriver: public AeroelasticityAnalysisDriverBase
    {
    public:
      
      ///  constructor
      RogerApproximationAeroelasticityDriver(const unsigned int ID,
                                       FESystem::FESystemController& controller);
      
      
      /// destructor
      virtual ~RogerApproximationAeroelasticityDriver();
      
      
    protected:
      
      /// this method initializes itself before the solutions are started for the driver
      virtual void initialize();

      /// this method solves for the current load case
      virtual void solveCurrentLoadCase();
      
      /// function will add the matrices and vectors for this analysis
      virtual void addMatricesAndVectors();
      
      
      void getSubmatrixAfterBoundaryCondition(SparseMatrix<double>& matrix,
                                              SparseMatrix<double>& matrix_sub,
                                              IS row_matrix_index_set, 
                                              IS col_matrix_index_set);

      
      void copyVectorToGlobalVector(NumericVector<double>& sub_vector,
                                    NumericVector<double>& global_vector,
                                    const std::vector<int>& local_to_global_map);

      
      void createSubMatrixIndices (Discipline::AnalysisDisciplineBase& discipline,
                                   IS* row_matrix_index_set, IS* col_matrix_index_set,
                                   std::vector<int>& local_to_global_map);
      
      void destroySubMatrixIndices(IS row_matrix_index_set, IS col_matrix_index_set);
      
    };
  
}


#endif //__fesystem_roger_approximation_aeroelasticity_driver_h__
