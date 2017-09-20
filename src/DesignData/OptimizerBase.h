// $Id: OptimizerBase.h,v 1.15.4.4 2007/05/11 05:16:54 manav Exp $

#ifndef __fesystem_optimizer_base_h__
#define __fesystem_optimizer_base_h__

// C++ includes
#include <string>
#include <map>
#include <memory>

// FESystem includes
#include "Utilities/NameEnumHandler.h"

namespace FESystem
{
  class FESystemController; 
}

#ifndef UNCONSTRAINED_MINIMIZATION_ENUM_ID
#define UNCONSTRAINED_MINIMIZATION_ENUM_ID 1
#else
#error
#endif

#ifndef UNCONSTRAINED_MINIMIZATION_ENUM_NAME
#define UNCONSTRAINED_MINIMIZATION_ENUM_NAME "UNCONSTRAINED_MINIMIZATION"
#else
#error
#endif



#ifndef CONSTRAINED_MINIMIZATION_ENUM_ID
#define CONSTRAINED_MINIMIZATION_ENUM_ID 2
#else
#error
#endif

#ifndef CONSTRAINED_MINIMIZATION_ENUM_NAME
#define CONSTRAINED_MINIMIZATION_ENUM_NAME "CONSTRAINED_MINIMIZATION"
#else
#error
#endif


#ifndef COMPLEMENTARITY_ENUM_ID
#define COMPLEMENTARITY_ENUM_ID 3
#else
#error
#endif

#ifndef COMPLEMENTARITY_ENUM_NAME
#define COMPLEMENTARITY_ENUM_NAME "COMPLEMENTARITY"
#else
#error
#endif


#ifndef NEWTONS_LINE_SEARCH_ENUM_ID
#define NEWTONS_LINE_SEARCH_ENUM_ID 1
#else
#error
#endif

#ifndef NEWTONS_LINE_SEARCH_ENUM_NAME
#define NEWTONS_LINE_SEARCH_ENUM_NAME "NEWTONS_LINE_SEARCH"
#else
#error
#endif


#ifndef NEWTONS_TRUST_REGION_ENUM_ID
#define NEWTONS_TRUST_REGION_ENUM_ID 2
#else
#error
#endif

#ifndef NEWTONS_TRUST_REGION_ENUM_NAME
#define NEWTONS_TRUST_REGION_ENUM_NAME "NEWTONS_TRUST_REGION" 
#else
#error
#endif


#ifndef LIMITED_MEMORY_VARIABLE_METRIC_ENUM_ID
#define LIMITED_MEMORY_VARIABLE_METRIC_ENUM_ID 3
#else
#error
#endif

#ifndef LIMITED_MEMORY_VARIABLE_METRIC_ENUM_NAME
#define LIMITED_MEMORY_VARIABLE_METRIC_ENUM_NAME "LIMITED_MEMORY_VARIABLE_METRIC" 
#else
#error
#endif


#ifndef FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID
#define FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID 4
#else
#error
#endif

#ifndef FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME
#define FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME "FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT" 
#else
#error
#endif


#ifndef POLAK_RIBIERE_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID
#define POLAK_RIBIERE_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID 5
#else
#error
#endif

#ifndef POLAK_RIBIERE_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME
#define POLAK_RIBIERE_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME "POLAK_RIBIERE_NONLINEAR_CONJUGATE_GRADIENT" 
#else
#error
#endif


#ifndef POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID
#define POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID 6
#else
#error
#endif

#ifndef POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME
#define POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME "POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT" 
#else
#error
#endif



#ifndef QUADRATIC_NEWTON_TRUST_REGION_ENUM_ID
#define QUADRATIC_NEWTON_TRUST_REGION_ENUM_ID 7
#else
#error
#endif

#ifndef QUADRATIC_NEWTON_TRUST_REGION_ENUM_NAME
#define QUADRATIC_NEWTON_TRUST_REGION_ENUM_NAME "QUADRATIC_NEWTON_TRUST_REGION" 
#else
#error
#endif


#ifndef QUADRATIC_INTERIOR_POINT_ENUM_ID
#define QUADRATIC_INTERIOR_POINT_ENUM_ID 8
#else
#error
#endif

#ifndef QUADRATIC_INTERIOR_POINT_ENUM_NAME
#define QUADRATIC_INTERIOR_POINT_ENUM_NAME "QUADRATIC_INTERIOR_POINT" 
#else
#error
#endif



#ifndef KUNH_TUCKER_ENUM_ID
#define KUNH_TUCKER_ENUM_ID 9
#else
#error
#endif

#ifndef KUNH_TUCKER_ENUM_NAME
#define KUNH_TUCKER_ENUM_NAME "KUNH_TUCKER" 
#else
#error
#endif


#ifndef INFEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_ID
#define INFEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_ID 9
#else
#error
#endif

#ifndef INFEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_NAME
#define INFEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_NAME "INFEASIBLE_SEMISMOOTH_LINESEARCH" 
#else
#error
#endif


#ifndef FEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_ID
#define FEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_ID 10
#else
#error
#endif

#ifndef FEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_NAME
#define FEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_NAME "FEASIBLE_SEMISMOOTH_LINESEARCH" 
#else
#error
#endif






namespace DesignData
{
  // forward declerations
  class OptimizerInfo; 

  DeclareEnumClass(OptimizationProblemClassEnum);
  DeclareEnumClass(OptimizationAlgorithmKindEnum);
  
  
  DeclareEnumName(UNCONSTRAINED_MINIMIZATION, OptimizationProblemClassEnum,
                  UNCONSTRAINED_MINIMIZATION_ENUM_ID,
                  UNCONSTRAINED_MINIMIZATION_ENUM_NAME);
  
  DeclareEnumName(CONSTRAINED_MINIMIZATION, OptimizationProblemClassEnum,
                  CONSTRAINED_MINIMIZATION_ENUM_ID,
                  CONSTRAINED_MINIMIZATION_ENUM_NAME);
  
  DeclareEnumName(COMPLEMENTARITY, OptimizationProblemClassEnum,
                  COMPLEMENTARITY_ENUM_ID,
                  COMPLEMENTARITY_ENUM_NAME);
  

  DeclareEnumName(NEWTONS_LINE_SEARCH, OptimizationAlgorithmKindEnum,
                  NEWTONS_LINE_SEARCH_ENUM_ID,
                  NEWTONS_LINE_SEARCH_ENUM_NAME);


  DeclareEnumName(NEWTONS_TRUST_REGION, OptimizationAlgorithmKindEnum,
                  NEWTONS_TRUST_REGION_ENUM_ID,
                  NEWTONS_TRUST_REGION_ENUM_NAME);
  
  
  DeclareEnumName(LIMITED_MEMORY_VARIABLE_METRIC, OptimizationAlgorithmKindEnum,
                  LIMITED_MEMORY_VARIABLE_METRIC_ENUM_ID,
                  LIMITED_MEMORY_VARIABLE_METRIC_ENUM_NAME);
  
  
  DeclareEnumName(FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT, OptimizationAlgorithmKindEnum,
                  FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID,
                  FLETCHER_REEVES_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME);
  
  
  DeclareEnumName(POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT, OptimizationAlgorithmKindEnum,
                  POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT_ENUM_ID,
                  POLAK_RIBIERE_PLUS_NONLINEAR_CONJUGATE_GRADIENT_ENUM_NAME);
  
  
  DeclareEnumName(QUADRATIC_NEWTON_TRUST_REGION, OptimizationAlgorithmKindEnum,
                  QUADRATIC_NEWTON_TRUST_REGION_ENUM_ID,
                  QUADRATIC_NEWTON_TRUST_REGION_ENUM_NAME);
  
  
  DeclareEnumName(QUADRATIC_INTERIOR_POINT, OptimizationAlgorithmKindEnum,
                  QUADRATIC_INTERIOR_POINT_ENUM_ID,
                  QUADRATIC_INTERIOR_POINT_ENUM_NAME);
  
  
  DeclareEnumName(KUNH_TUCKER, OptimizationAlgorithmKindEnum,
                  KUNH_TUCKER_ENUM_ID,
                  KUNH_TUCKER_ENUM_NAME);
  
  
  DeclareEnumName(INFEASIBLE_SEMISMOOTH_LINESEARCH, OptimizationAlgorithmKindEnum,
                  INFEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_ID,
                  INFEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_NAME);
  
  
  DeclareEnumName(FEASIBLE_SEMISMOOTH_LINESEARCH, OptimizationAlgorithmKindEnum,
                  FEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_ID,
                  FEASIBLE_SEMISMOOTH_LINESEARCH_ENUM_NAME);

  
  /// this is a base class for the various optimizers classes that are created
  class OptimizerBase
    {
    public:	
      
      /// creator function will need a pointer to the AnalysisDriver that owns the 
      /// instantiation of this class.
      OptimizerBase(const DesignData::OptimizerInfo& info);
      
      /// destructor fuction
      virtual ~OptimizerBase() = 0;
      
      // clears the data structures
      virtual void clear();
      
      /// @returns the solver type
      inline unsigned int getOptimizerClassEnumID() const;
      
      /// @returns the optimizer class class enum name
      inline std::string getOptimizerClassEnumName() const;
      
      /// attaches the fesystem controller object to the 
      void attachFESystemController(FESystem::FESystemController& controller);
      
      /// @returns a reference to the fesystem controller object 
      FESystem::FESystemController& getFESystemController();
      
      /// starts the optimization procedure
      virtual void optimize() = 0;
      
      
    protected:	
      
      /// monitor routine to update the solution vector after each interation
      void iterationMonitor();  
      
      
      /// queries the fesystem_controller for the objective function evaluation and its sensitivities
      void evaluateObjectiveFunctionValueAndGradients();
      
      
      /// queries the fesystem_controller for the constraint function evaluations and their sensitivities
      void evaluateConstraintFunctionValuesAndGradients();
      
      
      /// fesystem controller pointer
      FESystem::FESystemController* fesytem_controller;
      
      unsigned int optimizer_class_enum_ID;
    };
  

  /// @returns an optimizer object for the given details in the OptimizerInfo object
  std::auto_ptr<DesignData::OptimizerBase> 
  DesignData::getOptimizer(const DesignData::OptimizerInfo &info);
}



inline
unsigned int 
Optimizer::FESystemOptimizerBase::getOptimizerClassEnumID() const
{
  return this->optimizer_class_enum_ID;
}





inline 
std::string
Optimizer::FESystemOptimizerBase::getOptimizerClassEnumName() const
{
  return Optimizer::OptimizerClassEnum::enumName(this->optimizer_class_enum_ID);
}





#endif // __fesystem_optimizer_base_h__
