// $Id:$ 

#ifndef __fesystem_optimizer_interface_base_h__
#define __fesystem_optimizer_interface_base_h__


namespace FESystem
{
  // forward declerations
  class FESystem::OptimizerBase;
  class FESystem::FESystemController;
  
  class OptimizationInterface
    {
    public: 
      /// constructor
      OptimizationInterface();
      
      
      /// destructor
      ~OptimizationInterface();
      
      /// sets the optimizer object for this interface
      void setOptimizer(FESystem::OptimizerBase* opt); 
      
      
      /// @returns the objective function value at the design variable value defined in the 
      /// vector
      virtual void getObjectiveAndConstraintFunctionValueAndGradient
      (const FESystem::DesignVariableVector& vec, 
       FESystem::DesignObjectiveFunction& objective,
       FESystem::AutoPtrVec<FESystem::DesignConstraintFunction>& constr_vec) = 0;
      
      
      
    protected: 
      
      /// optimizer object
      FESystem::OptimizerBase* optimizer;
      
    };
}



#endif // __fesystem_fesystem_optimizer_interface_base_h__
