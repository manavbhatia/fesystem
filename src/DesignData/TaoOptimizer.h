// $Id:$

#ifndef __fesystem_tao_optimizer_h__
#define __fesystem_tao_optimizer_h__

// FESystem includes
#include "DesignData/OptimizerBase.h"


namespace DesignData
{
  
  
  class TaoOptimizer
    {
    public:
      TaoOptimizer(const DesignData::OptimizerDataInfo& data_info);
      
      ~TaoOptimizer();
      
      virtual void optimize();
      
      
    protected:
      
      TAO_APPLICATION tao_application; 
      
      TAO_SOLVER tao_solver; 
      
      TaoMethod method;
      
    }
}






#endif // __fesystem_tao_optimizer_h__ 
