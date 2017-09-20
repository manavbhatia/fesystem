// $Id:$
/*
 *  NonlinearStructuralAnalysisDriver.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 11/22/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

// C++ includes


// FESystem includes
#include "AnalysisDriver/UpdatedLagrangianNonlinearStructuralAnalysisDriver.h"
#include "Discipline/StructuralAnalysis.h"
#include "Utilities/Log.h"


Driver::UpdatedLagrangianNonlinearStructuralAnalysisDriver::
UpdatedLagrangianNonlinearStructuralAnalysisDriver(const unsigned int ID,
                                                   FESystem::FESystemController& controller):
Driver::NonlinearTransientAnalysisDriver(ID, controller)
{
  
}


Driver::UpdatedLagrangianNonlinearStructuralAnalysisDriver::
~UpdatedLagrangianNonlinearStructuralAnalysisDriver()
{
  
}



void 
Driver::UpdatedLagrangianNonlinearStructuralAnalysisDriver::
setNewIteration(const double time,
                const unsigned int iter_num,
                const std::vector<NumericVector<double>*>& sol_vec)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  // make sure that the numbers make some sense
  Assert(time == this->getSimulatedTime() && iter_num == this->getSimulatedIterationNumber(), 
         ExcInternalError());

  Discipline::StructuralAnalysis& discipline = dynamic_cast<Discipline::StructuralAnalysis&> 
  (this->getAnalysisDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num()));
  
  // store the iteration time and update the mesh using the displacement solution
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      this->time_values.push_back(time);
      // get the displacement from the last iteration
      // get the solution vectors from the database
      NumericVector<double>& scratch_vec1 = this->getVector("ScratchVector1");
      FESystemDatabase::TimeDependentDataInfo data_info;
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(0);
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo((this->getCurrentIterationNumber()-1), 
                                          this->getCurrentAnalysisTime());
      
      this->fesystem_controller.global_data_storage->fillVector(data_info, scratch_vec1);
      scratch_vec1.close();
      scratch_vec1.scale(-1.0);
      scratch_vec1.add(1.0, *(sol_vec[0]));
      scratch_vec1.close();
      
      discipline.setMeshNodalIncrementsFromDisplacementSolution(scratch_vec1);
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      Assert(time == this->time_values[iter_num-1], ExcInternalError());
      
      // get the solution vectors from the database
      NumericVector<double>& scratch_vec1 = this->getVector("ScratchVector1");
      NumericVector<double>& scratch_vec2 = this->getVector("ScratchVector2");
      FESystemDatabase::TimeDependentDataInfo data_info;
      
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(0);
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo(this->getCurrentIterationNumber(), 
                                          this->getCurrentAnalysisTime());      
      this->fesystem_controller.global_data_storage->fillVector(data_info, scratch_vec1);
      scratch_vec1.close();
      
      data_info.clear();
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(0);
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo((this->getCurrentIterationNumber()-1), 
                                          this->getCurrentAnalysisTime());      
      this->fesystem_controller.global_data_storage->fillVector(data_info, scratch_vec2);
      scratch_vec2.close();
      scratch_vec1.add(-1.0, scratch_vec2);
      scratch_vec1.close();
      
      // update the mesh nodal coordinates
      discipline.setMeshNodalIncrementsFromDisplacementSolution(scratch_vec1);
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
    
  // make sure that the number of solutions in the vector is same as the order of the
  // system
  unsigned int order = it->second->getTransientSystemOrder();
  
  Assert(sol_vec.size() == (order+1),
         ExcInternalError());
  
  FESystemDatabase::TimeDependentDataInfo data_info;
  
  // now store the solution vector
  for (unsigned int i=0; i <= order; i++)
    {
      data_info.clear();
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(i);
      if (this->current_analysis_kind_enum_ID == SENSITIVITY_ANALYSIS::num())
        data_info.setDVID(this->getCurrentDesignParameter().getID());	
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo(iter_num, time);
      
      this->fesystem_controller.global_data_storage->storeVector(data_info, *(sol_vec[i]));
      this->fesystem_controller.log_file->write2(data_info);
      this->writeSolutionVectorToLog(*(sol_vec[i]), it->second->getDisciplineEnumID());
    }
}


