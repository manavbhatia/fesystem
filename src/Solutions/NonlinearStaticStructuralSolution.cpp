// $Id: LinearStressSolution.C,v 1.4.6.3 2007-06-13 14:59:40 manav Exp $

/*
 *  NonlinearStaticStructuralSolution.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 11/21/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

// FESystem includes
#include "Solutions/NonlinearStaticStructuralSolution.h"
#include "FESystem/AnalysisCase.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/DisciplineInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignDatabase.h"
#include "Database/DataInfo.h"
#include "Interpolation/DirectInterpolation.h"


Solution::NonlinearStaticStructuralSolution::NonlinearStaticStructuralSolution():
Solution::SolutionBase(Solution::NONLINEAR_STATIC_STRUCTURAL_SOLUTION::num()),
coupled_thermoelastic(false),
interpolation_case_ID(FESystemNumbers::InvalidID)
{
  // add the participating discipline to this solution
  this->addParticipatingDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num());
}



Solution::NonlinearStaticStructuralSolution::~NonlinearStaticStructuralSolution()
{
  
}


unsigned int 
Solution::NonlinearStaticStructuralSolution::getDisciplineTransientNatureEnumID(const unsigned int discipline) const 
{
  unsigned int enum_ID = 0;
  
  switch (discipline)
  {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::TRANSIENT_SOLUTION::num();
      break;
      
    case THERMAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::STEADY_STATE_SOLUTION::num();
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  return enum_ID;
}



void
Solution::NonlinearStaticStructuralSolution::solve()
{
  // if a coupled thermal-structural response is requested, then
  // first perform a thermal analysis, import the loads to structural and
  // then perform the structural solution
  if (this->coupled_thermoelastic)
    {
      const Discipline::DisciplineInfo& thermal_info = 
      this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
      (Discipline::THERMAL_DISCIPLINE::num());
      
      switch (thermal_info.getAnalysisTypeEnumID())
      {
        case LINEAR_ANALYSIS_ENUM_ID:
          this->linearSolve(Discipline::THERMAL_DISCIPLINE::num());
          break;
          
        case NONLINEAR_ANALYSIS_ENUM_ID:
          this->nonlinearSolve(Discipline::THERMAL_DISCIPLINE::num());
          break;
          
        default:
          Assert(false, ExcInternalError());
          break;
      }
      
      // get the interpolation case from the analysis case
      const InterpolationCase& interpolation_case = 
      this->fesystem_controller->analysis_case->getInterpolationCase(this->interpolation_case_ID);
      std::auto_ptr<DirectInterpolation> interpolation
      (new DirectInterpolation(*(this->fesystem_controller), interpolation_case));
      interpolation->interpolate();
    }
  
  this->nonlinearTransientSolve(Discipline::STRUCTURAL_DISCIPLINE::num());
}



std::istream&
Solution::NonlinearStaticStructuralSolution::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0;
  
  FESystemIOUtility::readFromInput(input, Solution::SolutionEnum::enumName(this->solution_enum_ID));
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  FESystemIOUtility::readFromInput(input, "COUPLED_THERMAL_STRUCTURAL", 
                                   this->coupled_thermoelastic);
  if (this->coupled_thermoelastic)
    {
      this->addParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num());
      FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASE_ID", 
                                       this->interpolation_case_ID);
    }
  
  // Read in the list of load cases
  unsigned int case_ID =0;
  FESystemIOUtility::readFromInput(input, "N_ANALYSIS_LOAD_CASES", num);
  for (unsigned int case_it =0; case_it < num; case_it++)
    {
      case_ID = 0;
      input >> case_ID;
      
      this->load_case_IDs.push_back(case_ID);
    }
  
  FESystemIOUtility::readFromInput(input, "N_SOLVER_INFO_ID", num);
  unsigned int solver_class = 0, info_id = 0;
  for (unsigned int i=0; i < num; i++)
    {
      tag.clear();
      input >> tag;
      solver_class = Solver::SolverClassEnum::enumID(tag);
      input >> info_id;
      
      this->addSolverID(solver_class, info_id);
    }
  
  FESystemIOUtility::readFromInput(input, Solution::SolutionEnum::enumName(this->solution_enum_ID));
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}




FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
Solution::NonlinearStaticStructuralSolution::getTimeIndependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase> data_info_vec;
  
  switch (discipline_enum_ID)
  {
    case THERMAL_DISCIPLINE_ENUM_ID:
    {
      std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
      dv_vector(this->fesystem_controller->design_database->getParameters().release());
      
      data_info_vec.resize(this->load_case_IDs.size() * (1 + dv_vector->size()));
      
      FESystemDatabase::DataInfoBase *data_info  = NULL;
      
      unsigned int location_id = 0;
      
      // iterate over each load case and DV and add the solutions
      for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
        {
          data_info = new FESystemDatabase::TimeIndependentDataInfo;
          
          data_info->setDisciplineEnumID(discipline_enum_ID);
          data_info->setName("Solution");
          data_info->setLoadCase(this->load_case_IDs[i]);
          
          data_info_vec.reset(location_id, data_info);
          location_id++;
          data_info = NULL;
        }
      
      for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
        for (unsigned int j=0; j < dv_vector->size(); j++)
          {
            data_info = new FESystemDatabase::TimeIndependentDataInfo;
            
            data_info->setDisciplineEnumID(discipline_enum_ID);
            data_info->setName("Solution");
            data_info->setLoadCase(this->load_case_IDs[i]);
            data_info->setDVID((*dv_vector)[j]->getID());
            
            data_info_vec.reset(location_id, data_info);
            location_id++;
            data_info = NULL;
          }
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  return data_info_vec;
}





FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
Solution::NonlinearStaticStructuralSolution::getTimeDependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // this solution does not have any transients. Hence, return an empty vector
  FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet> data_info_vec;
  
  FESystemDatabase::TransientDataInfoSet* data_info_set = NULL;
  
  switch (discipline_enum_ID)
  {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
    {
      std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
      dv_vector(this->fesystem_controller->design_database->getParameters().release());
      
      // get the time values for the solution
      FESystemDatabase::TimeIndependentDataInfo data_info;
      FESystemNumerics::PetscSeqVector<double> vec;
      std::vector<double> time_values;
      
      // iterate over each load case and DV and add the solutions
      for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
        {
          data_info.clear();
          data_info.setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
          data_info.setName("TimeValues");
          data_info.setLoadCase(this->load_case_IDs[i]);
          this->fesystem_controller->global_data_storage->fillVector(data_info, vec);
          time_values.resize(vec.size());
          
          for (unsigned int k=0; k < vec.size(); k++)
            time_values[k] = vec.el(k);
          
          data_info_set = new FESystemDatabase::TransientDataInfoSet;
          
          data_info_set->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
          data_info_set->setName("Solution");
          data_info_set->setLoadCase(this->load_case_IDs[i]);
          data_info_set->setOrder(0);
          data_info_set->setTransientIterationTimeVector(time_values);
          data_info_vec.push_back(data_info_set);
          data_info_set = NULL;
          
          for (unsigned int j=0; j < dv_vector->size(); j++)
            {
              data_info_set = new FESystemDatabase::TransientDataInfoSet;
              
              data_info_set->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
              data_info_set->setName("Solution");
              data_info_set->setLoadCase(this->load_case_IDs[i]);
              data_info_set->setOrder(0);
              data_info_set->setDVID((*dv_vector)[j]->getID());
              data_info_set->setTransientIterationTimeVector(time_values);
              data_info_vec.push_back(data_info_set);
              data_info_set = NULL;
            }
        }
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  return data_info_vec;
}
