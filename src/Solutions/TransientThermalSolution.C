// $Id: TransientThermalSolution.C,v 1.1.2.1 2007-05-15 20:38:52 manav Exp $

// FESystem includes
#include "Solutions/TransientThermalSolution.h"
#include "FESystem/AnalysisCase.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/DisciplineInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignDatabase.h"
#include "Database/DataInfo.h"
#include "Numerics/PetscSeqVector.h"

Solution::TransientThermalSolution::TransientThermalSolution():
Solution::SolutionBase(Solution::TRANSIENT_THERMAL_SOLUTION::num())
{
  // add the participating discipline to this solution
  this->addParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num());
}



Solution::TransientThermalSolution::~TransientThermalSolution()
{
  
}



unsigned int 
Solution::TransientThermalSolution::getDisciplineTransientNatureEnumID
(const unsigned int discipline) const 
{
  unsigned int enum_ID = 0;
  
  switch (discipline)
    {
    case THERMAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::TRANSIENT_SOLUTION::num();
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return enum_ID;
}



void
Solution::TransientThermalSolution::solve()
{
  const Discipline::DisciplineInfo& info = 
  this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
  (Discipline::THERMAL_DISCIPLINE::num());
  
  switch (info.getAnalysisTypeEnumID())
    {
    case LINEAR_ANALYSIS_ENUM_ID:
      this->linearTransientSolve(Discipline::THERMAL_DISCIPLINE::num());
      break;
      
    case NONLINEAR_ANALYSIS_ENUM_ID:
      this->nonlinearTransientSolve(Discipline::THERMAL_DISCIPLINE::num());
      break;
      
    default:
      Assert(false, ExcInternalError());
      break;
    }
}




std::istream&
Solution::TransientThermalSolution::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0;
  
  FESystemIOUtility::readFromInput(input, Solution::SolutionEnum::enumName(this->solution_enum_ID));
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  
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
Solution::TransientThermalSolution::getTimeIndependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // unused parameter
  (void) discipline_enum_ID;
  
  FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase> data_info_vec;
  
  return data_info_vec;
}




FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
Solution::TransientThermalSolution::getTimeDependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // this solution does not have any transients. Hence, return an empty vector
  FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet> data_info_vec;

  FESystemDatabase::TransientDataInfoSet* data_info_set = NULL;
  
  switch (discipline_enum_ID)
    {
    case THERMAL_DISCIPLINE_ENUM_ID:
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
          data_info.setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
          data_info.setName("TimeValues");
          data_info.setLoadCase(this->load_case_IDs[i]);
          this->fesystem_controller->global_data_storage->fillVector(data_info, vec);
          time_values.resize(vec.size());
          
          for (unsigned int k=0; k < vec.size(); k++)
            time_values[k] = vec.el(k);
          
          data_info_set = new FESystemDatabase::TransientDataInfoSet;

          data_info_set->setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
          data_info_set->setName("Solution");
          data_info_set->setLoadCase(this->load_case_IDs[i]);
          data_info_set->setOrder(0);
          data_info_set->setTransientIterationTimeVector(time_values);
          data_info_vec.push_back(data_info_set);
          data_info_set = NULL;

          for (unsigned int j=0; j < dv_vector->size(); j++)
            {
            data_info_set = new FESystemDatabase::TransientDataInfoSet;
            
            data_info_set->setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
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

