// $Id: ThermalSolution.C,v 1.1.6.2 2007-05-08 05:19:13 manav Exp $

// FESystem includes
#include "Solutions/ThermalSolution.h"
#include "FESystem/AnalysisCase.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/DisciplineInfo.h"
#include "Utilities/InputOutputUtility.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignDatabase.h"
#include "Database/DataInfo.h"

Solution::ThermalSolution::ThermalSolution():
Solution::SolutionBase(Solution::THERMAL_SOLUTION::num())
{
  // add the participating discipline to this solution
  this->addParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num());
}



Solution::ThermalSolution::~ThermalSolution()
{
  
}



unsigned int 
Solution::ThermalSolution::getDisciplineTransientNatureEnumID(const unsigned int discipline) const 
{
  unsigned int enum_ID = 0;

  switch (discipline)
    {
    case THERMAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::STEADY_STATE_SOLUTION::num();
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return enum_ID;
}



void
Solution::ThermalSolution::solve()
{
  const Discipline::DisciplineInfo& info = 
    this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
    (Discipline::THERMAL_DISCIPLINE::num());
  
  switch (info.getAnalysisTypeEnumID())
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
}




std::istream&
Solution::ThermalSolution::readFromInputStream(std::istream& input)
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
Solution::ThermalSolution::getTimeIndependentNodalSolutionDataInfo
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

	    data_info->setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
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
	      
	      data_info->setDisciplineEnumID(Discipline::THERMAL_DISCIPLINE::num());
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
Solution::ThermalSolution::getTimeDependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // unused parameter
  (void) discipline_enum_ID;

  // this solution does not have any transients. Hence, return an empty vector
  FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet> vec;
  return vec;
}

