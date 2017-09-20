// $Id: AnalysisCase.C,v 1.8.4.2 2007-06-13 14:55:12 manav Exp $


// FESystem includes
#include "FESystem/AnalysisCase.h"
#include "Utilities/InputOutputUtility.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "OutputProcessor/OutputInfo.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Interpolation/InterpolationCase.h"
#include "Solutions/SolutionBase.h"

FESystem::AnalysisCase::AnalysisCase():
analysis_title (),
output_info(NULL)
{
  
}






FESystem::AnalysisCase::~AnalysisCase()
{
  // delete the info data structs
  if (this->output_info != NULL) 
    {
    delete this->output_info;
    this->output_info = NULL;
    }
  
//  if (this->interpolation_case != NULL)
//    {
//    delete this->interpolation_case;
//    this->interpolation_case = NULL;
//    }
  
  {
    // and also the solver info in the map
    FESystem::AnalysisCase::IDToSolverInfoMap::iterator it, end;
    it = this->ID_to_solver_info_map.begin();
    end = this->ID_to_solver_info_map.end();
    
    for (; it != end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
    this->ID_to_solver_info_map.clear();
  }
  
  {
    // and also the solver info in the map
    FESystem::AnalysisCase::DisciplineInfoMap::iterator it, end;
    it = this->discipline_info_map.begin();
    end = this->discipline_info_map.end();
    
    for (; it != end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
    this->discipline_info_map.clear();
  }

  {
    // and also the solver info in the map
    FESystem::AnalysisCase::SolutionBaseMap::iterator it, end;
    it = this->solution_base_map.begin();
    end = this->solution_base_map.end();
    
    for (; it != end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
    this->solution_base_map.clear();
  }
  
}


std::vector<Discipline::DisciplineInfo*>& 
FESystem::AnalysisCase::getAnalysisDisciplineInfos()
{
  static std::vector<Discipline::DisciplineInfo*> infos;
  
  if (infos.size() == 0)
  {
  FESystem::AnalysisCase::DisciplineInfoMap::iterator it, end;
  it = this->discipline_info_map.begin();
  end = this->discipline_info_map.end();
  
  for ( ; it != end; it++)
    infos.push_back(it->second);
  }
  
  return infos;
}



std::vector<Solution::SolutionBase*> 
FESystem::AnalysisCase::getSolutions()
{
  std::vector<Solution::SolutionBase*> vec;
  
  FESystem::AnalysisCase::SolutionBaseMap::iterator it, end;
  it = this->solution_base_map.begin();
  end = this->solution_base_map.end();
  
  for (; it != end; it++)
    vec.push_back(it->second);
  
  return vec;
}



std::istream& 
FESystem::AnalysisCase::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0;
  unsigned int enum_ID = 0;
  bool insert_success = false;

  // this will begin with BEGIN_ANALYSIS_CASE
  FESystemIOUtility::readFromInput(input, "ANALYSIS_CASE");	
  FESystemIOUtility::readFromInput(input, "BEGIN");	
	
	
  // read  in the title
  FESystemIOUtility::readFromInput(input, "ANALYSIS_TITLE", this->analysis_title);	
  
  // now read the analysis discipline
  FESystemIOUtility::readFromInput(input, "N_ANALYSIS_DISCIPLINE_INFO", num);
    
  for (unsigned int i = 0; i < num; i++)
    {
    tag.clear(); 
    FESystemIOUtility::peekFromInput(input, tag);
    
    enum_ID = Discipline::DisciplineInfoTypeEnum::enumID(tag);
    
    Discipline::DisciplineInfo *info =
      Discipline::createDisciplineInfo(enum_ID).release();
    
    input >> (*info);
    
    this->addAnalysisDisciplineInfo(info);
    }
  
  
  // now read the solutions to be performed
  FESystemIOUtility::readFromInput(input, "N_SOLUTIONS", num);
  
  for (unsigned int i = 0; i < num; i++)
    {
    tag.clear(); 
    FESystemIOUtility::peekFromInput(input, tag);
    
    enum_ID = Solution::SolutionEnum::enumID(tag);
    
    Solution::SolutionBase *solution =
      Solution::createSolutionBase(enum_ID).release();
    
    input >> (*solution);
    
    this->addSolutionBase(solution);
    }

  
  
  
  // read the solver info
  Solver::SolverInfo* solver_info = NULL;
  FESystemIOUtility::readFromInput(input, "SOLVER_INFO");	
  FESystemIOUtility::readFromInput(input, "BEGIN");	
  
  FESystemIOUtility::readFromInput(input, "N_SOLVER_INFO", num);
  
  for (unsigned int i=0; i < num; i++)
    {
    tag.clear();
    FESystemIOUtility::peekFromInput(input, tag);
    
    enum_ID = Solver::SolverInfoEnum::enumID(tag);
    
    solver_info = Solver::createSolverInfo(enum_ID).release();

    input >> (*solver_info);

    this->addSolverInfo(solver_info);
    }
  
  FESystemIOUtility::readFromInput(input, "SOLVER_INFO");	
  FESystemIOUtility::readFromInput(input, "END");	
  
  
  
  // read the interpolation cases 
  InterpolationCase* interpol_case = NULL;
  FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASES");	
  FESystemIOUtility::readFromInput(input, "BEGIN");	
  
  FESystemIOUtility::readFromInput(input, "N_INTERPOLATION_CASES", num);
  
  for (unsigned int i=0; i < num; i++)
    {
    interpol_case = new InterpolationCase;
    
    interpol_case->readFromInputStream(input);
    
    this->addInterpolationCase(interpol_case);
    }
  
  FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASES");	
  FESystemIOUtility::readFromInput(input, "END");

  
  // read solver parameters
  // first read the integer parameters
  FESystemIOUtility::readFromInput(input, "INTEGER_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  unsigned int n_params = 0, int_param=0;
  
  FESystemIOUtility::readFromInput(input, "N_INTEGER_PARAMETERS", n_params);
  
  for (unsigned int i=0; i < n_params; i++)
    {
    tag.clear();
    input >> tag;
    input >> int_param;
		
    insert_success =
      this->integer_params.insert(AnalysisCase::IntegerParameterMap::value_type
                                  (tag, int_param)).second;
    Assert(insert_success,
           AnalysisCase::ExcDuplicateParameterName(tag));
    }
  
  n_params = 0;
  
  FESystemIOUtility::readFromInput(input, "INTEGER_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "END");
  
	
  // first read the real parameters
  FESystemIOUtility::readFromInput(input, "REAL_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  double real_param=0;
  
  FESystemIOUtility::readFromInput(input, "N_REAL_PARAMETERS", n_params);
  
  for (unsigned int i=0; i < n_params; i++)
    {
    tag.clear();
    input >> tag;
    input >> real_param;
		
    insert_success =
      this->float_params.insert(AnalysisCase::RealParameterMap::value_type
                                (tag, real_param)).second;
    Assert(insert_success,
           AnalysisCase::ExcDuplicateParameterName(tag));
    }
  
  n_params = 0;
  
  FESystemIOUtility::readFromInput(input, "REAL_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "END");
	
  
  // first read the string parameters
  FESystemIOUtility::readFromInput(input, "STRING_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  std::string str_param;
  
  FESystemIOUtility::readFromInput(input, "N_STRING_PARAMETERS", n_params);
  
  for (unsigned int i=0; i < n_params; i++)
    {
    tag.clear();
    input >> tag;
    input >> str_param;
		
    insert_success =
      this->string_params.insert(AnalysisCase::StringParameterMap::value_type
                                 (tag, str_param)).second;
    Assert(insert_success,
           AnalysisCase::ExcDuplicateParameterName(tag));
    }
  
  n_params = 0;
  
  FESystemIOUtility::readFromInput(input, "STRING_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "END");
  
	
  // read in the output info
  this->output_info = new OutputInfo;
  input >> *(this->output_info);
  
	
  // this will end with END_ANALYSIS_CASE
  FESystemIOUtility::readFromInput(input, "ANALYSIS_CASE");
  FESystemIOUtility::readFromInput(input, "END");
  
		
  return input;
}




void
FESystem::AnalysisCase::addAnalysisDisciplineInfo(Discipline::DisciplineInfo* info)
{
  unsigned int ID = info->getDisciplineEnumID();
  
  bool insert_success = this->discipline_info_map.insert
    (FESystem::AnalysisCase::DisciplineInfoMap::value_type(ID, info)).second;
  
  Assert(insert_success, FESystemExceptions::ExcDuplicateID("Discipline", ID));
  
}



void
FESystem::AnalysisCase::addInterpolationCase(InterpolationCase* ipcase)
{
  unsigned int ID = ipcase->ID();
  
  bool insert_success = this->interpolation_case_map.insert
    (std::map<unsigned int, InterpolationCase*>::value_type(ID, ipcase)).second;
  
  Assert(insert_success, FESystemExceptions::ExcDuplicateID("InterpolationCase", ID));
  
}



void
FESystem::AnalysisCase::addSolverInfo(Solver::SolverInfo* info)
{
  unsigned int ID = info->getID();
  
  bool insert_success = this->ID_to_solver_info_map.insert
    (FESystem::AnalysisCase::IDToSolverInfoMap::value_type(ID, info)).second;
  
  Assert(insert_success, FESystemExceptions::ExcDuplicateID("Solver", ID));
  
}


void
FESystem::AnalysisCase::addSolutionBase(Solution::SolutionBase* info)
{
  unsigned int ID = info->getID();
  
  bool insert_success = this->solution_base_map.insert
    (FESystem::AnalysisCase::SolutionBaseMap::value_type(ID, info)).second;
  
  Assert(insert_success, FESystemExceptions::ExcDuplicateID("Solution", ID));
  
}
