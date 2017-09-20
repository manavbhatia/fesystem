// $Id: AnalysisCase.h,v 1.5.4.3 2007-06-13 14:55:12 manav Exp $

#ifndef __fesystem_analysis_case_h__
#define __fesystem_analysis_case_h__

// C++ includes
#include <iostream>
#include <vector>
#include <map>
#include <string>

// FESystem includes
#include "FESystem/FESystemExceptions.h"
#include "Discipline/AnalysisDisciplineBase.h"

// Forward declerations
class OutputInfo;

namespace Solver
{
  class SolverInfo;
}

namespace Discipline
{
  class DisciplineInfo;
}

namespace Solution
{
  class SolutionBase;
}


class InterpolationCase;


namespace FESystem
{
  class AnalysisCase
  {
public:
    /// constructor
    AnalysisCase();
    
    /// destructor
    ~AnalysisCase();
        
    /// @returns the title for analysis
    inline const std::string& getTitle() const;

    /// @returns the analysis discipline
    inline const Discipline::DisciplineInfo& getAnalysisDisciplineInfo
      (const unsigned int discipline) const;
    
    /// @returns all analysis disciplines
    std::vector<Discipline::DisciplineInfo*>& getAnalysisDisciplineInfos();
    
    /// @returns a vector containing the solutions for analysis
    std::vector<Solution::SolutionBase*> getSolutions();
    
    /// @returns the real parameter with the given name
    /// @param name name of the parameter
    inline double getRealParameter(const std::string& name) const;
    
    /// @returns the integer parameter with the given name
    /// @param name name of the parameter
    inline unsigned int getIntegerParameter(const std::string& name) const;
    
    /// @returns the file name
    /// @param name name of the parameter for which file name is to be obtained
    inline std::string getStringParameter(const std::string& name) const;
    
    /// reads from the input stream
    std::istream& readFromInputStream(std::istream& input);
	  
    /// @returns reference to the output info data structure
    inline const OutputInfo& getOutputInfo() const;
    
    
    /// @returns the solver info for the given ID
    inline const Solver::SolverInfo& getSolverInfo(const unsigned int ID) const;
    
    /// @returns the interpolation case referenced by the ID
    const InterpolationCase& getInterpolationCase(const unsigned int case_ID) const;
    
protected:

      /// adds the InterpolationCase to this object
      void addInterpolationCase(InterpolationCase* ipcase);

    /// adds the SolverInfo to this object
      void addSolverInfo(Solver::SolverInfo* info);

    /// adds the AnalysisDiscipline to this object
      void addAnalysisDisciplineInfo(Discipline::DisciplineInfo* info);
    
    /// adds the SolutionBase to this object
    void addSolutionBase(Solution::SolutionBase* info);
    
    DeclException1(ExcParameterNameDoesNotExist, std::string&, 
                   << "Parameter Name: " << arg1 << "does not exist.");
      
    DeclException1(ExcDuplicateParameterName, std::string&, 
                   << "Parameter Name: " << arg1 << "already exists.");

    /// title of the analysis
      std::string analysis_title;
    
    /// output format
    OutputInfo *output_info;
    
    /// interpolation case
    std::map<unsigned int, InterpolationCase*> interpolation_case_map;
        
    // local type definitions
    typedef std::map<unsigned int, Solver::SolverInfo*> IDToSolverInfoMap;
    typedef std::map<std::string, unsigned int> IntegerParameterMap;
    typedef std::map<std::string, double> RealParameterMap;
    typedef std::map<std::string, std::string> StringParameterMap;
    typedef std::map<unsigned int, Discipline::DisciplineInfo*> DisciplineInfoMap;
    typedef std::map<unsigned int, Solution::SolutionBase*> SolutionBaseMap;

    /// map of discipline info for each discipline
    DisciplineInfoMap discipline_info_map;
    
    
    /// map of ID to solver info
    IDToSolverInfoMap ID_to_solver_info_map;
    
    /// map of solutions
    SolutionBaseMap solution_base_map;
    
    /// integer parameter map
    IntegerParameterMap integer_params;
    
    /// real parameter map
    RealParameterMap float_params;
    
    /// string parameters
    StringParameterMap string_params;
  };
}






inline
const std::string& 
FESystem::AnalysisCase::getTitle() const
{
  return this->analysis_title;
}




inline 
const Discipline::DisciplineInfo&
FESystem::AnalysisCase::getAnalysisDisciplineInfo(const unsigned int discipline) const
{
  FESystem::AnalysisCase::DisciplineInfoMap::const_iterator it, end;
  it = this->discipline_info_map.find(discipline);
  end = this->discipline_info_map.end();
  
  Assert(it != end, FESystemExceptions::ExcIDDoesNotExist("Discipline", discipline));
  
  return *(it->second);
}






inline
double 
FESystem::AnalysisCase::getRealParameter(const std::string& name) const
{
  FESystem::AnalysisCase::RealParameterMap::const_iterator 
  param_it = this->float_params.find(name);
	
  Assert(param_it != this->float_params.end(),
         FESystem::AnalysisCase::ExcParameterNameDoesNotExist(name));
	
  return param_it->second;
}




inline
unsigned int
FESystem::AnalysisCase::getIntegerParameter(const std::string& name) const
{
  FESystem::AnalysisCase::IntegerParameterMap::const_iterator 
  param_it = this->integer_params.find(name);
	
  Assert(param_it != this->integer_params.end(),
         FESystem::AnalysisCase::ExcParameterNameDoesNotExist(name));
	
  return param_it->second;
}





inline 
std::string 
FESystem::AnalysisCase::getStringParameter(const std::string& name) const
{
  FESystem::AnalysisCase::StringParameterMap::const_iterator 
  param_it = this->string_params.find(name);
	
  Assert(param_it != this->string_params.end(),
         FESystem::AnalysisCase::ExcParameterNameDoesNotExist(name));
	
  return param_it->second;
}



inline 
const OutputInfo& 
FESystem::AnalysisCase::getOutputInfo() const
{
  Assert(this->output_info != NULL,
         ExcEmptyObject());
  
  return *(this->output_info);
}




inline 
const Solver::SolverInfo& 
FESystem::AnalysisCase::getSolverInfo(const unsigned int ID) const
{
  FESystem::AnalysisCase::IDToSolverInfoMap::const_iterator 
  param_it = this->ID_to_solver_info_map.find(ID);
	
  Assert(param_it != this->ID_to_solver_info_map.end(),
         ExcInternalError());
  
	Assert(param_it->second != NULL, 
         ExcEmptyObject());
  
  return *(param_it->second);
}



inline 
const InterpolationCase&  
FESystem::AnalysisCase::getInterpolationCase(const unsigned int case_ID) const
{ 
  std::map<unsigned int, InterpolationCase*>::const_iterator it, end;
  it = this->interpolation_case_map.find(case_ID);
  end = this->interpolation_case_map.end();
  
  AssertThrow(it != end, ExcInternalError());
  
  return *(it->second); 
} 



#endif // __fesystem_analysis_case_h__
