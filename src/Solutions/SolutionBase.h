// $Id: SolutionBase.h,v 1.7.6.3 2008-08-21 20:42:16 manav Exp $

#ifndef __fesystem_solution_base_h__
#define __fesystem_solution_base_h__


// C++ includes
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <set>

// FESystem includes
#include "Utilities/NameEnumHandler.h"
#include "Utilities/AutoptrVector.h"
#include "FESystem/FESystemNumbers.h"

// forward declerations
namespace FESystem
{
  class FESystemController;
}


namespace FESystemDatabase
{
  class DataInfoBase;
  class TimeDependentDataInfo;
  class TransientDataInfoSet;
}

namespace Solution
{

  enum SolType
    {
      STATIC,
      EIGEN
    };
  
  
  DeclareEnumClass(SolutionEnum);

  
  /// this class provides the base functionality for all solutions
  class SolutionBase
  {
public:
    /// constructor
    SolutionBase(const unsigned int sol_enum_ID);
    
    /// destructor
    virtual ~SolutionBase();
    
    /// @returns the ID of this object
    inline unsigned int getID() const;
    
    /// @returns the solution enum ID for this solution
    unsigned int getSolutionEnumID() const;
    
    /// @returns the solution enum name for this solution
    std::string getSolutionEnumName() const;
    
    /// @returns true if the specified discipline participates in the solution
    bool checkParticipatingDiscipline(const unsigned int discipline_enum_id);

    /// @returns the transient nature of the given discipline
    virtual unsigned int getDisciplineTransientNatureEnumID(const unsigned int discipline) const = 0;
    

    /// @returns the vector of SolNameInfo structures for the nodal solutions, and for the
    /// specified disciplines
    virtual FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
      getTimeIndependentNodalSolutionDataInfo(const unsigned int discipline_enum_ID) = 0;
    
    /// @returns the vector of SolNameInfo structures for the nodal solutions, and for the
    /// specified disciplines
    virtual FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
      getTimeDependentNodalSolutionDataInfo(const unsigned int discipline_enum_ID) = 0;
    
    
    /// the vector of analysis load cases
    inline  const std::vector<unsigned int>& getLoadCaseIDs() const;

    
    /// attaches the FESystemController to this object for solution purposes
    void attachFESystemController(FESystem::FESystemController& controller);
          
    
    /// @returns the SolverInfo ID for the solver class 
    /// @params solver_class enumeration ID of the solver class 
    unsigned int getSolverInfoID(const unsigned int solver_class) const;

    /// this method will solve for all the load cases and sensitivity analysis
    virtual void solve() = 0;


    /// @this will read from the input stream
    virtual std::istream& readFromInputStream(std::istream& input) = 0;

    
    /// overloaded operator to read from input stream
    friend std::istream& operator>>(std::istream& input, 
                                    Solution::SolutionBase& solution)
    {
      solution.readFromInputStream(input);
      return input;
    }
    
        
protected:
    

    /// adds a participating discipline to this solution
    void addParticipatingDiscipline(const unsigned int discipline_enum_ID);
    
    /// adds the solver ID for the solver class 
    void addSolverID(const unsigned int solver_class, 
		     const unsigned int solver_info_ID);


    /// this method will solve for all the load cases and sensitivity analysis
    void linearSolve(const unsigned int discipline_enum);
    
    
    /// this method will solve for all the load cases and sensitivity analysis
    void nonlinearSolve(const unsigned int discipline_enum);

    
    /// this method will solve for all load cases and sensitivity analysis
    void eigenSolve(const unsigned int discipilne_enum);


    
    /// this method will solve for all the load cases and sensitivity analysis
    void linearTransientSolve(const unsigned int discipline_enum);
    
    
    /// this method will solve for all the load cases and sensitivity analysis
    void nonlinearTransientSolve(const unsigned int discipline_enum);

    
    /// reference to the FESystemController class
    FESystem::FESystemController* fesystem_controller;
    
    /// enum ID of the solution
    const unsigned int solution_enum_ID;
    
    /// ID of this solution
    unsigned int ID;
    
    /// solver info IDs per solver class
    std::map<unsigned int, unsigned int> solver_info_ID_map;
    
    /// vector of load case IDs that will be used in the analysis
    std::vector<unsigned int> load_case_IDs;
    
    /// participating disciplines for this soltion
    std::set<unsigned int> discipline_enum_IDs;
  };
  
  /// @returns a pointer to the solution 
  std::auto_ptr<Solution::SolutionBase> createSolutionBase(const unsigned int enum_ID);
}



inline
unsigned int 
Solution::SolutionBase::getID() const
{
  return this->ID;
}


inline
const std::vector<unsigned int>& 
Solution::SolutionBase::getLoadCaseIDs() const
{
  return this->load_case_IDs;
}



#endif // __fesystem_solution_base_h__
