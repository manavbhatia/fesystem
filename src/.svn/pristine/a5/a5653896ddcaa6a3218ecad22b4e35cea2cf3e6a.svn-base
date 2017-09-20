// $Id: FESystemSolverBase.h,v 1.15.4.4 2007-05-11 05:16:54 manav Exp $

#ifndef __fesystem_solver_base_h__
#define __fesystem_solver_base_h__

// C++ includes
#include <string>
#include <map>
#include <memory>

// FESystem includes
#include "Utilities/NameEnumHandler.h"

namespace Driver
{
  class AnalysisDriver;
}



namespace Solver
{
  // Forward declerations
  class LinearSolverInfo;
  class NonlinearSolverInfo;
  class EigenSolverInfo;
  class SolverInfo;
  class LinearTransientSolverInfo;
  class NonlinearTransientSolverInfo;
  class LinearSolver;
  class NonlinearSolver;
  class EigenSolverBase;
  class LinearTransientSolverBase;
  class NonlinearTransientSolverBase;
  
  DeclareEnumClass(SolverClassEnum);
  
  
  /// this class will provide an interface to various solution algorithms
  /// available in the solver package. It can also be used to create an 
  /// algorithm that is unavailable in the solver package, and can still have 
  /// a uniform interface through the class methods
  /// each solver object will need a set of matrices and vectors for the solution
  /// process, which will vary depending upon the solver type. Before each solver
  /// step, pointers to these matrices and vectors will have to be set.
  class FESystemSolverBase
    {
public:	
      
      /// creator function will need a pointer to the AnalysisDriver that owns the 
      /// instantiation of this class.
      FESystemSolverBase(const Solver::SolverInfo& info);
      
      /// destructor fuction
      virtual ~FESystemSolverBase() = 0;

      // clears the data structures
      virtual void clear();
      
      /// returns the solver type
      inline unsigned int getSolverClassEnumID() const;
      
      inline std::string getSolverClassEnumName() const;
      
      // attaches the analysis driver for use in solution process
      void attachAnalysisDriver(Driver::AnalysisDriver* driver);
      
      // @returns a reference to the analysis driver
      Driver::AnalysisDriver& getAnalysisDriver();

protected:	

//        /// method to clear the data structures after a use
//        void clearFESystemSolverBase();
        
        /// analysis driver pointer
        Driver::AnalysisDriver* analysis_driver;
      
      unsigned int solver_class_enum_ID;
    };
  
  
  
  std::auto_ptr<Solver::LinearSolver>
    createLinearSolver(const Solver::LinearSolverInfo& info);
  
  
  
  std::auto_ptr<Solver::NonlinearSolver> 
    createNonlinearSolver(const Solver::NonlinearSolverInfo& info,
                          const Solver::LinearSolverInfo& linear_info);
  
  
  std::auto_ptr<Solver::EigenSolverBase> 
    createEigenSolver(const Solver::EigenSolverInfo& eigen_info,
                      const Solver::LinearSolverInfo& linear_info);


  std::auto_ptr<Solver::LinearTransientSolverBase> 
    createLinearTransientSolver(const Solver::LinearTransientSolverInfo& transient_info,
                                const Solver::LinearSolverInfo& linear_info);

  std::auto_ptr<Solver::NonlinearTransientSolverBase> 
    createNonlinearTransientSolver(const Solver::NonlinearTransientSolverInfo& transient_info,
                                   const Solver::NonlinearSolverInfo& nonlinear_info,
                                   const Solver::LinearSolverInfo& linear_info);
  
}



inline
unsigned int 
Solver::FESystemSolverBase::getSolverClassEnumID() const
{
  return this->solver_class_enum_ID;
}





inline 
std::string
Solver::FESystemSolverBase::getSolverClassEnumName() const
{
  return Solver::SolverClassEnum::enumName(this->solver_class_enum_ID);
}





#endif // __fesystem_solver_base_h__
