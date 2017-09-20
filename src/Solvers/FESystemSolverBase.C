// $Id: FESystemSolverBase.C,v 1.20.4.6 2008-02-25 04:31:50 manav Exp $

//C++ includes
#include <memory>


// FESystem includes
#include "AnalysisDriver/AnalysisDriver.h"
#include "Solvers/FESystemSolverBase.h"
#include "FESystem/FESystemExceptions.h"
#include "FESystem/FESystemNumbers.h"
#include "Solvers/SolverInfo.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/EigenSolverInfo.h"
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/PetscLinearSolver.h"
#include "Solvers/NewtonNonlinearSolver.h"
#include "Solvers/PetscNonlinearSolver.h"
#include "Solvers/SlepcEigenSolver.h"
#include "Solvers/ArpackEigenSolver.h"
#include "Solvers/NewmarkTransientSolverInfo.h"
#include "Solvers/NewmarkTransientSolver.h"
#include "Solvers/EulerTransientSolverInfo.h"
#include "Solvers/EulerTransientSolver.h"

Solver::FESystemSolverBase::FESystemSolverBase(const Solver::SolverInfo& info):
analysis_driver(NULL),
solver_class_enum_ID(info.getSolverClassEnumID())
{
}





Solver::FESystemSolverBase::~FESystemSolverBase()
{
}



void
Solver::FESystemSolverBase::clear()
{
  this->analysis_driver = NULL;
}


void
Solver::FESystemSolverBase::attachAnalysisDriver(Driver::AnalysisDriver* driver)
{
  Assert(driver != NULL, ExcInternalError());
  
  this->analysis_driver = driver;
}



Driver::AnalysisDriver&
Solver::FESystemSolverBase::getAnalysisDriver()
{
  Assert(this->analysis_driver != NULL, ExcInternalError());
  
  return *(this->analysis_driver);
}



std::auto_ptr<Solver::LinearSolver>
Solver::createLinearSolver(const Solver::LinearSolverInfo& info)
{
  std::auto_ptr<Solver::LinearSolver> solver_ptr(NULL);
  
  unsigned int solver_enum_ID = info.getLinearSolverEnumID();
  switch(solver_enum_ID)
    {
    case PETSC_LINEAR_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::PetscLinearSolver(info));
      }
      break;
      
    default:
      Assert(false, ExcInternalError());
    };
  
  return solver_ptr;
}



std::auto_ptr<Solver::NonlinearSolver> 
Solver::createNonlinearSolver(const Solver::NonlinearSolverInfo& info,
                              const Solver::LinearSolverInfo& linear_info)
{
  std::auto_ptr<Solver::NonlinearSolver> solver_ptr(NULL);
  
  unsigned int nonlinear_solver_type = info.getNonlinearSolverKindEnumID();
  switch(nonlinear_solver_type)
    {
    case FESYSTEM_NEWTON_NONLINEAR_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::NewtonNonlinearSolver(info, linear_info));
      }
      break;
      
    case NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH_ENUM_ID:
    case NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH_ENUM_ID:
    case NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH_ENUM_ID:
    case NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH_ENUM_ID:
      {
          //solver_ptr.reset(new Solver::PetscNonlinearSolver(info, linear_info));
      }
      break;

    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (Solver::NonlinearSolverKindEnum::enumName(nonlinear_solver_type)));
    };
  
  return solver_ptr;
}


std::auto_ptr<Solver::EigenSolverBase> 
Solver::createEigenSolver(const Solver::EigenSolverInfo& eigen_info, 
                          const Solver::LinearSolverInfo& linear_info)
{
  std::auto_ptr<Solver::EigenSolverBase> solver_ptr(NULL);
  
  unsigned int eigen_solver_type = eigen_info.getEigenSolverKindEnumID();
  switch(eigen_solver_type)
    {
    case EPS_LAPACK_EIGEN_SOLVER_ENUM_ID:
    case EPS_POWER_EIGEN_SOLVER_ENUM_ID:
    case EPS_SUBSPACE_EIGEN_SOLVER_ENUM_ID:
    case EPS_ARNOLDI_EIGEN_SOLVER_ENUM_ID:
    case EPS_LANCZOS_EIGEN_SOLVER_ENUM_ID:
    case EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::SlepcEigenSolver(eigen_info, linear_info));
      }
      break;
      
    case ARPACK_EIGEN_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::ArpackEigenSolver(eigen_info, linear_info));
      }
      break;

    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (Solver::EigenSolverKindEnum::enumName(eigen_solver_type)));
    };
  
  return solver_ptr;
}


std::auto_ptr<Solver::LinearTransientSolverBase> 
Solver::createLinearTransientSolver(const Solver::LinearTransientSolverInfo& transient_info,
                            const Solver::LinearSolverInfo& linear_info)
{
  std::auto_ptr<Solver::LinearTransientSolverBase> solver_ptr(NULL);

  switch (transient_info.getTransientSolverKindEnumID())
    {
    case NEWMARK_LINEAR_TRANSIENT_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::NewmarkLinearTransientSolver(transient_info, linear_info));
      }
      break;
      
      case EULER_LINEAR_TRANSIENT_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::EulerLinearTransientSolver(transient_info, linear_info));
      }
        break;

      default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (Solver::TransientSolverKindEnum::enumName
              (transient_info.getTransientSolverKindEnumID())));
    };
  
  return solver_ptr;
}




std::auto_ptr<Solver::NonlinearTransientSolverBase> 
Solver::createNonlinearTransientSolver(const Solver::NonlinearTransientSolverInfo& transient_info,
                                       const Solver::NonlinearSolverInfo& nonlinear_info,
                                       const Solver::LinearSolverInfo& linear_info)
{
  std::auto_ptr<Solver::NonlinearTransientSolverBase> solver_ptr(NULL);
  
  switch (transient_info.getTransientSolverKindEnumID())
    {
    case NEWMARK_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::NewmarkNonlinearTransientSolver
                         (transient_info, nonlinear_info, linear_info));
      }
      break;
      
      case EULER_NONLINEAR_TRANSIENT_SOLVER_ENUM_ID:
      {
        solver_ptr.reset(new Solver::EulerNonlinearTransientSolver
                         (transient_info, nonlinear_info, linear_info));
      }
        break;

      default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (Solver::TransientSolverKindEnum::enumName
              (transient_info.getTransientSolverKindEnumID())));
    };
  
  return solver_ptr;
}

