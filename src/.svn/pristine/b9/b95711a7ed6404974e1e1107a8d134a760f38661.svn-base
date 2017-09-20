// $Id: SolverInfo.C,v 1.4.6.2 2008-02-25 04:31:50 manav Exp $

// FESystem includes
#include "Solvers/SolverInfo.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/EigenSolverInfo.h"
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/NewmarkTransientSolverInfo.h"
#include "Solvers/EulerTransientSolverInfo.h"



Solver::SolverInfo::SolverInfo(const unsigned int info_enum_ID,
                               const unsigned int solver_class):
ID(FESystemNumbers::InvalidID),
solver_info_enum_ID(info_enum_ID),
solver_class_enum_ID(solver_class),
initialized(false)
{
  
}


Solver::SolverInfo::~SolverInfo()
{
  
}

std::istream&
Solver::operator>>(std::istream& input, Solver::SolverInfo& info)
{
  info.readFromInputStream(input);
  return input;
}



std::auto_ptr<Solver::SolverInfo>
Solver::createSolverInfo(const unsigned int enum_ID)
{
  std::auto_ptr<Solver::SolverInfo> info(NULL);
  
  switch (enum_ID)
    {
    case LINEAR_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::LinearSolverInfo());
      break;
      
    case NONLINEAR_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::NonlinearSolverInfo());
      break;
      
    case EIGEN_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::EigenSolverInfo());
      break;
      
    case LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::LinearTransientSolverInfo());
      break;

    case NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::NonlinearTransientSolverInfo());
      break;

    case NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::NewmarkLinearTransientSolverInfo());
      break;
      
    case NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID:
      info.reset(new Solver::NewmarkNonlinearTransientSolverInfo());
      break;

      case EULER_LINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID:
        info.reset(new Solver::EulerLinearTransientSolverInfo());
        break;
        
      case EULER_NONLINEAR_TRANSIENT_SOLVER_INFO_ENUM_ID:
        info.reset(new Solver::EulerNonlinearTransientSolverInfo());
        break;
        
    default:
      Assert(false, ExcInternalError());
    }
  
  return info;
}

