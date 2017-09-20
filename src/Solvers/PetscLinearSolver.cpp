// $Id:$
/*
 *  PetscLinearSolver.cpp
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 12/4/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

// C++ includes


// FESystem include
#include "Solvers/PetscLinearSolver.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Utilities/TimeLogs.h"
#include "Utilities/ParallelUtility.h"



Solver::PetscLinearSolver::PetscLinearSolver(const LinearSolverInfo& info):
Solver::LinearSolver(info, PETSC_LINEAR_SOLVER::num()),
ksp(PETSC_NULL),
ksptype(PETSC_NULL),
pc(PETSC_NULL),
pctype(PETSC_NULL),
petsc_draw(PETSC_NULL)
{
    // set the pc type
  switch(this->pc_type_enum_ID)
  {
    case PC_LU_ENUM_ID:
      this->pctype = PCLU;
      break;
      
    case PC_JACOBI_ENUM_ID:
      this->pctype = PCJACOBI;
      break;
      
    case PC_BJACOBI_ENUM_ID:
      this->pctype = PCBJACOBI;
      break;
      
    case PC_SOR_ENUM_ID:
      this->pctype = PCSOR;
      break;
      
    case PC_EISENSTAT_ENUM_ID:
      this->pctype = PCEISENSTAT;
      break;
      
    case PC_ICC_ENUM_ID:
      this->pctype = PCICC;
      break;
      
    case PC_ILU_ENUM_ID:
      this->pctype = PCILU;
      break;
      
    case PC_ASM_ENUM_ID:
      this->pctype = PCASM;
      break;
      
    case PC_KSP_ENUM_ID:
      this->pctype = PCKSP;
      break;
      
    case PC_CHOLESKY_ENUM_ID:
      this->pctype = PCCHOLESKY;
      break;
      
    case PC_NONE_ENUM_ID:
      this->pctype = PCNONE;
      break;
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (Solver::PCTypeEnum::enumName(this->pc_type_enum_ID)));
  }
  
// set the ksp type
  switch(this->ksp_type_enum_ID)
  {
    case KSP_RICHARDSON_ENUM_ID:
      this->ksptype = KSPRICHARDSON;
      break;
      
    case KSP_CHEBYCHEV_ENUM_ID:
      this->ksptype = KSPCHEBYCHEV;
      break;
      
    case KSP_CG_ENUM_ID:
      this->ksptype = KSPCG;
      break;
      
    case KSP_GMRES_ENUM_ID:
      this->ksptype = KSPGMRES;
      break;
      
    case KSP_TCQMR_ENUM_ID:
      this->ksptype = KSPTCQMR;
      break;
      
    case KSP_BCGS_ENUM_ID:
      this->ksptype = KSPBCGS;
      break;
    case KSP_CGS_ENUM_ID:
      this->ksptype = KSPCGS;
      break;
      
    case KSP_TFQMR_ENUM_ID:
      this->ksptype = KSPTFQMR;
      break;
      
    case KSP_CR_ENUM_ID:
      this->ksptype = KSPCR;
      break;
      
    case KSP_LSQR_ENUM_ID:
      this->ksptype = KSPLSQR;
      break;
      
    case KSP_BICG_ENUM_ID:
      this->ksptype = KSPBICG;
      break;

    case KSP_PRE_ONLY_ENUM_ID:
      this->ksptype = KSPPREONLY;
      break;
      
    default:
      Assert(false,
             FESystemExceptions::ExcEnumCaseNotImplemented
             (Solver::KSPTypeEnum::enumName(this->ksp_type_enum_ID)));
  }
  
  this->createSolverContext();
  
}





Solver::PetscLinearSolver::~PetscLinearSolver()
{
    PetscErrorCode ierr=0;
  
  // the petsc draw context is to be deleted only on the 0 processor
  if (FESystem::local_processor == 0)
    {
      if (this->petsc_draw != PETSC_NULL)
        {
          ierr = KSPMonitorLGDestroy(this->petsc_draw);
          CHKERRABORT(FESystem::COMM_WORLD, ierr);
          this->petsc_draw = PETSC_NULL;
        }
    }

  // delete the PC and KSP contexts if they have been created
  if (this->pc != PETSC_NULL)
  {
    ierr = PCDestroy(this->pc);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    this->pc = PETSC_NULL;
  }
  
  if (this->ksp != PETSC_NULL)
  {
    ierr = KSPDestroy(this->ksp);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    this->ksp = PETSC_NULL;
  }
  
}



void
Solver::PetscLinearSolver::createSolverContext()
{
  PetscErrorCode ierr = 0;

  // the petsc draw context is to be deleted only on the 0 processor
  if (FESystem::local_processor == 0)
    {
      if (this->petsc_draw != PETSC_NULL)
        {
          ierr = KSPMonitorLGDestroy(this->petsc_draw);
          CHKERRABORT(FESystem::COMM_WORLD, ierr);
          this->petsc_draw = PETSC_NULL;
        }
    }

  if (this->pc != PETSC_NULL)
    {
      ierr = PCDestroy(this->pc);
      CHKERRABORT(FESystem::COMM_WORLD, ierr);
      this->pc = PETSC_NULL;
    }

  if (this->ksp != PETSC_NULL)
    {
      ierr = KSPDestroy(this->ksp); 
      CHKERRABORT(FESystem::COMM_WORLD, ierr); 
      this->ksp = PETSC_NULL;
    }

  {
    ierr = KSPCreate(FESystem::COMM_WORLD, &(this->ksp));
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    
    ierr = KSPSetType(this->ksp, this->ksptype);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    
    if (this->pc_type_enum_ID != Solver::PC_NONE::num())
      {
        ierr = PCCreate(FESystem::COMM_WORLD, &(this->pc));
        CHKERRABORT(FESystem::COMM_WORLD, ierr);
        
        ierr = PCSetType(this->pc, this->pctype);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);
        
        ierr = KSPSetPC(this->ksp, this->pc);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);
      }
  }
  
//  ierr = KSPMonitorLGTrueResidualNormCreate
//  (FESystem::COMM_WORLD, PETSC_NULL, "residual norm", 0, 0, 500, 500, &(this->petsc_draw));
//  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}



KSP 
Solver::PetscLinearSolver::getKSP()
{
  return this->ksp;
}




PC 
Solver::PetscLinearSolver::getPC()
{
  return this->pc;
}




void 
Solver::PetscLinearSolver::clear()
{
  
  Solver::LinearSolver::clear();
}






void 
Solver::PetscLinearSolver::solve(NumericVector<double>& rhs,
                                 NumericVector<double>& sol)
{
  
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("DirectLinearSolver::solve()", "Solver");
  
  PetscErrorCode ierr = 0;
  
  // next, solve
  PetscVector<double>& rhs_ref = dynamic_cast<PetscVector<double>& >(rhs);
  PetscVector<double>& sol_ref = dynamic_cast<PetscVector<double>& >(sol);
  
  MatStructure petsc_pc_option;

  // if the same matrix option between solves is false, then new matrices must 
  // be specified, otherwise, there should be no change between the solves
  if (this->same_matrices_between_solves)
  {
    if (this->new_PC_matrix || this->new_system_matrix)
    {
      // if same matrix is to be used, and new ones have been specified, then
      // the number of time these have been used should be zero
      Assert (this->n_solves_after_new_matrices == 0, ExcInternalError());
    }
    else 
    {
      // if new matrices have not been given, then these must have been used
      Assert (this->n_solves_after_new_matrices > 0, ExcInternalError());
    }
    petsc_pc_option = SAME_PRECONDITIONER;
  }
  else
  {
    // if new matrices have to be given at each iteration, then 
    // the number of solves after new matrices should be 0
    Assert (this->n_solves_after_new_matrices == 0, ExcInternalError());
    petsc_pc_option = DIFFERENT_NONZERO_PATTERN;
  }
  
  // attach the matrices to the solver before solving if new ones have been specified
  if (this->new_system_matrix || this->new_PC_matrix)
  {
    this->createSolverContext();
    Assert(this->system_matrix != NULL, ExcInternalError());
    this->system_matrix->close();
    Mat sys_mat = (dynamic_cast<PetscMatrix<double>* >(this->system_matrix))->mat();
    if (this->pc_type_enum_ID != Solver::PC_NONE::num())
    {
      Assert(this->preconditioner_matrix != NULL, ExcInternalError());
      this->preconditioner_matrix->close();
      Mat pc_mat = (dynamic_cast<PetscMatrix<double>* >(this->preconditioner_matrix))->mat();
      ierr = KSPSetOperators(this->ksp, sys_mat, pc_mat, petsc_pc_option);
    }
    else
      ierr = KSPSetOperators(this->ksp, sys_mat, sys_mat, petsc_pc_option);

    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    ierr = KSPSetUp(this->ksp);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
  }

  rhs_ref.close();
  sol_ref.close();
    
  ierr = KSPSolve(this->ksp, rhs_ref.vec(), sol_ref.vec());
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->n_solves_after_new_matrices += 1;
  // the matrices now are no longer new
  this->new_system_matrix = false;
  this->new_PC_matrix = false;
  
  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
  ("DirectLinearSolver::solve()", "Solver");
  
}





