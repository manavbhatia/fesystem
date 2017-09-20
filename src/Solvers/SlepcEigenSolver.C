// $Id: SlepcEigenSolver.C,v 1.1.4.6 2008-08-25 04:40:03 manav Exp $


// FESystem includes
#include "Solvers/SlepcEigenSolver.h"
#include "Solvers/EigenSolverInfo.h"
#include "FESystem/FESystemController.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/PetscLinearSolver.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Utilities/TimeLogs.h"
#include "Utilities/ParallelUtility.h"

// libMesh includes
#include "numerics/petsc_matrix.h"
#include "numerics/petsc_vector.h"


Solver::SlepcEigenSolver::SlepcEigenSolver(const Solver::EigenSolverInfo& eigen_info,
                                 const Solver::LinearSolverInfo& linear_info):
Solver::EigenSolverBase(eigen_info, linear_info),
solution_completed(false),
n_converged_eigen_pairs(FESystemNumbers::InvalidID)
{
  PetscErrorCode ierr = 0;

  // create the solver
  ierr = EPSCreate(FESystem::COMM_WORLD, &(this->eps));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




Solver::SlepcEigenSolver::~SlepcEigenSolver()
{
  // destroy the solver context
  PetscErrorCode ierr = 0;
  ierr = EPSDestroy(this->eps);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




void
Solver::SlepcEigenSolver::clear()
{
  this->solution_completed = false;
  this->n_converged_eigen_pairs = FESystemNumbers::InvalidID;
  this->eigen_val_permutation.clear();
  
  Solver::FESystemSolverBase::clear();
}



unsigned int 
Solver::SlepcEigenSolver::getNConvergedEigenPairs()
{
  // make sure that the solution has been performed.
  Assert(this->solution_completed, ExcInvalidState());

  return this->n_converged_eigen_pairs;
}






void
Solver::SlepcEigenSolver::setNEigenSetsToCompute(const unsigned int n_eigs)
{
  PetscErrorCode ierr = 0;
  
  // set the number of eigen values to compute
  ierr = EPSSetDimensions(this->eps, 
                          n_eigs,
                          PETSC_DECIDE, PETSC_DECIDE); 
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}



void
Solver::SlepcEigenSolver::setEigenSpectrumEnd(const unsigned int eigen_spectrum_end)
{
  PetscErrorCode ierr = 0;
  
  switch(eigen_spectrum_end)
    {
    case LARGEST_MAGNITUDE_ENUM_ID:
      this->eps_spectrum_end = EPS_LARGEST_MAGNITUDE;
      break;

    case SMALLEST_MAGNITUDE_ENUM_ID:
      this->eps_spectrum_end = EPS_SMALLEST_MAGNITUDE;
      break;
      
    case LARGEST_REAL_ENUM_ID:
      this->eps_spectrum_end = EPS_LARGEST_REAL;
      break;

    case SMALLEST_REAL_ENUM_ID:
      this->eps_spectrum_end = EPS_SMALLEST_REAL;
      break;

    case LARGEST_IMAGINARY_ENUM_ID:
      this->eps_spectrum_end = EPS_LARGEST_IMAGINARY;
      break;

    case SMALLEST_IMAGINARY_ENUM_ID:
      this->eps_spectrum_end = EPS_SMALLEST_IMAGINARY;
      break;

    default: 
      AssertThrow(false, ExcInternalError());
    }
    
  ierr = EPSSetWhichEigenpairs(this->eps, this->eps_spectrum_end); 
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}





void
Solver::SlepcEigenSolver::setEigenShift(const unsigned int eigen_shift_kind,
                                   const double shift_value)
{
  PetscErrorCode ierr = 0;
  bool apply_shift = true;
  STType type;
    
  switch(eigen_shift_kind)
    {
    case NO_SHIFT_ENUM_ID:
      apply_shift = false;
      break;
      
    case ORIGIN_SHIFT_ENUM_ID:
      type = STSHIFT;
      break;
      
    case SPECTRUM_FOLD_ENUM_ID:
      type = STFOLD;
      break;
      
    case SHIFT_AND_INVERT_ENUM_ID:
      type = STSINV;
      break;
      
    case CAYLEY_SHIFT_ENUM_ID:
      type = STCAYLEY;
      break;
      
    default: 
      AssertThrow(false, ExcInternalError());
    }
  
  ST st;
  ierr = EPSGetST(this->eps, &st); 
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  if (apply_shift)
    {
    ierr = STSetType(st, type); 
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    
    ierr = STSetShift(st, shift_value); 
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    }
  else
    {
    // Slepc handles the no shift as STSHIFT with shift of 0.0
    ierr = STSetType(st, STSHIFT); 
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    
    ierr = STSetShift(st, 0.0); 
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    }
}



void
Solver::SlepcEigenSolver::setSolverType(const unsigned int enum_ID)
{
  PetscErrorCode ierr = 0;
  
  // and set the solver type
  switch (enum_ID)
    {
    case EPS_LAPACK_EIGEN_SOLVER_ENUM_ID:
      {
        ierr = EPSSetType(this->eps, EPSLAPACK);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case EPS_POWER_EIGEN_SOLVER_ENUM_ID:
      {
        ierr = EPSSetType(this->eps, EPSPOWER);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case EPS_SUBSPACE_EIGEN_SOLVER_ENUM_ID:
      {
        ierr = EPSSetType(this->eps, EPSSUBSPACE);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case EPS_ARNOLDI_EIGEN_SOLVER_ENUM_ID:
      {
        ierr = EPSSetType(this->eps, EPSARNOLDI);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case EPS_LANCZOS_EIGEN_SOLVER_ENUM_ID:
      {
        ierr = EPSSetType(this->eps, EPSLANCZOS);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_ID:
      {
        ierr = EPSSetType(this->eps, EPSKRYLOVSCHUR);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
    }
  
}

void
Solver::SlepcEigenSolver::setProblemKind(const unsigned int enum_ID)
{
  PetscErrorCode ierr = 0;

  // now set the  problem type 
  switch (enum_ID)
    {
    case HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        ierr = EPSSetProblemType(this->eps, EPS_HEP);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        ierr = EPSSetProblemType(this->eps, EPS_NHEP);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        ierr = EPSSetProblemType(this->eps, EPS_GHEP);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        ierr = EPSSetProblemType(this->eps, EPS_GNHEP);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
    }
}



void 
Solver::SlepcEigenSolver::setLinearSolverDetails()
{
  PetscErrorCode ierr = 0;

  // get ST context 
  ST st;
  ierr = EPSGetST(this->eps, &st); 
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->petsc_linear_solver.reset(new Solver::PetscLinearSolver(this->linear_solver_info));
  
  KSP ksp = this->petsc_linear_solver->getKSP();

  ierr = STSetKSP(st, ksp);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);


}







void 
Solver::SlepcEigenSolver::setSolverOptions()
{
  // the options for this solver are set based on how Slepc handles them in a sequence. 
  // the sequence in 2.3.2 is
  // type, problem type, class, orthogonalization, orthogonalization refinement, 
  // tolerances, dimensions, monitor, which eigenpairs, st_options
  // and the st options are handles in the following sequence
  // type, shift, mat mode, mat structure
  
  // now all of these options are handles in this solver, but the ones that are, will be set
  // in the same sequence.

  // type
  unsigned int solver_type = this->eigen_solver_info.getEigenSolverKindEnumID();
  this->setSolverType(solver_type);
  
  // problem type
  unsigned int problem_type_enum_ID = this->eigen_solver_info.getEigenProblemKindEnumID();
  this->setProblemKind(problem_type_enum_ID);
  
  // class not handled
  
  // orthogonalization not handled
  
  // orthogonalization refinement not handled
  
  // tolerances not handled
  
  // dimensions
  unsigned int n_eigen_sets = this->eigen_solver_info.getNEigenPairsToSolve();
  this->setNEigenSetsToCompute(n_eigen_sets);
  
  // monitor not handled
  
  // which eigenpairs
  unsigned int eigen_spectrum_end = this->eigen_solver_info.getEigenSpectrumToCompute();
  this->setEigenSpectrumEnd(eigen_spectrum_end);
  
  // ST type and shift
  unsigned int eigen_shift_kind = this->eigen_solver_info.getEigenSolverShiftKind();
  double shift_value = this->eigen_solver_info.getEigenSolverShiftValue();
  this->setEigenShift(eigen_shift_kind, shift_value);
  
  // set the details of the ksp and pc for the linear solver used in the solver
  this->setLinearSolverDetails();
}





void 
Solver::SlepcEigenSolver::solve(SparseMatrix<double>* A_mat, 
			   SparseMatrix<double>* B_mat)
{
  this->analysis_driver->getFESystemController().performance_logging->setEvent
  ("SlepcEigenSolver::solve()", "Solver");
  
  PetscErrorCode ierr = 0;
  
  Assert(A_mat != NULL, ExcInternalError());
  
  A_mat->close();
  PetscMatrix<double> *petsc_A_mat = dynamic_cast<PetscMatrix<double>*>(A_mat);

  unsigned int problem_type_enum_ID = this->eigen_solver_info.getEigenProblemKindEnumID();

  switch(problem_type_enum_ID)
    {
    case HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        // if the problem is not generalized eigen problem, then only the A_matrix 
        // needs to be specified
        Assert(B_mat == NULL, ExcInternalError());
        ierr = EPSSetOperators(this->eps, petsc_A_mat->mat(), PETSC_NULL);
        CHKERRABORT(FESystem::COMM_WORLD, ierr);  
      }
      break;

    case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        // if the problem is a generalized eigen problem, then both the A and B matrices 
        // need to be specified
        Assert(B_mat != NULL, ExcInternalError());
        B_mat->close();
        PetscMatrix<double> *petsc_B_mat = dynamic_cast<PetscMatrix<double>*>(B_mat);
        
        if (this->eigen_solver_info.ifSwapMatrices())
          ierr = EPSSetOperators(this->eps, petsc_B_mat->mat(), petsc_A_mat->mat());
        else
          ierr = EPSSetOperators(this->eps, petsc_A_mat->mat(), petsc_B_mat->mat());

        CHKERRABORT(FESystem::COMM_WORLD, ierr);
      }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
    }
  

  // before solution, set all the options of the solver
  this->setSolverOptions();
  
 // ierr = EPSSetType(eps, EPSLANCZOS); CHKERRABORT(FESystem::COMM_WORLD, ierr);  
//  
//  ierr = EPSSetProblemType(eps, EPS_GHEP); CHKERRABORT(FESystem::COMM_WORLD, ierr);  
//  
//  ierr = EPSSetDimensions(eps, 4, PETSC_DECIDE); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ST st;
//  ierr = EPSGetST(eps, &st); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = STSetType(st, STSINV); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = STSetShift(st, 40000.0); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = EPSGetST(eps, &st); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  KSP ksp;
//  PC pc;
//  
//  ierr = STGetKSP(st, &ksp); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = KSPGetPC(ksp, &pc); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = KSPSetType(ksp, KSPPREONLY); CHKERRABORT(FESystem::COMM_WORLD, ierr);
//  
//  ierr = PCSetType(pc, PCLU); CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // solve 
  ierr = EPSSolve(this->eps);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // get the reason for convergence of the solver
  EPSConvergedReason reason;
  ierr = EPSGetConvergedReason(this->eps, &reason);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  
  // get the number of converged eigenpairs
  ierr = EPSGetConverged(this->eps, &this->n_converged_eigen_pairs);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
	
  
  switch(reason)
    {
    case EPS_CONVERGED_TOL:
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD," EigenSolver converged up to tolerance.");
      }
      break;
      
    case EPS_DIVERGED_ITS:
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD,
                                            " EigenSolver did not converge in maximum iterations.");
      }      
      break;
      
    case EPS_DIVERGED_BREAKDOWN:
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD," EigenSolver broke down.");
      }      
      break;
      
    case EPS_DIVERGED_NONSYMMETRIC:
      {
        FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD," EigenSolver operator is non-symmetric.");
      }
      break;
      
    default:
      abort();
      break;
    }
  
  // once the solution is done, sort the eigen values according to the 
  // criteria for solution
  int *permute_vec = new int[this->n_converged_eigen_pairs];
  double *val_r_vec = new double[this->n_converged_eigen_pairs],
    *val_i_vec = new double[this->n_converged_eigen_pairs];
  
  double r_val, i_val;
  for (int i=0; i<this->n_converged_eigen_pairs; i++)
    {
    r_val = 0.0; i_val = 0.0;
    ierr = EPSGetValue(this->eps, i, &r_val, &i_val);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    val_r_vec[i] = r_val;
    val_i_vec[i] = i_val;
    }
  
  ierr = EPSSortEigenvalues(this->n_converged_eigen_pairs,
                            val_r_vec, val_i_vec,
                            this->eps_spectrum_end,
                            this->n_converged_eigen_pairs,
                            permute_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  // now store the permutation indices
  this->eigen_val_permutation.resize(this->n_converged_eigen_pairs);
  for (int i=0; i<this->n_converged_eigen_pairs; i++)
    this->eigen_val_permutation[i] = permute_vec[i];
  
  // now delete the vectors created
  delete[] permute_vec;
  delete[] val_r_vec;
  delete[] val_i_vec;
  
  this->solution_completed = true;

  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
    ("SlepcEigenSolver::solve()", "Solver");
}



void 
Solver::SlepcEigenSolver::getEigenValue(unsigned int eigen_index,
                                   double* real_value,
                                   double* img_value)
{
  // make sure that the requested index is less than the number of converged 
  // eigenpairs
  Assert(this->solution_completed, ExcInvalidState());
  Assert(this->n_converged_eigen_pairs > 0, ExcInternalError());
  Assert((int) eigen_index < this->n_converged_eigen_pairs, ExcInternalError());
  
  // now get the eigen pair
  PetscErrorCode ierr = 0;
  if (this->eigen_solver_info.ifSwapMatrices())
    {
    double rval=0.0, ival=0.0, mag=0.0;
    ierr = EPSGetValue(this->eps, this->eigen_val_permutation[eigen_index],
                       &rval, &ival);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    mag = sqrt(rval * rval + ival * ival);
    *real_value = rval / mag;
    *img_value = ival / mag;
    }
  else
    {
    ierr = EPSGetValue(this->eps, this->eigen_val_permutation[eigen_index],
                       real_value, img_value);
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    }
}




void
Solver::SlepcEigenSolver::getInvariantSubspace(std::vector<NumericVector<double>*>& vectors) 
{
  Assert(this->solution_completed, ExcInvalidState()); 
  Assert(vectors.size() == this->getNConvergedEigenPairs(), ExcInternalError());
  PetscErrorCode ierr = 0; 
  Vec pvec = PETSC_NULL;
  
  std::vector<Vec> vecs(this->getNConvergedEigenPairs());
  for (unsigned int i=0; i<this->getNConvergedEigenPairs(); i++)
    {
      pvec = dynamic_cast<PetscVector<double>*>(vectors[i])->vec();
      vecs[i] = pvec;
    }
  
  ierr = EPSGetInvariantSubspace(this->eps, &(vecs[0]));
  
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
}




void 
Solver::SlepcEigenSolver::getEigenPair(unsigned int pair_index,
                                  double* real_value,
                                  double* img_value, 
                                  NumericVector<double>& real_vec,
                                  NumericVector<double>& img_vec)
{
  // make sure that the requested index is less than the number of converged 
  // eigenpairs 
  Assert(this->solution_completed, ExcInvalidState()); 
  Assert(this->n_converged_eigen_pairs > 0, ExcInternalError()); 
  Assert((int) pair_index < this->n_converged_eigen_pairs, ExcInternalError()); 
  
  PetscVector<double>& real_vec_ref = dynamic_cast<PetscVector<double>& >(real_vec); 
  PetscVector<double>& img_vec_ref = dynamic_cast<PetscVector<double>& >(img_vec); 
  
  // now get the eigen pair
  PetscErrorCode ierr = 0; 
  if (this->eigen_solver_info.ifSwapMatrices())
    {
    double rval=0.0, ival=0.0, mag=0.0;
    ierr = EPSGetEigenpair(this->eps, this->eigen_val_permutation[pair_index],
                           &rval, &ival,
                           real_vec_ref.vec(), img_vec_ref.vec());
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    mag = rval * rval + ival * ival;
      *real_value = rval/mag;
      *img_value = -ival/mag;
    }
  else
    {
    ierr = EPSGetEigenpair(this->eps, this->eigen_val_permutation[pair_index],
                           real_value, img_value,
                           real_vec_ref.vec(), img_vec_ref.vec());
    CHKERRABORT(FESystem::COMM_WORLD, ierr);
    }
}


double
Solver::SlepcEigenSolver::getResidualForEigenPair(const unsigned int i)
{
  PetscErrorCode ierr = 0; 
  double return_val = 0.0;
  ierr = EPSComputeResidualNorm(this->eps, this->eigen_val_permutation[i], &return_val);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  return return_val;
}
