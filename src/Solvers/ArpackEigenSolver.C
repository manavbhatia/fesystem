// $Id: ArpackEigenSolver.C,v 1.1.4.6 2007-06-13 15:01:46 manav Exp $

// FESystem includes
#include "Solvers/ArpackEigenSolver.h"
#include "Solvers/EigenSolverInfo.h"
#include "FESystem/FESystemController.h"
#include "Solvers/LinearSolverInfo.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Utilities/TimeLogs.h"

// libMesh includes
#include "numerics/petsc_matrix.h"
#include "numerics/petsc_vector.h"


Solver::ArpackEigenSolver::ArpackEigenSolver(const Solver::EigenSolverInfo& eigen_info,
                                             const Solver::LinearSolverInfo& linear_info):
Solver::EigenSolverBase(eigen_info, linear_info),
solution_completed(false),
n_converged_eigen_pairs(FESystemNumbers::InvalidID),
ido(0),
n(0),
nev(0),
tol(0.0),
ncv(0),
ldv(0),
lworkl(0),
info(0),
rvec(0),
ldz(0),
sigmar(0.0),
sigmai(0.0)
{
  AssertThrow(FESystem::total_processors == 1, ExcInternalError());
}




Solver::ArpackEigenSolver::~ArpackEigenSolver()
{
  this->clear();
}





void
Solver::ArpackEigenSolver::intermediateClear()
{
  this->solution_completed = false;
  
  /// number of converged eigenpairs
  this->n_converged_eigen_pairs = FESystemNumbers::InvalidID;
  
  this->ido = FESystemNumbers::InvalidID;
  this->n = FESystemNumbers::InvalidID;
  this->nev = FESystemNumbers::InvalidID;
  this->tol = 0.0;
  this->ncv = FESystemNumbers::InvalidID;
  this->ldv = FESystemNumbers::InvalidID;
  this->lworkl = FESystemNumbers::InvalidID;
  this->info = FESystemNumbers::InvalidID;
  this->rvec = FESystemNumbers::InvalidID;
  this->ldz = FESystemNumbers::InvalidID;
  this->sigmar = 0.0;
  this->sigmai = 0.0;
  this->bmat.clear();
  this->which.clear();
  this->select.clear();
  this->HowMny.clear();
  this->Z.clear();
  this->dr.clear();
  this->di.clear();
  this->V.clear();
  this->resid.clear();
  this->iparam.clear();
  this->ipntr.clear();
  this->workd.clear();
  this->workl.clear();
  this->workev.clear();
  this->operator_matrix.reset();
  //this->linear_solver.reset();
}


void
Solver::ArpackEigenSolver::clear()
{
  this->intermediateClear();
  if (this->linear_solver.get() != NULL)
    this->linear_solver->clear();
  
  Solver::EigenSolverBase::clear();
}



void 
Solver::ArpackEigenSolver::prepare(const unsigned int dimensions)
{
  // make sure that the solution has been performed.
  Assert(!this->solution_completed, ExcInvalidState());
  this->ido = 0;
  this->n = dimensions;
  this->nev = this->eigen_solver_info.getNEigenPairsToSolve();
  this->tol = this->eigen_solver_info.getTolerance();
  // a higher than recommended value of ncv is kept. 
  this->ncv = 3 * this->nev+1;
  if (this->ncv > this->n) this->ncv = this->n;
  this->ldv = this->n;
  this->info = 0;
  if (this->eigen_solver_info.calculateEigenVectors())
    {
    this->rvec = 1;
    this->HowMny = "A";
      // the size for hermitian is nev*n + 1, and for non-hermitian it 
      // is (nev+1)*n+1. Since both have the same number of rows, it is 
      // safe to use just the bigger of the two number of columns. 
      this->Z.resize((this->nev+1) * this->n + 1);
      std::fill(this->Z.begin(), this->Z.end(), 0.0);
    }
  else
    this->rvec = 0;
  // select should contain the selected eigen vectors to be computed, but we are computing 
  // all by default
  this->select.resize(this->ncv + 1);
  std::fill(this->select.begin(), this->select.end(), 0);
  this->ldz = this->n;
  

  switch (this->eigen_solver_info.getEigenProblemKindEnumID())
    {
    case HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      this->bmat = "I";
      break;
      
    case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      this->bmat = "G";
      break;

    default:
      AssertThrow(false, ExcInternalError());
      break;
    }


  switch (this->eigen_solver_info.getEigenProblemKindEnumID())
  {
    case HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    {
      this->lworkl = this->ncv * (this->ncv + 8);
      this->ipntr.resize(12);
      std::fill(this->ipntr.begin(), this->ipntr.end(), 0);


      switch(this->eigen_solver_info.getEigenSpectrumToCompute())
      {
        case LARGEST_MAGNITUDE_ENUM_ID:
          this->which = "LM";
          break;
          
        case SMALLEST_MAGNITUDE_ENUM_ID:
          this->which = "SM";
          break;
          
        case LARGEST_REAL_ENUM_ID:
          this->which = "LA";
          break;
          
        case SMALLEST_REAL_ENUM_ID:
          this->which = "SA";
          break;
          
        case LARGEST_IMAGINARY_ENUM_ID:
        case SMALLEST_IMAGINARY_ENUM_ID:
        default: 
          AssertThrow(false, ExcInternalError());
      }
    }
      break;
  

    case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    {
      this->lworkl = 2*( this->ncv * (3 * this->ncv + 6));
      this->ipntr.resize(15);
      std::fill(this->ipntr.begin(), this->ipntr.end(), 0);
      this->workev.resize(3* this->ncv+1);
      std::fill(this->workev.begin(), this->workev.end(), 0.0);
      
      switch(this->eigen_solver_info.getEigenSpectrumToCompute())
      {
        case LARGEST_MAGNITUDE_ENUM_ID:
          this->which = "LM";
          break;
          
        case SMALLEST_MAGNITUDE_ENUM_ID:
          this->which = "SM";
          break;
          
        case LARGEST_REAL_ENUM_ID:
          this->which = "LR";
          break;
          
        case SMALLEST_REAL_ENUM_ID:
          this->which = "SR";
          break;
          
        case LARGEST_IMAGINARY_ENUM_ID:
          this->which = "LI";
          break;
          
        case SMALLEST_IMAGINARY_ENUM_ID:
          this->which = "SI";
          break;
          
        default: 
          AssertThrow(false, ExcInternalError());
      }
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
  }


  this->dr.resize(this->nev+1+1);
  std::fill(this->dr.begin(), this->dr.end(), 0.0);
  this->di.resize(this->nev+1+1);
  std::fill(this->di.begin(), this->di.end(), 0.0);
  this->V.resize(this->ncv * this->n + 1);
  std::fill(this->V.begin(), this->V.end(), 0.0);
  this->resid.resize(this->n+1);
  std::fill(this->resid.begin(), this->resid.end(), 0.0);
  this->iparam.resize(12);
  std::fill(this->iparam.begin(), this->iparam.end(), 0);
  this->workd.resize(3 * this->n + 1);
  std::fill(this->workd.begin(), this->workd.end(), 0.0);
  this->workl.resize(this->lworkl + 1);   
  std::fill(this->workl.begin(), this->workl.end(), 0.0);

  //  tell the solver that exact shifts are used
  this->iparam[1] = 1;
  
  // set the max allowable iterations
  this->iparam[3] = this->eigen_solver_info.getMaxAllowableIterations();
  
  // now set the problem shift type
  switch (this->eigen_solver_info.getEigenSolverShiftKind())
    {
    case NO_SHIFT_ENUM_ID:
      if (this->bmat == "I")
        this->iparam[7] = 1;
      else 
        this->iparam[7] = 2;
      break;
      
    case SHIFT_AND_INVERT_ENUM_ID:
      this->iparam[7] = 3;
      break;
      
    case CAYLEY_SHIFT_ENUM_ID:
      {
        switch (this->eigen_solver_info.getEigenProblemKindEnumID())
        {
        case HERMITIAN_EIGENPROBLEM_ENUM_ID:
        case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
          this->iparam[7] = 5;
          break;
          
        case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
        case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
        default:
          AssertThrow(false, ExcInternalError());
          break;
        }
      }
      break;
      
    case SPECTRUM_FOLD_ENUM_ID:
    case ORIGIN_SHIFT_ENUM_ID:
    default: 
      AssertThrow(false, ExcInternalError());
    }
  
  // only sigmar is used, since the imaginary part of the shift is not allowed
  if (this->iparam[7] > 2)
    this->sigmar = this->eigen_solver_info.getEigenSolverShiftValue();
}



void 
Solver::ArpackEigenSolver::initLinearSolver(SparseMatrix<double>* A_mat,
                                            SparseMatrix<double>* B_mat)
{
  Assert(this->analysis_driver != NULL, ExcInternalError());
  
  if (this->linear_solver.get() == NULL)
    {
    this->linear_solver.reset
    (dynamic_cast<Solver::LinearSolver*>
     (Solver::createLinearSolver(this->linear_solver_info).release()));
    this->linear_solver->attachAnalysisDriver(this->analysis_driver);
    }
  
  Assert (A_mat != NULL, ExcEmptyObject());
  
  Assert(this->operator_matrix.get() == NULL, ExcInternalError());
  this->operator_matrix.reset(SparseMatrix<double>::build().release());
  
  switch (this->eigen_solver_info.getEigenSolverShiftKind())
    {
    case NO_SHIFT_ENUM_ID:
      if (this->bmat == "G")
        {
        Assert (B_mat != NULL, ExcEmptyObject());
        this->operator_matrix->duplicate_matrix(*B_mat, true);
        this->operator_matrix->close();
        this->linear_solver->setSystemMatrix(*this->operator_matrix);
          this->linear_solver->setPreconditionerMatrix(*this->operator_matrix);
          this->linear_solver->useSameMatricesBetweenSolves();
        }
      break;
      
    case SHIFT_AND_INVERT_ENUM_ID:
      switch (this->eigen_solver_info.getEigenProblemKindEnumID())
        {
        case HERMITIAN_EIGENPROBLEM_ENUM_ID:
        case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
         {
            this->operator_matrix->duplicate_matrix(*A_mat, true);
            this->operator_matrix->shift(- this->sigmar);
            this->operator_matrix->close();
            this->linear_solver->setSystemMatrix(*this->operator_matrix);
            this->linear_solver->setPreconditionerMatrix(*this->operator_matrix);
            this->linear_solver->useSameMatricesBetweenSolves();
          }
          break;
          
        case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
        case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
          {
            Assert(B_mat != NULL, ExcEmptyObject());
            this->operator_matrix->duplicate_matrix(*A_mat, true);
            this->operator_matrix->add(-this->sigmar, *B_mat);
            this->operator_matrix->close();
            this->linear_solver->setSystemMatrix(*this->operator_matrix);
            this->linear_solver->setPreconditionerMatrix(*this->operator_matrix);
            this->linear_solver->useSameMatricesBetweenSolves();
          }
          break;
          
        default:
          AssertThrow(false, ExcInternalError());
          break;
        }
      break;
      
    case CAYLEY_SHIFT_ENUM_ID:
    case SPECTRUM_FOLD_ENUM_ID:
    case ORIGIN_SHIFT_ENUM_ID:
    default: 
      AssertThrow(false, ExcInternalError());
    }
}




unsigned int 
Solver::ArpackEigenSolver::getNConvergedEigenPairs()
{
  // make sure that the solution has been performed.
  Assert(this->solution_completed, ExcInvalidState());
  
  return this->iparam[5];
}



void 
Solver::ArpackEigenSolver::solve(SparseMatrix<double>* A_mat, 
                                 SparseMatrix<double>* B_mat)
{
  
  if (this->solution_completed)
    this->intermediateClear();
  Assert(A_mat != NULL, ExcInternalError());
  Assert(A_mat->m() == A_mat->n(), ExcInternalError());
  
  this->analysis_driver->getFESystemController().performance_logging->setEvent
    ("ArpackEigenSolver::solve()", "Solver");

  A_mat->close();
  SparseMatrix<double>* solver_A_mat = NULL, *solver_B_mat = NULL;
  unsigned int problem_type_enum_ID = this->eigen_solver_info.getEigenProblemKindEnumID();

  switch(problem_type_enum_ID)
    {
    case HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        Assert(B_mat == NULL, ExcInternalError());
        AssertThrow(!this->eigen_solver_info.ifSwapMatrices(), 
                    ExcInternalError());
        solver_A_mat = A_mat;
      }
      break;

    case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
      {
        Assert(B_mat != NULL, ExcInternalError());
        Assert (B_mat->m() == B_mat->n(), ExcInternalError());
        Assert (A_mat->m() == B_mat->m(), ExcInternalError());
        B_mat->close();
        
        if (this->eigen_solver_info.ifSwapMatrices())
          {
          solver_A_mat = B_mat;
          solver_B_mat = A_mat;
          }
        else
          {
          solver_A_mat = A_mat;
          solver_B_mat = B_mat;
          }
        
  
      }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
    }
  

  
  
  // before solution, set all the options of the solver
  this->prepare(solver_A_mat->m());
  this->initLinearSolver(solver_A_mat, solver_B_mat);

  // create the vectors for matrix operations
  std::auto_ptr<NumericVector<double> > 
    vec1(NumericVector<double>::build().release()), vec2(NumericVector<double>::build().release());
  
  vec1->init(solver_A_mat->m());
  vec2->init(solver_A_mat->m());
  
  // solve 
  while (this->ido != 99)
    {
    switch (problem_type_enum_ID)
      {
      // this is for real symmetric
      case HERMITIAN_EIGENPROBLEM_ENUM_ID:
      case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
        FC_FUNC(dsaupd, DSAUPD)(&ido, const_cast<char*>(bmat.c_str()), 
                        &n, const_cast<char*>(which.c_str()), 
                        &nev, &tol,
                        &resid[1], &ncv,
                        &V[1], &ldv,
                        &iparam[1], &ipntr[1],
                        &workd[1], &workl[1],
                        &lworkl, &info);
        break;
        
          // this is for real un-symmetric
        case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
        case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
          FC_FUNC(dnaupd, DNAUPD)(&ido, const_cast<char*>(bmat.c_str()), 
                                  &n, const_cast<char*>(which.c_str()), 
                                  &nev, &tol,
                                  &resid[1], &ncv,
                                  &V[1], &ldv,
                                  &iparam[1], &ipntr[1],
                                  &workd[1], &workl[1],
                                  &lworkl, &info);
          break;

        default:
        AssertThrow(false, ExcInternalError());
        break;
      }


    // check the error, if any
    this->checkError(this->info);

    // copy the vector into vec1
    for (unsigned int i=0; i<solver_A_mat->m(); i++)
      vec1->set(i, workd[ipntr[1]+i]);
    
    vec2->zero();
    
    // this is for real symmetric problems only
    switch (this->ido)
      {
      case -1:
      case 1:
        {
          // compute  Y = OP * X  for case -1
          // compute  Y = OP * X without performing the B * X operation in 3, 4 and 5 for case 1
          switch (this->iparam[7])
            {
            case 1:
              {
                // OP = A, B = I
                solver_A_mat->multiply_vector(*vec1, *vec2);
              }
            break;

            case 2:
              {
                // OP = M^(-1) A , B = M
                solver_A_mat->multiply_vector(*vec1, *vec2);
                vec1->zero();
                this->linear_solver->solve(*vec2, *vec1);
                *vec2 = *vec1;
              }
            break;

            case 3:
              {
                // OP = (A - sigma I)^(-1) , B = I
                // OP = (A - sigma M)^(-1) M , B = M  (for generalized problem)
                switch (problem_type_enum_ID)
                  {
                  case HERMITIAN_EIGENPROBLEM_ENUM_ID:
                  case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
                    this->linear_solver->solve(*vec1, *vec2);
                    break;
                    
                  case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
                  case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
                    {
                      if (this->ido == -1)
                        {
                        solver_B_mat->multiply_vector(*vec1, *vec2);
                        vec1->zero();
                        this->linear_solver->solve(*vec2, *vec1);
                        *vec2 = *vec1;
                        }
                      else 
                        {
                        for (unsigned int i=0; i<solver_A_mat->m(); i++)
                          vec1->set(i, workd[ipntr[3]+i]);
                        this->linear_solver->solve(*vec1, *vec2);
                        }
                    }
                    break;
                    
                  default:
                    AssertThrow(false, ExcInternalError());
                    break;
                  }
              }
            break;

            case 4: // not handled for now
            case 5: // not handled for now
            default:
              Assert(false, ExcInternalError());
            }
        }
        break;
        
      case 2:
        {
          // compute  Y = B * X 
          switch (this->iparam[7])
            {
            case 1:
              {
                // OP = A, B = I
                *vec2 = *vec1;
              }
              break;
              
            case 2:
              {
                // OP = M^(-1) A , B = M
                solver_B_mat->multiply_vector(*vec1, *vec2);
              }
              break;
              
            case 3:
              {
                // OP = (A - sigma I)^(-1) , B = I
                // OP = (A - sigma M)^(-1) M , B = M  (for generalized problem)
                switch (problem_type_enum_ID)
                  {
                  case HERMITIAN_EIGENPROBLEM_ENUM_ID:
                  case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
                    *vec2 = *vec1;
                    break;
                    
                  case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
                  case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
                    solver_B_mat->multiply_vector(*vec1, *vec2);
                    break;
                    
                  default:
                    AssertThrow(false, ExcInternalError());
                    break;
                  }
              }
              break;
              
            case 4: // not handled for now
            case 5: // not handled for now
            default:
              Assert(false, ExcInternalError());
            }
          
        }
        break;
        
      case 3:
        {
          Assert (false, ExcInternalError());
        }
        break;

      case 99:
        break;
        
      default:
        Assert (false, ExcInternalError());
      }
    
    // after the matrix operations, place the vector for the solver to access
    for (unsigned int i=0; i < solver_A_mat->m(); i++)
      this->workd[ipntr[2]+i] = (*vec2)(i);
    }
    

  // once the solution is done, calculate the eigenvectors and eigenvalues
  switch (problem_type_enum_ID)
  {
          // this is for real symmetric
      case HERMITIAN_EIGENPROBLEM_ENUM_ID:
      case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
          dseupd_(&rvec, const_cast<char*>(HowMny.c_str()), 
                  &select[1], &dr[1],
                  &Z[1], &ldz, 
                  &sigmar, const_cast<char*>(bmat.c_str()),
                  &n, const_cast<char*>(which.c_str()),
                  &nev, &tol, 
                  &resid[1], &ncv, 
                  &V[1], &ldv, 
                  &iparam[1],
                  &ipntr[1], &workd[1], 
                  &workl[1], &lworkl, &info );
      break;
      
          // this is for real un-symmetric
      case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    {
      // the Z vector is overwritten with the first NEV+1 vectors of 
      // V, since we need the ritz vectors of the problem, instead of the
      // the Arnoldi basis
//      for (int j=0; j<(this->nev+1); j++)
//        for (int i=0; i < this->n; i++)
//          this->Z[ j*this->n + i + 1] = this->V[ j*this->n + i + 1];
      
      FC_FUNC(dneupd, DNEUPD)(&rvec, const_cast<char*>(HowMny.c_str()), 
                              &select[1], &dr[1], &di[1],
                              &Z[1], &ldz, &sigmar, &sigmai, &workev[1],
                              const_cast<char*>(bmat.c_str()),
                              &n, const_cast<char*>(which.c_str()),
                              &nev, &tol, 
                              &resid[1], &ncv, 
                              &V[1], &ldv, 
                              &iparam[1],
                              &ipntr[1], &workd[1], 
                              &workl[1], &lworkl, &info );
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
  }

  
  // check the error, if any
  this->checkError(this->info);

  // if everything is done and is complete, then the solution will be available 
  // in the vectors
  this->solution_completed = true;

  this->analysis_driver->getFESystemController().performance_logging->unsetEvent
    ("ArpackEigenSolver::solve()", "Solver");
}



void 
Solver::ArpackEigenSolver::getEigenValue(unsigned int eigen_index,
                                         double* real_value,
                                         double* img_value)
{
  Assert (this->solution_completed, ExcInvalidState());
  // make sure the eigen value being asked for was calculated
  Assert (eigen_index < this->eigen_solver_info.getNEigenPairsToSolve(), ExcInternalError());
  Assert ((int) eigen_index < this->iparam[5], ExcInternalError());
  
  
  *real_value = this->dr[eigen_index+1];
  *img_value = this->di[eigen_index+1];
  if (this->eigen_solver_info.ifSwapMatrices())
    {
      double magnitude = (*real_value)*(*real_value) + (*img_value)*(*img_value);
      (*real_value) /= magnitude;
      (*img_value) /= magnitude;
    }
}



void 
Solver::ArpackEigenSolver::getEigenPair(unsigned int pair_index,
                                        double* real_value,
                                        double* img_value, 
                                        NumericVector<double>& real_vec,
                                        NumericVector<double>& img_vec)
{
  Assert (this->solution_completed, ExcInvalidState());
  // make sure the eigen value being asked for was calculated
  Assert (pair_index < this->eigen_solver_info.getNEigenPairsToSolve(), ExcInternalError());
  Assert ((int) pair_index < this->iparam[5], ExcInternalError());
  
  *real_value = this->dr[pair_index+1];
  *img_value = this->di[pair_index+1];
  if (this->eigen_solver_info.ifSwapMatrices())
    {
      double magnitude = (*real_value)*(*real_value) + (*img_value)*(*img_value);
      (*real_value) /= magnitude;
      (*img_value) /= magnitude;
    }
  
  // next, copy the vectors
  // make sure that the sizes are appropriate
  Assert ((int) real_vec.size() == this->n, ExcInternalError());
  Assert ((int) img_vec.size() == this->n, ExcInternalError());

  unsigned int problem_type_enum_ID = this->eigen_solver_info.getEigenProblemKindEnumID();

  switch (problem_type_enum_ID)
  {
      // this is for real symmetric
    case HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    {
      img_vec.zero();  
      for (int i=0; i < this->n; i++)
          real_vec.set(i, this->Z[(pair_index)*this->n+i+1]);
    }
      break;
      
      // this is for real un-symmetric
    case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID:
    {
      unsigned int index = 0;
      if (pair_index % 2 ==0) 
        index = pair_index / 2;
      else 
        index = (pair_index-1) / 2;
        
      for (int i=0; i < this->n; i++)
        {
          real_vec.set(i, this->V[(2*index) * this->n + i + 1]);
          img_vec.set(i, this->V[(2*index+1) * this->n + i + 1]);
        }
    }
      break;
      
    default:
      AssertThrow(false, ExcInternalError());
      break;
  }
}


double
Solver::ArpackEigenSolver::getResidualForEigenPair(const unsigned int eigen_index)
{
  Assert (this->solution_completed, ExcInvalidState());
  // make sure the eigen value being asked for was calculated
  Assert (eigen_index < this->eigen_solver_info.getNEigenPairsToSolve(), ExcInternalError());
  Assert ((int) eigen_index < this->iparam[5], ExcInternalError());
  
  // to be implemented
  Assert(false, ExcInternalError());
  return 0.0;
}




void 
Solver::ArpackEigenSolver::checkError(const unsigned int ierr)
{
  switch (ierr)
    {
    case 0:
      {
        //    =  0   : Normal exit.
        // nothing to be done, just return
      }
      break;
      
    case 1:
      {
        //    =  1   : Maximum number of iterations taken. All possible
        //    eigenvalues of OP has been found. iparam[5]
        //    returns the number of converged Ritz values.
        std::cout << "Eigensolver reached maximum iterations." << std::endl;
        std::cout << "Number of converged eigenvalues: "
        << this->iparam[5] << std::endl;
      }
      break;

    case 3:
      {
        //    =  3   : No shifts could be applied during a cycle of the
        //    Implicitly restarted Arnoldi iteration. One
        //    possibility is to increase the size of NCV relative
        //    to nev. See remark 4 below.
        Assert(false, ExcInternalError());
      }
      break;

    case -1:
      {
        //    = -1   : n must be positive.
        Assert(false, ExcInternalError());
      }
      break;
      
    case -2:
      {
        //    = -2   : nev must be positive.
        Assert(false, ExcInternalError());
      }
      break;

    case -3:
      {
        //    = -3   : ncv must satisfy nev < ncv <= n.
        Assert(false, ExcInternalError());
      }
      break;

    case -4:
      {
        //    = -4   : The maximum number of Arnoldi update iterations allowed
        //    must be greater than zero.
        Assert(false, ExcInternalError());
      }
      break;

    case -5:
      {
        //    = -5   : which must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
        Assert(false, ExcInternalError());
      }
      break;

    case -6:
      {
        //    = -6   : bmat must be one of 'I' or 'G'.
        Assert(false, ExcInternalError());
      }
      break;

    case -7:
      {
        //    = -7   : Length of private work array workl is not sufficient.
        Assert(false, ExcInternalError());
      }
      break;

    case -8:
      {
        //    = -8   : Error return from trid. eigenvalue calculation;
        //    Informational error from LAPACK routine dsteqr.
        Assert(false, ExcInternalError());
      }
      break;

    case -9:
      {
        //      = -9   : Starting vector is zero.
        Assert(false, ExcInternalError());
      }
      break;

    case -10:
      {
        //      = -10  : iparam[7] must be 1,2,3,4,5.
        Assert(false, ExcInternalError());
      }
      break;

    case -11:
      {
        //      = -11  : iparam[7] = 1 and bmat = 'G' are incompatible.
        Assert(false, ExcInternalError());
      }
      break;

    case -12:
      {
        //      = -12  : iparam[1] must be equal to 0 or 1.
        //    = -12: nev and which = 'BE' are incompatible.
        Assert(false, ExcInternalError());
      }
      break;

    case -13:
      {
        //      = -13  : nev and which = 'BE' are incompatible.
        Assert(false, ExcInternalError());
      }
      break;

    case -14:
      {
        //    = -14: dsaupp did not find any eigenvalues to sufficient
        //    accuracy.
        Assert(false, ExcInternalError());
      }
      break;

    case -15:
      {
        //    = -15: HowMny must be one of 'A' or 'S' if rvec = true.
        Assert(false, ExcInternalError());
      }
      break;

    case -16:
      {
        //      = -16: HowMny = 'S' not yet implemented.
        Assert(false, ExcInternalError());
      }
      break;

    case -9999:
      {
        //      = -9999: Could not build an Arnoldi factorization. iparam[5]
        //      returns the size of the current Arnoldi factorization.
        //      The user is advised to check that enough workspace
        //      and array storage has been allocated.
        Assert(false, ExcInternalError());
      }
      break;
      
      
    default:
      Assert(false, ExcInternalError());
    }
}

