// $Id: SlepcEigenSolver.h,v 1.1.4.3 2007-03-14 22:05:03 manav Exp $

#ifndef __fesystem_slepc_eigen_solver_h__
#define __fesystem_slepc_eigen_solver_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/EigenSolver.h"

// SLEPc includes
#include "slepceps.h"

// Forward Declerations
template <typename DataType> class NumericVector;
template <typename DataType> class SparseMatrix;


#ifndef EPS_LAPACK_EIGEN_SOLVER_ENUM_ID
#define EPS_LAPACK_EIGEN_SOLVER_ENUM_ID 1
#else
#error
#endif

#ifndef EPS_LAPACK_EIGEN_SOLVER_ENUM_NAME
#define EPS_LAPACK_EIGEN_SOLVER_ENUM_NAME "EPS_LAPACK_EIGEN_SOLVER"
#else
#error
#endif


#ifndef EPS_POWER_EIGEN_SOLVER_ENUM_ID
#define EPS_POWER_EIGEN_SOLVER_ENUM_ID 2
#else
#error
#endif

#ifndef EPS_POWER_EIGEN_SOLVER_ENUM_NAME
#define EPS_POWER_EIGEN_SOLVER_ENUM_NAME "EPS_POWER_EIGEN_SOLVER"
#else
#error
#endif


#ifndef EPS_SUBSPACE_EIGEN_SOLVER_ENUM_ID
#define EPS_SUBSPACE_EIGEN_SOLVER_ENUM_ID 3
#else
#error
#endif

#ifndef EPS_SUBSPACE_EIGEN_SOLVER_ENUM_NAME
#define EPS_SUBSPACE_EIGEN_SOLVER_ENUM_NAME "EPS_SUBSPACE_EIGEN_SOLVER"
#else
#error
#endif


#ifndef EPS_ARNOLDI_EIGEN_SOLVER_ENUM_ID
#define EPS_ARNOLDI_EIGEN_SOLVER_ENUM_ID 4
#else
#error
#endif

#ifndef EPS_ARNOLDI_EIGEN_SOLVER_ENUM_NAME
#define EPS_ARNOLDI_EIGEN_SOLVER_ENUM_NAME "EPS_ARNOLDI_EIGEN_SOLVER"
#else
#error
#endif


#ifndef EPS_LANCZOS_EIGEN_SOLVER_ENUM_ID
#define EPS_LANCZOS_EIGEN_SOLVER_ENUM_ID 5
#else
#error
#endif

#ifndef EPS_LANCZOS_EIGEN_SOLVER_ENUM_NAME
#define EPS_LANCZOS_EIGEN_SOLVER_ENUM_NAME "EPS_LANCZOS_EIGEN_SOLVER"
#else
#error
#endif


#ifndef EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_ID
#define EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_ID 6
#else
#error
#endif

#ifndef EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_NAME
#define EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_NAME "EPS_KRYLOV_SCHUR_EIGEN_SOLVER"
#else
#error
#endif




namespace Solver
{
  
  // Forward declerations
  class PetscLinearSolver;
      
  DeclareEnumName(EPS_LAPACK_EIGEN_SOLVER, Solver::EigenSolverKindEnum,
                  EPS_LAPACK_EIGEN_SOLVER_ENUM_ID,
                  EPS_LAPACK_EIGEN_SOLVER_ENUM_NAME);
  
  DeclareEnumName(EPS_POWER_EIGEN_SOLVER, Solver::EigenSolverKindEnum,
                  EPS_POWER_EIGEN_SOLVER_ENUM_ID,
                  EPS_POWER_EIGEN_SOLVER_ENUM_NAME);
  
  DeclareEnumName(EPS_SUBSPACE_EIGEN_SOLVER, Solver::EigenSolverKindEnum,
                  EPS_SUBSPACE_EIGEN_SOLVER_ENUM_ID,
                  EPS_SUBSPACE_EIGEN_SOLVER_ENUM_NAME);
  
  DeclareEnumName(EPS_ARNOLDI_EIGEN_SOLVER, Solver::EigenSolverKindEnum,
                  EPS_ARNOLDI_EIGEN_SOLVER_ENUM_ID,
                  EPS_ARNOLDI_EIGEN_SOLVER_ENUM_NAME);
  
  DeclareEnumName(EPS_LANCZOS_EIGEN_SOLVER, Solver::EigenSolverKindEnum,
                  EPS_LANCZOS_EIGEN_SOLVER_ENUM_ID,
                  EPS_LANCZOS_EIGEN_SOLVER_ENUM_NAME);
  
  DeclareEnumName(EPS_KRYLOV_SCHUR_EIGEN_SOLVER, Solver::EigenSolverKindEnum,
                  EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_ID,
                  EPS_KRYLOV_SCHUR_EIGEN_SOLVER_ENUM_NAME);

  
  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class SlepcEigenSolver : public EigenSolverBase
  {
public:

    /// constructor
    SlepcEigenSolver(const Solver::EigenSolverInfo& eigen_info,
                const Solver::LinearSolverInfo& linear_info);
    
    // destructor
    virtual ~SlepcEigenSolver();
    
    /// this method clears the data structures of this object. This should be called 
    /// each time the used finishes using this object.
    virtual void clear();
    
    
    /// @returns the number of converged eigen pairs
    virtual unsigned int getNConvergedEigenPairs();
    


    /// this method returns the eigen valie
    /// @param pair_index index of the requested eigenvalue
    /// @param real_value real part of the eigenvalue
    /// @param img_value imaginary part of the eigenvalue
    virtual void getEigenValue(unsigned int eigen_index,
                      double* real_value,
                      double* img_value);
    
    
    /// this method returns the eigen pair 
    /// @param pair_index index of the requested eigenpair
    /// @param real_value real part of the eigenvalue
    /// @param img_value imaginary part of the eigenvalue
    /// @param real_vec real part of the eigenvector
    /// @param img_vec imaginary part of the eigenvector
    virtual void getEigenPair(unsigned int pair_index,
                      double* real_value,
                      double* img_value, 
                      NumericVector<double>& real_vec,
                      NumericVector<double>& img_vec);

    virtual void getInvariantSubspace(std::vector<NumericVector<double>*>& vectors);

    /// this method computes and returns the residual error for the specified 
    /// eigen value index. This is only for Hermitian problems
    virtual double getEigenPairResidualError(const unsigned int index)
    {
      (void) index;
      Assert(false, ExcInternalError());
      return (double) NULL;
    }

    
    /// this method computes and returns the relative error for the specified eigen value 
    /// index. This is only for Hermitian problems
    virtual double getEigenPairRelativeError(const unsigned int index)
    {
      (void) index;
      Assert(false, ExcInternalError());
      return (double) NULL;
    }

    
    /// method to solve the eigen system
    virtual void solve(SparseMatrix<double>* A_mat, SparseMatrix<double>* B_mat = NULL);

    /// this method calculates a residual for the given eigen value and returns it
    virtual double getResidualForEigenPair(const unsigned int i);
    
protected:
    
      /// sets the solver options
      void setSolverOptions();
    
      /// sets the problem type 
      void setSolverType(const unsigned int enum_ID);

    /// this sets the problem kind
      void setProblemKind(const unsigned int enum_ID);
      
      /// sets the number of eigen sets to compute
      void setNEigenSetsToCompute(const unsigned int n_eigs);

      /// sets the eigen spectrum end to solve for
      void setEigenSpectrumEnd(const unsigned int eigen_spectrum_end);
      
      /// sets the shift for the solver
      void setEigenShift(const unsigned int eigen_shift_kind,
                         const double shift_value);

      /// sets the linear solver details
      void setLinearSolverDetails();

    /// if the solution has been performed
    bool solution_completed;
      
    /// number of converged eigenpairs
    int n_converged_eigen_pairs;

    /// an eigenvalue problem solver context
    EPS eps;
    
    /// which eigenspectrum end to solve for
    EPSWhich eps_spectrum_end;
    
    /// vector of sorted eigenvalue indexes, for permutation
    std::vector<int> eigen_val_permutation;
    
    std::auto_ptr<Solver::PetscLinearSolver> petsc_linear_solver;
  };
  
}


#endif // __fesystem_slepc_eigen_solver_h__
