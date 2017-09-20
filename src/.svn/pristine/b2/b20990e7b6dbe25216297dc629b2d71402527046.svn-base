// $Id: EigenSolver.h,v 1.11.4.2 2007-05-11 05:16:54 manav Exp $

#ifndef __fesystem_eigen_solver_h__
#define __fesystem_eigen_solver_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Solvers/FESystemSolverBase.h"

// Forward Declerations
template <typename DataType> class NumericVector;
template <typename DataType> class SparseMatrix;

namespace Solver
{
  class LinearSolverInfo;
  class EigenSolverInfo;
}



#ifndef EIGEN_SOLVER_ENUM_ID
#define EIGEN_SOLVER_ENUM_ID 3
#else
#error
#endif

#ifndef EIGEN_SOLVER_ENUM_NAME
#define EIGEN_SOLVER_ENUM_NAME "EIGEN_SOLVER"
#else
#error
#endif



#ifndef HERMITIAN_EIGENPROBLEM_ENUM_ID
#define HERMITIAN_EIGENPROBLEM_ENUM_ID 1
#else
#error
#endif

#ifndef HERMITIAN_EIGENPROBLEM_ENUM_NAME
#define HERMITIAN_EIGENPROBLEM_ENUM_NAME "HERMITIAN_EIGENPROBLEM"
#else
#error
#endif


#ifndef GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID
#define GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID 2
#else
#error
#endif

#ifndef GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_NAME
#define GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_NAME "GENERALIZED_HERMITIAN_EIGENPROBLEM"
#else
#error
#endif


#ifndef NON_HERMITIAN_EIGENPROBLEM_ENUM_ID
#define NON_HERMITIAN_EIGENPROBLEM_ENUM_ID 3
#else
#error
#endif

#ifndef NON_HERMITIAN_EIGENPROBLEM_ENUM_NAME
#define NON_HERMITIAN_EIGENPROBLEM_ENUM_NAME "NON_HERMITIAN_EIGENPROBLEM"
#else
#error
#endif


#ifndef GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID
#define GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID 4
#else
#error
#endif

#ifndef GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_NAME
#define GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_NAME "GENERALIZED_NON_HERMITIAN_EIGENPROBLEM"
#else
#error
#endif



#ifndef LARGEST_MAGNITUDE_ENUM_ID 
#define LARGEST_MAGNITUDE_ENUM_ID 1
#else
#error
#endif

#ifndef LARGEST_MAGNITUDE_ENUM_NAME 
#define LARGEST_MAGNITUDE_ENUM_NAME "LARGEST_MAGNITUDE"
#else
#error
#endif

#ifndef SMALLEST_MAGNITUDE_ENUM_ID 
#define SMALLEST_MAGNITUDE_ENUM_ID 2
#else
#error
#endif

#ifndef SMALLEST_MAGNITUDE_ENUM_NAME 
#define SMALLEST_MAGNITUDE_ENUM_NAME "SMALLEST_MAGNITUDE"
#else
#error
#endif


#ifndef LARGEST_REAL_ENUM_ID 
#define LARGEST_REAL_ENUM_ID 3
#else
#error
#endif

#ifndef LARGEST_REAL_ENUM_NAME 
#define LARGEST_REAL_ENUM_NAME "LARGEST_REAL"
#else
#error
#endif

#ifndef SMALLEST_REAL_ENUM_ID 
#define SMALLEST_REAL_ENUM_ID 4
#else
#error
#endif

#ifndef SMALLEST_REAL_ENUM_NAME 
#define SMALLEST_REAL_ENUM_NAME "SMALLEST_REAL"
#else
#error
#endif



#ifndef LARGEST_IMAGINARY_ENUM_ID 
#define LARGEST_IMAGINARY_ENUM_ID 5
#else
#error
#endif

#ifndef LARGEST_IMAGINARY_ENUM_NAME 
#define LARGEST_IMAGINARY_ENUM_NAME "LARGEST_IMAGINARY"
#else
#error
#endif

#ifndef SMALLEST_IMAGINARY_ENUM_ID 
#define SMALLEST_IMAGINARY_ENUM_ID 6
#else
#error
#endif

#ifndef SMALLEST_IMAGINARY_ENUM_NAME 
#define SMALLEST_IMAGINARY_ENUM_NAME "SMALLEST_IMAGINARY"
#else
#error
#endif


#ifndef NO_SHIFT_ENUM_ID 
#define NO_SHIFT_ENUM_ID 1
#else
#error
#endif

#ifndef NO_SHIFT_ENUM_NAME 
#define NO_SHIFT_ENUM_NAME "NO_SHIFT"
#else
#error
#endif


#ifndef ORIGIN_SHIFT_ENUM_ID 
#define ORIGIN_SHIFT_ENUM_ID 2
#else
#error
#endif

#ifndef ORIGIN_SHIFT_ENUM_NAME 
#define ORIGIN_SHIFT_ENUM_NAME "ORIGIN_SHIFT"
#else
#error
#endif


#ifndef SPECTRUM_FOLD_ENUM_ID 
#define SPECTRUM_FOLD_ENUM_ID 3
#else
#error
#endif

#ifndef SPECTRUM_FOLD_ENUM_NAME 
#define SPECTRUM_FOLD_ENUM_NAME "SPECTRUM_FOLD"
#else
#error
#endif


#ifndef SHIFT_AND_INVERT_ENUM_ID 
#define SHIFT_AND_INVERT_ENUM_ID 4
#else
#error
#endif

#ifndef SHIFT_AND_INVERT_ENUM_NAME 
#define SHIFT_AND_INVERT_ENUM_NAME "SHIFT_AND_INVERT"
#else
#error
#endif



#ifndef CAYLEY_SHIFT_ENUM_ID 
#define CAYLEY_SHIFT_ENUM_ID 5
#else
#error
#endif

#ifndef CAYLEY_SHIFT_ENUM_NAME 
#define CAYLEY_SHIFT_ENUM_NAME "CAYLEY_SHIFT"
#else
#error
#endif



namespace Solver
{
  DeclareEnumName(EIGEN_SOLVER, Solver::SolverClassEnum,
                  EIGEN_SOLVER_ENUM_ID,
                  EIGEN_SOLVER_ENUM_NAME);
  
  DeclareEnumClass(EigenProblemKindEnum);

  DeclareEnumName(HERMITIAN_EIGENPROBLEM, Solver::EigenProblemKindEnum,
                  HERMITIAN_EIGENPROBLEM_ENUM_ID,
                  HERMITIAN_EIGENPROBLEM_ENUM_NAME);

  DeclareEnumName(GENERALIZED_HERMITIAN_EIGENPROBLEM, Solver::EigenProblemKindEnum,
                  GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID,
                  GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_NAME);

  DeclareEnumName(NON_HERMITIAN_EIGENPROBLEM, Solver::EigenProblemKindEnum,
                  NON_HERMITIAN_EIGENPROBLEM_ENUM_ID,
                  NON_HERMITIAN_EIGENPROBLEM_ENUM_NAME);

  DeclareEnumName(GENERALIZED_NON_HERMITIAN_EIGENPROBLEM, Solver::EigenProblemKindEnum,
                  GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID,
                  GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_NAME);
  
  
  DeclareEnumClass(EigenSolutionSpectrumKindEnum);
  
  DeclareEnumName(LARGEST_MAGNITUDE, Solver::EigenSolutionSpectrumKindEnum,
                  LARGEST_MAGNITUDE_ENUM_ID,
                  LARGEST_MAGNITUDE_ENUM_NAME);
  
  DeclareEnumName(SMALLEST_MAGNITUDE, Solver::EigenSolutionSpectrumKindEnum,
                  SMALLEST_MAGNITUDE_ENUM_ID,
		  SMALLEST_MAGNITUDE_ENUM_NAME);
  
  DeclareEnumName(LARGEST_REAL, Solver::EigenSolutionSpectrumKindEnum,
                  LARGEST_REAL_ENUM_ID,
                  LARGEST_REAL_ENUM_NAME);

  DeclareEnumName(SMALLEST_REAL, Solver::EigenSolutionSpectrumKindEnum,
                  SMALLEST_REAL_ENUM_ID,
		  SMALLEST_REAL_ENUM_NAME);

  DeclareEnumName(LARGEST_IMAGINARY, Solver::EigenSolutionSpectrumKindEnum,
                  LARGEST_IMAGINARY_ENUM_ID,
                  LARGEST_IMAGINARY_ENUM_NAME);

  DeclareEnumName(SMALLEST_IMAGINARY, Solver::EigenSolutionSpectrumKindEnum,
                  SMALLEST_IMAGINARY_ENUM_ID,
		  SMALLEST_IMAGINARY_ENUM_NAME);
  
  
  
  DeclareEnumClass(EigenSolutionShiftKindEnum);
  
  
  DeclareEnumName(NO_SHIFT, Solver::EigenSolutionShiftKindEnum,
                  NO_SHIFT_ENUM_ID,
                  NO_SHIFT_ENUM_NAME);

  
  DeclareEnumName(ORIGIN_SHIFT, Solver::EigenSolutionShiftKindEnum,
                  ORIGIN_SHIFT_ENUM_ID,
                  ORIGIN_SHIFT_ENUM_NAME);

  DeclareEnumName(SPECTRUM_FOLD, Solver::EigenSolutionShiftKindEnum,
                  SPECTRUM_FOLD_ENUM_ID,
                  SPECTRUM_FOLD_ENUM_NAME);
  
  DeclareEnumName(SHIFT_AND_INVERT, Solver::EigenSolutionShiftKindEnum,
                  SHIFT_AND_INVERT_ENUM_ID,
                  SHIFT_AND_INVERT_ENUM_NAME);

  DeclareEnumName(CAYLEY_SHIFT, Solver::EigenSolutionShiftKindEnum,
                  CAYLEY_SHIFT_ENUM_ID,
                  CAYLEY_SHIFT_ENUM_NAME);

  DeclareEnumClass(EigenSolverKindEnum);

  /// this class provides an interface to a solver for the solution of 
  /// eigensystem. This inherits from the FESystemSolverBase
  /// class. It takes a pointer to the analysis driver object that owns 
  /// an instantiation of this class
  class EigenSolverBase : public FESystemSolverBase
  {
public:

    /// constructor
    EigenSolverBase(const Solver::EigenSolverInfo& eigen_info,
                    const Solver::LinearSolverInfo& linear_info);
    
    // destructor
    virtual ~EigenSolverBase();
    
    
    /// @returns the number of converged eigen pairs
    virtual unsigned int getNConvergedEigenPairs() =0;
    


    /// this method returns the eigen valie
    /// @param pair_index index of the requested eigenvalue
    /// @param real_value real part of the eigenvalue
    /// @param img_value imaginary part of the eigenvalue
    virtual void getEigenValue(unsigned int eigen_index,
                               double* real_value,
                               double* img_value) = 0;
    
    
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
                              NumericVector<double>& img_vec) = 0;

    virtual void getInvariantSubspace(std::vector<NumericVector<double>*>& vectors) = 0;

    
    /// this method computes and returns the residual error for the specified 
    /// eigen value index. This is only for Hermitian problems
    virtual double getEigenPairResidualError(const unsigned int index) = 0;

    
    /// this method computes and returns the relative error for the specified eigen value 
    /// index. This is only for Hermitian problems
    virtual double getEigenPairRelativeError(const unsigned int index) = 0;

    
    /// method to solve the eigen system
    virtual void solve(SparseMatrix<double>* A_mat, SparseMatrix<double>* B_mat = NULL) = 0;

    /// this method calculates a residual for the given eigen value and returns it
    virtual double getResidualForEigenPair(const unsigned int i) = 0;
    
protected:
    
        /// EigenSolverInfo that stores information about this solver
    const Solver::EigenSolverInfo& eigen_solver_info;

    /// LinearSolverInfo that stores information about this solver
    const Solver::LinearSolverInfo& linear_solver_info;
  };
  
}


#endif // __fesystem_eigen_solver_h__
