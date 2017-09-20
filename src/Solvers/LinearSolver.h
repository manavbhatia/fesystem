// $Id: LinearSolver.h,v 1.5.6.3 2007-05-11 05:16:54 manav Exp $

#ifndef __fesystem_linear_solver_h__
#define __fesystem_linear_solver_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Solvers/FESystemSolverBase.h"

// petsc includes
#include "petscksp.h"

// Forward Declerations
template <typename DataType> class NumericVector;
template <typename DataType> class SparseMatrix;




#ifndef LINEAR_SOLVER_ENUM_ID
#define LINEAR_SOLVER_ENUM_ID 1
#else
#error
#endif

#ifndef LINEAR_SOLVER_ENUM_NAME
#define LINEAR_SOLVER_ENUM_NAME "LINEAR_SOLVER"
#else
#error
#endif


#ifndef PC_LU_ENUM_ID
#define PC_LU_ENUM_ID 1
#else
#error
#endif

#ifndef PC_LU_ENUM_NAME
#define PC_LU_ENUM_NAME "PC_LU"
#else
#error
#endif


#ifndef PC_JACOBI_ENUM_ID
#define PC_JACOBI_ENUM_ID 2
#else
#error
#endif

#ifndef PC_JACOBI_ENUM_NAME
#define PC_JACOBI_ENUM_NAME "PC_JACOBI"
#else
#error
#endif


#ifndef PC_BJACOBI_ENUM_ID
#define PC_BJACOBI_ENUM_ID 3
#else
#error
#endif

#ifndef PC_BJACOBI_ENUM_NAME
#define PC_BJACOBI_ENUM_NAME "PC_BJACOBI"
#else
#error
#endif


#ifndef PC_SOR_ENUM_ID
#define PC_SOR_ENUM_ID 4
#else
#error
#endif

#ifndef PC_SOR_ENUM_NAME
#define PC_SOR_ENUM_NAME "PC_SOR"
#else
#error
#endif


#ifndef PC_EISENSTAT_ENUM_ID
#define PC_EISENSTAT_ENUM_ID 5
#else
#error
#endif

#ifndef PC_EISENSTAT_ENUM_NAME
#define PC_EISENSTAT_ENUM_NAME "PC_EISENSTAT"
#else
#error
#endif


#ifndef PC_ICC_ENUM_ID
#define PC_ICC_ENUM_ID 6
#else
#error
#endif

#ifndef PC_ICC_ENUM_NAME
#define PC_ICC_ENUM_NAME "PC_ICC"
#else
#error
#endif


#ifndef PC_ILU_ENUM_ID
#define PC_ILU_ENUM_ID 7
#else
#error
#endif

#ifndef PC_ILU_ENUM_NAME
#define PC_ILU_ENUM_NAME "PC_ILU"
#else
#error
#endif


#ifndef PC_ASM_ENUM_ID
#define PC_ASM_ENUM_ID 8
#else
#error
#endif

#ifndef PC_ASM_ENUM_NAME
#define PC_ASM_ENUM_NAME "PC_ASM"
#else
#error
#endif


#ifndef PC_KSP_ENUM_ID
#define PC_KSP_ENUM_ID 9
#else
#error
#endif

#ifndef PC_KSP_ENUM_NAME
#define PC_KSP_ENUM_NAME "PC_KSP"
#else
#error
#endif



#ifndef PC_CHOLESKY_ENUM_ID
#define PC_CHOLESKY_ENUM_ID 10
#else
#error
#endif

#ifndef PC_CHOLESKY_ENUM_NAME
#define PC_CHOLESKY_ENUM_NAME "PC_CHOLESKY"
#else
#error
#endif


#ifndef PC_NONE_ENUM_ID
#define PC_NONE_ENUM_ID 11
#else
#error
#endif

#ifndef PC_NONE_ENUM_NAME
#define PC_NONE_ENUM_NAME "PC_NONE"
#else
#error
#endif


#ifndef KSP_RICHARDSON_ENUM_ID
#define KSP_RICHARDSON_ENUM_ID 1
#else
#error
#endif

#ifndef KSP_RICHARDSON_ENUM_NAME
#define KSP_RICHARDSON_ENUM_NAME "KSP_RICHARDSON"
#else
#error
#endif


#ifndef KSP_CHEBYCHEV_ENUM_ID
#define KSP_CHEBYCHEV_ENUM_ID 2
#else
#error
#endif

#ifndef KSP_CHEBYCHEV_ENUM_NAME
#define KSP_CHEBYCHEV_ENUM_NAME "KSP_CHEBYCHEV"
#else
#error
#endif


#ifndef KSP_CG_ENUM_ID
#define KSP_CG_ENUM_ID 3
#else
#error
#endif

#ifndef KSP_CG_ENUM_NAME
#define KSP_CG_ENUM_NAME "KSP_CG"
#else
#error
#endif


#ifndef KSP_GMRES_ENUM_ID
#define KSP_GMRES_ENUM_ID 4
#else
#error
#endif

#ifndef KSP_GMRES_ENUM_NAME
#define KSP_GMRES_ENUM_NAME "KSP_GMRES"
#else
#error
#endif


#ifndef KSP_TCQMR_ENUM_ID
#define KSP_TCQMR_ENUM_ID 5
#else
#error
#endif

#ifndef KSP_TCQMR_ENUM_NAME
#define KSP_TCQMR_ENUM_NAME "KSP_TCQMR"
#else
#error
#endif


#ifndef KSP_BCGS_ENUM_ID
#define KSP_BCGS_ENUM_ID 6
#else
#error
#endif

#ifndef KSP_BCGS_ENUM_NAME
#define KSP_BCGS_ENUM_NAME "KSP_BCGS"
#else
#error
#endif



#ifndef KSP_CGS_ENUM_ID
#define KSP_CGS_ENUM_ID 7
#else
#error
#endif

#ifndef KSP_CGS_ENUM_NAME
#define KSP_CGS_ENUM_NAME "KSP_CGS"
#else
#error
#endif


#ifndef KSP_TFQMR_ENUM_ID
#define KSP_TFQMR_ENUM_ID 8
#else
#error
#endif

#ifndef KSP_TFQMR_ENUM_NAME
#define KSP_TFQMR_ENUM_NAME "KSP_TFQMR"
#else
#error
#endif


#ifndef KSP_CR_ENUM_ID
#define KSP_CR_ENUM_ID 9
#else
#error
#endif

#ifndef KSP_CR_ENUM_NAME
#define KSP_CR_ENUM_NAME "KSP_CR"
#else
#error
#endif


#ifndef KSP_LSQR_ENUM_ID
#define KSP_LSQR_ENUM_ID 10
#else
#error
#endif

#ifndef KSP_LSQR_ENUM_NAME
#define KSP_LSQR_ENUM_NAME "KSP_LSQR"
#else
#error
#endif


#ifndef KSP_BICG_ENUM_ID
#define KSP_BICG_ENUM_ID 11
#else
#error
#endif

#ifndef KSP_BICG_ENUM_NAME
#define KSP_BICG_ENUM_NAME "KSP_BICG"
#else
#error
#endif

#ifndef KSP_PRE_ONLY_ENUM_ID
#define KSP_PRE_ONLY_ENUM_ID 12
#else
#error
#endif

#ifndef KSP_PRE_ONLY_ENUM_NAME
#define KSP_PRE_ONLY_ENUM_NAME "KSP_PRE_ONLY"
#else
#error
#endif

//#ifndef ATTACH_AND_FACTORIZE_MATRIX_ENUM_ID
//#define ATTACH_AND_FACTORIZE_MATRIX_ENUM_ID 1
//#else
//#error
//#endif
//
//#ifndef ATTACH_AND_FACTORIZE_MATRIX_ENUM_NAME
//#define ATTACH_AND_FACTORIZE_MATRIX_ENUM_NAME "ATTACH_AND_FACTORIZE_MATRIX"
//#else
//#error
//#endif
//
//
//#ifndef ATTACH_FACTORIZED_MATRIX_ENUM_ID
//#define ATTACH_FACTORIZED_MATRIX_ENUM_ID 2
//#else
//#error
//#endif
//
//#ifndef ATTACH_FACTORIZED_MATRIX_ENUM_NAME
//#define ATTACH_FACTORIZED_MATRIX_ENUM_NAME "ATTACH_FACTORIZED_MATRIX"
//#else
//#error
//#endif




/// this class provides an interface to a solver for the solution of 
/// a linear system of equations. This inherits from the FESystemSolverBase
/// class. It takes a pointer to the analysis driver object that owns 
/// an instantiation of this class
namespace Solver
{
  DeclareEnumName(LINEAR_SOLVER, Solver::SolverClassEnum,
                  LINEAR_SOLVER_ENUM_ID,
                  LINEAR_SOLVER_ENUM_NAME);

  DeclareEnumClass(LinearSolverKindEnum);
  
  DeclareEnumClass(PCTypeEnum);

  DeclareEnumName(PC_LU, Solver::PCTypeEnum,
                  PC_LU_ENUM_ID,
                  PC_LU_ENUM_NAME);
  DeclareEnumName(PC_JACOBI, Solver::PCTypeEnum,
                  PC_JACOBI_ENUM_ID,
                  PC_JACOBI_ENUM_NAME);
  DeclareEnumName(PC_BJACOBI, Solver::PCTypeEnum,
                  PC_BJACOBI_ENUM_ID,
                  PC_BJACOBI_ENUM_NAME);
  DeclareEnumName(PC_SOR, Solver::PCTypeEnum,
                  PC_SOR_ENUM_ID,
                  PC_SOR_ENUM_NAME);
  DeclareEnumName(PC_EISENSTAT, Solver::PCTypeEnum,
                  PC_EISENSTAT_ENUM_ID,
                  PC_EISENSTAT_ENUM_NAME);
  DeclareEnumName(PC_ICC, Solver::PCTypeEnum,
                  PC_ICC_ENUM_ID,
                  PC_ICC_ENUM_NAME);
  DeclareEnumName(PC_ILU, Solver::PCTypeEnum,
                  PC_ILU_ENUM_ID,
                  PC_ILU_ENUM_NAME);
  DeclareEnumName(PC_ASM, Solver::PCTypeEnum,
                  PC_ASM_ENUM_ID,
                  PC_ASM_ENUM_NAME);
  DeclareEnumName(PC_KSP, Solver::PCTypeEnum,
                  PC_KSP_ENUM_ID,
                  PC_KSP_ENUM_NAME);
  DeclareEnumName(PC_CHOLESKY, Solver::PCTypeEnum,
                  PC_CHOLESKY_ENUM_ID,
                  PC_CHOLESKY_ENUM_NAME);
  DeclareEnumName(PC_NONE, Solver::PCTypeEnum,
                  PC_NONE_ENUM_ID,
                  PC_NONE_ENUM_NAME);
  
  DeclareEnumClass(KSPTypeEnum);
  
  DeclareEnumName(KSP_RICHARDSON, Solver::KSPTypeEnum,
                  KSP_RICHARDSON_ENUM_ID,
                  KSP_RICHARDSON_ENUM_NAME);
  DeclareEnumName(KSP_CHEBYCHEV, Solver::KSPTypeEnum,
                  KSP_CHEBYCHEV_ENUM_ID,
                  KSP_CHEBYCHEV_ENUM_NAME);
  DeclareEnumName(KSP_CG, Solver::KSPTypeEnum,
                  KSP_CG_ENUM_ID,
                  KSP_CG_ENUM_NAME);
  DeclareEnumName(KSP_GMRES, Solver::KSPTypeEnum,
                  KSP_GMRES_ENUM_ID,
                  KSP_GMRES_ENUM_NAME);
  DeclareEnumName(KSP_TCQMR, Solver::KSPTypeEnum,
                  KSP_TCQMR_ENUM_ID,
                  KSP_TCQMR_ENUM_NAME);
  DeclareEnumName(KSP_BCG, Solver::KSPTypeEnum,
                  KSP_BCGS_ENUM_ID,
                  KSP_BCGS_ENUM_NAME);
  DeclareEnumName(KSP_CGS, Solver::KSPTypeEnum,
                  KSP_CGS_ENUM_ID,
                  KSP_CGS_ENUM_NAME);
  DeclareEnumName(KSP_TFQMR, Solver::KSPTypeEnum,
                  KSP_TFQMR_ENUM_ID,
                  KSP_TFQMR_ENUM_NAME);
  DeclareEnumName(KSP_CR, Solver::KSPTypeEnum,
                  KSP_CR_ENUM_ID,
                  KSP_CR_ENUM_NAME);
  DeclareEnumName(KSP_LSQR, Solver::KSPTypeEnum,
                  KSP_LSQR_ENUM_ID,
                  KSP_LSQR_ENUM_NAME);
  DeclareEnumName(KSP_BICG, Solver::KSPTypeEnum,
                  KSP_BICG_ENUM_ID,
                  KSP_BICG_ENUM_NAME);
  DeclareEnumName(KSP_PRE_ONLY, Solver::KSPTypeEnum,
                  KSP_PRE_ONLY_ENUM_ID,
                  KSP_PRE_ONLY_ENUM_NAME);
  
//  DeclareEnumClass(LinearSolverOptionEnum);
//  
//  
//  DeclareEnumName(ATTACH_AND_FACTORIZE_MATRIX, Solver::LinearSolverOptionEnum,
//                  ATTACH_AND_FACTORIZE_MATRIX_ENUM_ID,
//                  ATTACH_AND_FACTORIZE_MATRIX_ENUM_NAME);
//  
//  DeclareEnumName(ATTACH_FACTORIZED_MATRIX, Solver::LinearSolverOptionEnum,
//                  ATTACH_FACTORIZED_MATRIX_ENUM_ID,
//                  ATTACH_FACTORIZED_MATRIX_ENUM_NAME);
  
  
  
  
  
  class LinearSolver : public FESystemSolverBase
  {
public:
    
    /// constructor
    LinearSolver(const Solver::LinearSolverInfo& info,
                 const unsigned int solver_kind);
    
    
    /// destructor
    ~LinearSolver();
    
    /// method to clear the data structures after use
    virtual void clear();

    inline unsigned int getLinearSolverKindEnumID() const;

    inline std::string getLinearSolverKindEnumName() const;

    inline const Solver::LinearSolverInfo& getSolverInfo() const; 
      
    virtual void setSystemMatrix(SparseMatrix<double>& matrix);

    virtual void setPreconditionerMatrix(SparseMatrix<double>& matrix);

    
    /// this is a function that implements the solution steps. Before this 
    /// function is called, the matrices and vector should already have been 
    /// set.
    /// @param rhs right hand side of the system of equations
    /// @param sol vector in which solution will be stored
    virtual void solve(NumericVector<double>& rhs,
                       NumericVector<double>& sol) = 0;
    
    bool useSameMatricesBetweenSolves();
      
protected:
        
    const unsigned int linear_solver_kind_enum_ID;
    const Solver::LinearSolverInfo& solver_info;
    
    const unsigned int ksp_type_enum_ID;
    const unsigned int pc_type_enum_ID;
    
    /// the preconditioner matrix
    SparseMatrix<double> *system_matrix;

    /// bool to check if the system matrix is new
    bool new_system_matrix;
    
    /// the preconditioner matrix
    SparseMatrix<double> *preconditioner_matrix; 
    
    /// bool to check if the PC matrix is new
    bool new_PC_matrix;
    
    /// flag that the user must set to state if the same marices should
    /// be assumed between solves
    bool same_matrices_between_solves;
    
    /// flag that the user must set to state if the same marices should
    /// be assumed between solves
    unsigned int n_solves_after_new_matrices;
    
  };
}








inline 
unsigned int
Solver::LinearSolver::getLinearSolverKindEnumID() const
{
  return this->linear_solver_kind_enum_ID;
}






inline
std::string
Solver::LinearSolver::getLinearSolverKindEnumName() const
{
  return Solver::LinearSolverKindEnum::enumName(this->linear_solver_kind_enum_ID);
}




inline
const Solver::LinearSolverInfo& 
Solver::LinearSolver::getSolverInfo() const
{
  return this->solver_info;
}



#endif // __fesystem_linear_solver_h__
