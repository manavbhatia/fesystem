// $Id: NonlinearSolver.h,v 1.6.4.3 2007-05-14 16:45:07 manav Exp $

#ifndef __fesystem_nonlinear_solver_h__
#define __fesystem_nonlinear_solver_h__


// FESystem includes
#include "Solvers/FESystemSolverBase.h"
#include "AnalysisDriver/NonLinearAnalysisDriver.h"

// Forward Declerations
template <typename DataType> class NumericVector;
template <typename DataType> class SparseMatrix;


#ifndef NONLINEAR_SOLVER_ENUM_ID
#define NONLINEAR_SOLVER_ENUM_ID 2
#else
#error
#endif

#ifndef NONLINEAR_SOLVER_ENUM_NAME
#define NONLINEAR_SOLVER_ENUM_NAME "NONLINEAR_SOLVER"
#else
#error
#endif


/// this class provides an interface to a solver for the solution of 
/// a linear system of equations. This inherits from the FESystemSolverBase
/// class. It takes a pointer to the analysis driver object that owns 
/// an instantiation of this class
namespace Solver
{
  
  // Forward decleration
  class LinearSolverInfo;
  class NonlinearSolverInfo;
  
  DeclareEnumClass(NonlinearSolverKindEnum);
  
  
  DeclareEnumName(NONLINEAR_SOLVER, Solver::SolverClassEnum,
                  NONLINEAR_SOLVER_ENUM_ID,
                  NONLINEAR_SOLVER_ENUM_NAME);
  
  class NonlinearSolver : public Solver::FESystemSolverBase
    {
public:
      /// constructor
      NonlinearSolver(const Solver::NonlinearSolverInfo& nonlinear_info,
                      const Solver::LinearSolverInfo& linear_info);
      
      
      /// destructor
      ~NonlinearSolver();
      
      /// method to clear the data structures of this class
      virtual void clear();
      
      inline unsigned int getNonlinearSolverKindEnumID() const;
      
      inline std::string getNonlinearSolverKindEnumName() const;
      
      
      inline void attachMatrixAndVector(SparseMatrix<double>& joc_mat, 
                                        NumericVector<double>& res_vec);
      
      /// function to perform solution steps.
      /// @param solution vector in which the solution will be returned. 
      /// This should be initialized to the 
      /// initial guess
      virtual void solve(NumericVector<double>& solution) = 0;
      
      /// this returns the itertion number of the solution process
      virtual unsigned int getCurrentIterationNumber() const = 0;
      
      /// evaluates the residual vector at the specified solution
      /// @param sol vector at which the residual should be evaluated
      /// @param res vector in which the residual is stored
      void evaluateResidual(NumericVector<double>& sol,
                                   NumericVector<double>& res);
      
      /// evaluates the residual vector at the specified solution
      /// @param sol vector at which the residual should be evaluated
      /// @param jac matrix in which the jacobian is stored
      void evaluateJacobian(NumericVector<double>& sol,
                                   SparseMatrix<double>& jac);
      
      /// this function can be used to set the functions that are called by 
      /// the nonlinear solver to calculate the jacobian and the residual. 
      /// if this is not done, then the solver uses the default functions 
      /// that ask the analysis driver for the information. 
      void setFunctions(void (*res_func)(NumericVector<double>& sol,
                                         NumericVector<double>& res,
                                         void* ctx),
                        void (*jac_func)(NumericVector<double>& sol,
                                         SparseMatrix<double>& jac,
                                         void* ctx),
                        void* ctx);
      
protected:
              
      const unsigned int nonlinear_solver_kind_enum_ID;
      const Solver::NonlinearSolverInfo& nonlinear_solver_info;
      const Solver::LinearSolverInfo& linear_solver_info;
      
      SparseMatrix<double>* jacobian;
      
      NumericVector<double>* residual;
      
      void (* jacobian_function) (NumericVector<double>& sol,
                                  SparseMatrix<double>& jac,
                                  void* ctx);

      void (* residual_function) (NumericVector<double>& sol,
                                  NumericVector<double>& res,
                                  void* ctx);
      
      void * context;
};
}





inline 
unsigned int
Solver::NonlinearSolver::getNonlinearSolverKindEnumID() const
{
  return this->nonlinear_solver_kind_enum_ID;
}






inline
std::string
Solver::NonlinearSolver::getNonlinearSolverKindEnumName() const
{
  return Solver::NonlinearSolverKindEnum::enumName(this->nonlinear_solver_kind_enum_ID);
}




inline
void Solver::NonlinearSolver::attachMatrixAndVector(SparseMatrix<double>& jac_mat, 
                                                    NumericVector<double>& res_vec)
{
  this->jacobian = &jac_mat;
  this->residual = &res_vec;
}



#endif // __fesystem_nonlinear_solver_h__
