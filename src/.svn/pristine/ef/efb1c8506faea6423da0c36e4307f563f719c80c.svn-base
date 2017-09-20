// $Id:$
/*
 *  PetscLinearSolver.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 12/4/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */

// $Id: PetscLinearSolver.h,v 1.5.4.3 2008-08-21 00:56:46 manav Exp $

#ifndef __fesystem_petsc_linear_solver_h__
#define __fesystem_petsc_linear_solver_h__

// FESystem includes
#include "Solvers/LinearSolver.h"
#include "FESystem/FESystemController.h"

// libMesh includes
#include "numerics/petsc_matrix.h"
#include "numerics/petsc_vector.h"




#ifndef PETSC_LINEAR_SOLVER_ENUM_ID
#define PETSC_LINEAR_SOLVER_ENUM_ID 1
#else
#error
#endif

#ifndef PETSC_LINEAR_SOLVER_ENUM_NAME
#define PETSC_LINEAR_SOLVER_ENUM_NAME "PETSC_LINEAR_SOLVER"
#else
#error
#endif


namespace Solver
{
  
  DeclareEnumName(PETSC_LINEAR_SOLVER, LinearSolverKindEnum,
                  PETSC_LINEAR_SOLVER_ENUM_ID,
                  PETSC_LINEAR_SOLVER_ENUM_NAME);
  
  
  class PetscLinearSolver: public LinearSolver
    {
    public:
      PetscLinearSolver(const LinearSolverInfo& info);
      
      ~PetscLinearSolver();
      
      /// method to clear the data structures
      virtual void clear();
      
      /// this is a function that implements the solution steps. Before this 
      /// function is called, the matrices should already have been set.
      /// @param rhs right hand side of the system of equations
      /// @param sol vector in which solution will be stored
      virtual void solve(NumericVector<double>& rhs,
                         NumericVector<double>& sol);
      
      
      /// @returns the Petsc KSP object
      KSP getKSP();
      
      /// @returns the Petsc PC object
      PC getPC();
    
    protected:
      
      // creates the solver context, or recreates if it already exists
      void createSolverContext();
      
      DeclException1(ExcInvalidKSPType, std::string, 
                     << "Invalid KSP type : " << arg1 << " specified");
     
      /// the krylov solver context
      KSP ksp;
      
      /// the type of KSP solver being used
      KSPType ksptype;
      
      /// the PC context that will handle all the factoring of the matrices
      PC pc;
      
      /// type of the PC to be used
      PCType pctype;
      
      PetscDrawLG petsc_draw;
      
    };
}





#endif //__fesystem_petsc_linear_solver_h__



