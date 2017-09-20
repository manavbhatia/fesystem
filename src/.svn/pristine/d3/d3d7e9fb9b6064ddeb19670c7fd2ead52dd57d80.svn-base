/*
 *  PistonTheory.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */


#ifndef __fesystem_piston_theory_h__
#define __fesystem_piston_theory_h__

// FESystem includes
#include "Discipline/AerodynamicDisciplineBase.h"
#include "Loads/load.h"


// forward declerations
namespace FESystem{
  class FESystemController;
}


#ifndef PISTON_THEORY_ENUM_ID
#define PISTON_THEORY_ENUM_ID 4
#else
#error
#endif

#ifndef PISTON_THEORY_ENUM_NAME
#define PISTON_THEORY_ENUM_NAME "PISTON_THEORY"
#else
#error
#endif


#ifndef PISTON_THEORY_SURFACE_ENUM_ID
#define PISTON_THEORY_SURFACE_ENUM_ID 12
#else
#error
#endif

#ifndef PISTON_THEORY_SURFACE_ENUM_NAME
#define PISTON_THEORY_SURFACE_ENUM_NAME "PISTON_THEORY_SURFACE"
#else
#error
#endif


DeclareEnumName(PISTON_THEORY_SURFACE, LoadNameEnum,
                PISTON_THEORY_SURFACE_ENUM_ID,
                PISTON_THEORY_SURFACE_ENUM_NAME);




namespace Discipline {
  
  class PistonTheoryInfo;
  
  DeclareEnumName(PISTON_THEORY, Discipline::AnalysisDisciplineEnum,
                  PISTON_THEORY_ENUM_ID,
                  PISTON_THEORY_ENUM_NAME);
  
  class PistonTheory: public AerodynamicDisciplineBase
    {
    public:
      /// constructor
      PistonTheory(FESystem::FESystemController& controller,
                   const Discipline::PistonTheoryInfo& info);
      
      /// destructor
      virtual ~PistonTheory();
      
      /// method clears the data structures of this object
      virtual void clear();
      
      /// returns if same system matrix can be used for all load cases
      virtual inline bool sameSystemMatrixForLoadCases();
      
      /// @returns the eigenproblem kind for the discipline.
      virtual unsigned int getEigenProblemKindEnumID() ;
      
      /// @returns the order of the transient nature of the problem
      virtual inline unsigned int getTransientSystemOrder();
      
      /// returns the enum ID of the boundary condition for this disicipline
      virtual std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
      getBoundaryConditionLoadInfo() const;
      
      /// @returns true if the specified matrix has exists for the discipline
      virtual bool disciplineHasMatrix(const unsigned int qty_enum_ID) const;
      
      /// @returns the data info object for the specified quantity, for the current analysis
      /// in progress
      virtual std::auto_ptr<FESystemDatabase::DataInfoBase> 
      getDataInfoForQty(const unsigned int qty_enum_ID);
      
      /// attach structural discipline
      void attachStructuralDiscipline(Discipline::AnalysisDisciplineBase& discipline);
      
    protected:
      
      /// function to add analysis variables
      virtual void addAnalysisVariables();
      
      
      /// function will fill the input vector with the element quantities required for the calculation
      /// of the specified global quantity. This is an abstract function and will be implemented in the
      /// derived class
      /// @param global quantitiy
      /// @param element pointer
      /// @param vector in which the quantities will be added
      void getElemQty(FESystemElem::FESystemElemBase* elem,
                      std::map<unsigned int, double>& real_elem_data,
                      std::map<unsigned int, DenseMatrix<double> >& mat_elem_data,
                      std::map<unsigned int, DenseVector<double> >& vec_elem_data,
                      const bool sensitivity,
                      const unsigned int DV_ID);
      
      
      /// function to calculate the conductance matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateTransientMatrix(const unsigned int order, bool = false);
      
      /// function to calculate the conductance matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateK(bool = false);
      
      /// function to calculate the Jacobian matrix. First arguement is the matrix to store 
      ///  the result in
      virtual bool calculateJac();
      
      /// function to calculate the A matrix for eigen problem analysis
      virtual bool calculateEigenProblemAMatrix(bool sensitivity = false);
      
      /// function to calculate the B matrix for eigen problem analysis
      virtual bool calculateEigenProblemBMatrix(bool sensitivity = false);
      
      /// function to calculate the force vector. First arguement is the vector to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateF(bool = false);
      
      /// this method will perform some discipline specific operations once the 
      /// global quantities have been calculated
      virtual void performAdditionalOperationsOnGlobalQuantities
      (std::map<unsigned int, SparseMatrix<double>* >& matrix_map,
       std::map<unsigned int, NumericVector<double>* >& vector_map);
      
      /// this method fills the vector with the solution vector necessary for filling the quantity
      /// for the current analysis. It checks the current status of the discipline, and 
      /// apropriately loads the solution vector
      virtual const NumericVector<double>& getSolutionVector(const unsigned int transient_order,
                                                             const bool sensitivity);
      
      
      /// the structural discipline for which this analysis is being performed
      Discipline::AnalysisDisciplineBase* structural_discipline;
      
    };

  bool PistonTheory::sameSystemMatrixForLoadCases()
  {
    return false;
  }


  unsigned int PistonTheory::getTransientSystemOrder()
  {
    return 1;
  }
}


#endif // __fesystem_piston_theory_h__

