// $Id: StructuralAnalysis.h,v 1.11.6.7 2008-08-21 00:37:18 manav Exp $

#ifndef __fesystem_structural_analysis_h__
#define __fesystem_structural_analysis_h__

// C++ includes 


// FESystem includes 
#include "Discipline/AnalysisDisciplineBase.h"
#include "Utilities/NameEnumHandler.h"
#include "Loads/load.h"


// libMesh includes
#include "base/libmesh_common.h"

// Forward Decleration
namespace FESystem
{
  class FESystemController;
}


namespace FESystemElem
{
  class FESystemElemBase;
}


#ifndef STRUCTURAL_DISCIPLINE_ENUM_ID
#define STRUCTURAL_DISCIPLINE_ENUM_ID 1
#else
#error
#endif

#ifndef STRUCTURAL_DISCIPLINE_ENUM_NAME
#define STRUCTURAL_DISCIPLINE_ENUM_NAME "STRUCTURAL_DISCIPLINE"
#else
#error
#endif



#ifndef NODAL_FORCE_ENUM_ID
#define NODAL_FORCE_ENUM_ID 7
#else
#error
#endif

#ifndef NODAL_FORCE_ENUM_NAME
#define NODAL_FORCE_ENUM_NAME "NODAL_FORCE"
#else
#error
#endif



#ifndef NODAL_TEMPERATURE_ENUM_ID
#define NODAL_TEMPERATURE_ENUM_ID 8
#else
#error
#endif

#ifndef NODAL_TEMPERATURE_ENUM_NAME
#define NODAL_TEMPERATURE_ENUM_NAME "NODAL_TEMPERATURE"
#else
#error
#endif



#ifndef SURFACE_PRESSURE_ENUM_ID
#define SURFACE_PRESSURE_ENUM_ID 9
#else
#error
#endif

#ifndef SURFACE_PRESSURE_ENUM_NAME
#define SURFACE_PRESSURE_ENUM_NAME "SURFACE_PRESSURE"
#else
#error
#endif


#ifndef DISPLACEMENT_BOUNDARY_CONDITION_ENUM_ID
#define DISPLACEMENT_BOUNDARY_CONDITION_ENUM_ID 10
#else
#error
#endif

#ifndef DISPLACEMENT_BOUNDARY_CONDITION_ENUM_NAME
#define DISPLACEMENT_BOUNDARY_CONDITION_ENUM_NAME "DISPLACEMENT_BOUNDARY_CONDITION"
#else
#error
#endif



#ifndef SELF_WEIGHT_ENUM_ID
#define SELF_WEIGHT_ENUM_ID 11
#else
#error
#endif

#ifndef SELF_WEIGHT_ENUM_NAME
#define SELF_WEIGHT_ENUM_NAME "SELF_WEIGHT"
#else
#error
#endif



DeclareEnumName(NODAL_FORCE, LoadNameEnum,
		NODAL_FORCE_ENUM_ID,
		NODAL_FORCE_ENUM_NAME);

DeclareEnumName(NODAL_TEMPERATURE, LoadNameEnum,
		NODAL_TEMPERATURE_ENUM_ID,
		NODAL_TEMPERATURE_ENUM_NAME);
  
DeclareEnumName(SURFACE_PRESSURE, LoadNameEnum,
		SURFACE_PRESSURE_ENUM_ID,
		SURFACE_PRESSURE_ENUM_NAME);


DeclareEnumName(DISPLACEMENT_BOUNDARY_CONDITION, LoadNameEnum,
                DISPLACEMENT_BOUNDARY_CONDITION_ENUM_ID,
                DISPLACEMENT_BOUNDARY_CONDITION_ENUM_NAME);

DeclareEnumName(SELF_WEIGHT, LoadNameEnum,
                SELF_WEIGHT_ENUM_ID,
                SELF_WEIGHT_ENUM_NAME);


/// This class provides methods to calculate the matrices and vectors related to the 
/// structural discipline. 
namespace Discipline
{
  
  class StructuralDisciplineInfo;

  
  DeclareEnumName(STRUCTURAL_DISCIPLINE, Discipline::AnalysisDisciplineEnum,
                  STRUCTURAL_DISCIPLINE_ENUM_ID,
                  STRUCTURAL_DISCIPLINE_ENUM_NAME);
  

  

  
  
  class StructuralAnalysis: public Discipline::AnalysisDisciplineBase
    {
public:
      
      /// constructor. It takes two arguements,
      /// the arguement is a pointer to the AnalysisDriver that owns this object
      StructuralAnalysis(FESystem::FESystemController& controller,
                         const Discipline::StructuralDisciplineInfo& info);
      
      /// destructor
      ~StructuralAnalysis();
      
      /// returns if same system matrix can be used for all load cases
      virtual bool sameSystemMatrixForLoadCases();
      
      /// @returns the eigenproblem kind for the discipline.
      virtual unsigned int getEigenProblemKindEnumID() ;
      
      /// @returns the order of the transient nature of the problem
      virtual unsigned int getTransientSystemOrder();

      /// returns the enum ID of the boundary condition for this disicipline
      virtual std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
        getBoundaryConditionLoadInfo() const;
      
      /// @returns true if the specified matrix has exists for the discipline
      virtual bool disciplineHasMatrix(const unsigned int qty_enum_ID) const;

      /// @returns the data info object for the specified quantity, for the current analysis
      /// in progress
      virtual std::auto_ptr<FESystemDatabase::DataInfoBase> 
        getDataInfoForQty(const unsigned int qty_enum_ID);

      /// This method updates the nodal locations of the mesh from a displacement vector. 
      /// The vector provides the values of the displacements relative to the current
      /// mesh geometry, so that the displacements are simply added to the current 
      /// nodal locations to get the new mesh. 
      void setMeshNodalIncrementsFromDisplacementSolution(const NumericVector<double>& vec);
      
protected:
        
        /// function to add analysis variables
        void addAnalysisVariables();
      
            
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
    };
}







#endif // __fesystem_structural_analysis_h__