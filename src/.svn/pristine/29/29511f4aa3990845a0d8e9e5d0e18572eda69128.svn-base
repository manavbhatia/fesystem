// $Id: AnalysisDriver.h,v 1.24.6.6 2008/06/03 05:19:20 manav Exp $

#ifndef __fesystem_analysis_driver_h__
#define __fesystem_analysis_driver_h__

// C++ includes
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>


// FESystem includes
#include "Utilities/NameEnumHandler.h"
#include "FESystem/FESystemController.h"
#include "Database/GlobalDataStorage.h"
#include "Utilities/AutoptrVector.h"


// libMesh includes
#include "fe/fe_type.h"


// Forward decleraions
namespace FESystem
{
  class FESystemController;
}

namespace DesignData
{
  class DesignParameter;
}

namespace Solver
{
  class FESystemSolverBase;
}

namespace Discipline
{
  class DisciplineInfo;
}

namespace Discipline
{
  class AnalysisDisciplineBase;
}


namespace Solution
{
  class SolutionBase;
}

#ifndef ANALYSIS_ENUM_ID
#define ANALYSIS_ENUM_ID 1
#else
#error
#endif

#ifndef ANALYSIS_ENUM_NAME
#define ANALYSIS_ENUM_NAME "ANALYSIS"
#else
#error
#endif

#ifndef SENSITIVITY_ANALYSIS_ENUM_ID
#define SENSITIVITY_ANALYSIS_ENUM_ID 2
#else
#error
#endif

#ifndef SENSITIVITY_ANALYSIS_ENUM_NAME
#define SENSITIVITY_ANALYSIS_ENUM_NAME "SENSITIVITY_ANALYSIS"
#else
#error
#endif


#ifndef POST_PROCESS_ENUM_ID
#define POST_PROCESS_ENUM_ID 3
#else
#error
#endif

#ifndef POST_PROCESS_ENUM_NAME
#define POST_PROCESS_ENUM_NAME "POST_PROCESS"
#else
#error
#endif



#ifndef SYSTEM_MATRIX_ENUM_ID
#define SYSTEM_MATRIX_ENUM_ID 1
#else
#error
#endif

#ifndef SYSTEM_MATRIX_ENUM_NAME
#define SYSTEM_MATRIX_ENUM_NAME "SYSTEM_MATRIX"
#else
#error
#endif

#ifndef FORCE_VECTOR_ENUM_ID
#define FORCE_VECTOR_ENUM_ID 2
#else
#error
#endif

#ifndef FORCE_VECTOR_ENUM_NAME
#define FORCE_VECTOR_ENUM_NAME "FORCE_VECTOR"
#else
#error
#endif

#ifndef JACOBIAN_MATRIX_ENUM_ID
#define JACOBIAN_MATRIX_ENUM_ID 3
#else
#error
#endif

#ifndef JACOBIAN_MATRIX_ENUM_NAME
#define JACOBIAN_MATRIX_ENUM_NAME "JACOBIAN_MATRIX"
#else
#error
#endif

//#ifndef RESIDUAL_VECTOR_ENUM_ID
//#define RESIDUAL_VECTOR_ENUM_ID 4
//#else
//#error
//#endif
//
//#ifndef RESIDUAL_VECTOR_ENUM_NAME
//#define RESIDUAL_VECTOR_ENUM_NAME "RESIDUAL_VECTOR"
//#else
//#error
//#endif

#ifndef SYSTEM_MATRIX_SENSITIVITY_ENUM_ID
#define SYSTEM_MATRIX_SENSITIVITY_ENUM_ID 5
#else
#error
#endif

#ifndef SYSTEM_MATRIX_SENSITIVITY_ENUM_NAME
#define SYSTEM_MATRIX_SENSITIVITY_ENUM_NAME "SYSTEM_MATRIX_SENSITIVITY"
#else
#error
#endif

#ifndef FORCE_VECTOR_SENSITIVITY_ENUM_ID
#define FORCE_VECTOR_SENSITIVITY_ENUM_ID 6
#else
#error
#endif

#ifndef FORCE_VECTOR_SENSITIVITY_ENUM_NAME
#define FORCE_VECTOR_SENSITIVITY_ENUM_NAME "FORCE_VECTOR_SENSITIVITY"
#else
#error
#endif


#ifndef EIGENPROBLEM_A_MATRIX_ENUM_ID
#define EIGENPROBLEM_A_MATRIX_ENUM_ID 7
#else
#error
#endif


#ifndef EIGENPROBLEM_A_MATRIX_ENUM_NAME
#define EIGENPROBLEM_A_MATRIX_ENUM_NAME "EIGENPROBLEM_A_MATRIX"
#else
#error
#endif



#ifndef EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID
#define EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID 8
#else
#error
#endif

#ifndef EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_NAME
#define EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_NAME "EIGENPROBLEM_A_MATRIX_SENSITIVITY"
#else
#error
#endif



#ifndef EIGENPROBLEM_B_MATRIX_ENUM_ID
#define EIGENPROBLEM_B_MATRIX_ENUM_ID 9
#else
#error
#endif


#ifndef EIGENPROBLEM_B_MATRIX_ENUM_NAME
#define EIGENPROBLEM_B_MATRIX_ENUM_NAME "EIGENPROBLEM_B_MATRIX"
#else
#error
#endif



#ifndef EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID
#define EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID 10
#else
#error
#endif

#ifndef EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_NAME
#define EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_NAME "EIGENPROBLEM_B_MATRIX_SENSITIVITY"
#else
#error
#endif


#ifndef MODEL_MASS_ENUM_ID
#define MODEL_MASS_ENUM_ID 15
#else
#error
#endif

#ifndef MODEL_MASS_ENUM_NAME
#define MODEL_MASS_ENUM_NAME "MODEL_MASS"
#else
#error
#endif


#ifndef MODEL_MASS_SENSITIVITY_ENUM_ID
#define MODEL_MASS_SENSITIVITY_ENUM_ID 16
#else
#error
#endif

#ifndef MODEL_MASS_SENSITIVITY_ENUM_NAME
#define MODEL_MASS_SENSITIVITY_ENUM_NAME "MODEL_MASS_SENSITIVITY"
#else
#error
#endif


namespace Driver
{
  
  DeclareEnumClass(AnalysisDriverTypeEnum);
  
  DeclareEnumClass(AnalysisTypeEnum);
  
  DeclareEnumName(ANALYSIS, AnalysisTypeEnum,
                  ANALYSIS_ENUM_ID, ANALYSIS_ENUM_NAME);
  
  DeclareEnumName(SENSITIVITY_ANALYSIS, AnalysisTypeEnum,
                  SENSITIVITY_ANALYSIS_ENUM_ID, 
                  SENSITIVITY_ANALYSIS_ENUM_NAME);
  
  DeclareEnumName(POST_PROCESS, AnalysisTypeEnum,
                  POST_PROCESS_ENUM_ID, 
                  POST_PROCESS_ENUM_NAME);
  
  DeclareEnumClass(AnalysisDriverQtyEnum);
  
  DeclareEnumName(SYSTEM_MATRIX, AnalysisDriverQtyEnum,
                  SYSTEM_MATRIX_ENUM_ID, 
                  SYSTEM_MATRIX_ENUM_NAME);
  
  DeclareEnumName(FORCE_VECTOR, AnalysisDriverQtyEnum,
                  FORCE_VECTOR_ENUM_ID, 
                  FORCE_VECTOR_ENUM_NAME);
  
  //  DeclareEnumName(RESIDUAL_VECTOR, AnalysisDriverQtyEnum,
  //                  RESIDUAL_VECTOR_ENUM_ID, 
  //                  RESIDUAL_VECTOR_ENUM_NAME);
  
  DeclareEnumName(JACOBIAN_MATRIX, AnalysisDriverQtyEnum,
                  JACOBIAN_MATRIX_ENUM_ID, 
                  JACOBIAN_MATRIX_ENUM_NAME);
  
  DeclareEnumName(SYSTEM_MATRIX_SENSITIVITY, AnalysisDriverQtyEnum,
                  SYSTEM_MATRIX_SENSITIVITY_ENUM_ID, 
                  SYSTEM_MATRIX_SENSITIVITY_ENUM_NAME);
  
  DeclareEnumName(FORCE_VECTOR_SENSITIVITY, AnalysisDriverQtyEnum,
                  FORCE_VECTOR_SENSITIVITY_ENUM_ID, 
                  FORCE_VECTOR_SENSITIVITY_ENUM_NAME);
  
  DeclareEnumName(EIGENPROBLEM_A_MATRIX, AnalysisDriverQtyEnum,
                  EIGENPROBLEM_A_MATRIX_ENUM_ID, 
                  EIGENPROBLEM_A_MATRIX_ENUM_NAME);
  
  DeclareEnumName(EIGENPROBLEM_A_MATRIX_SENSITIVITY, AnalysisDriverQtyEnum,
                  EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_ID, 
                  EIGENPROBLEM_A_MATRIX_SENSITIVITY_ENUM_NAME);
  
  DeclareEnumName(EIGENPROBLEM_B_MATRIX, AnalysisDriverQtyEnum,
                  EIGENPROBLEM_B_MATRIX_ENUM_ID, 
                  EIGENPROBLEM_B_MATRIX_ENUM_NAME);
  
  DeclareEnumName(EIGENPROBLEM_B_MATRIX_SENSITIVITY, AnalysisDriverQtyEnum,
                  EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_ID, 
                  EIGENPROBLEM_B_MATRIX_SENSITIVITY_ENUM_NAME);

  DeclareEnumName(MODEL_MASS, AnalysisDriverQtyEnum,
                  MODEL_MASS_ENUM_ID, 
                  MODEL_MASS_ENUM_NAME);

  DeclareEnumName(MODEL_MASS_SENSITIVITY, AnalysisDriverQtyEnum,
                  MODEL_MASS_SENSITIVITY_ENUM_ID, 
                  MODEL_MASS_SENSITIVITY_ENUM_NAME);
  
  
  class AnalysisDriver
    {
    public:
      
      /// consturctor
      /// @param controller fesystem controller that owns this object
      /// @param discipline discipline for this driver
      AnalysisDriver(const unsigned int ID,
                     FESystem::FESystemController& controller, 
                     const unsigned int driver_type_id);
      
      
      /// destructor
      virtual ~AnalysisDriver();
      
      /// @returns the FESystemController class that owns this object
      inline FESystem::FESystemController& getFESystemController();
      
      /// @returns the ID of this driver
      inline unsigned int getID() const;
      
      /// @returns the analysis driver type enum ID
      inline unsigned int getAnalysisDriverTypeEnumID();
      
      /// @returns the discipline
      inline Discipline::AnalysisDisciplineBase& getAnalysisDiscipline();
      
      /// @returns the solver
      inline Solver::FESystemSolverBase& getSolver();
      
      /// ataches the solver to this driver for analysis solutions
      void attachSolver(Solver::FESystemSolverBase* solver);
      
      /// attaches  a solution to this driver
      inline void attachSolution(Solution::SolutionBase* sol);
      
      /// Adds the vector \p vec_name. All vectors are inititialized to zero.
      NumericVector<double> & addVector (const std::string& vec_name);
      
      
      /// @returns a writeable reference to the vector
      /// named \p vec_name.
      inline NumericVector<double> & getVector (const std::string& vec_name);
      
      
      /// @returns \p true if this \p AnalysisDriver has a vector associated with the
      /// given name, \p false otherwise.
      inline bool haveVector (const std::string& vec_name) const;
      
      
      /// Adds the additional matrix \p mat_name to this analysis. 
      SparseMatrix<double> & addMatrix (const std::string& mat_name);
      
      
      /// @returns a const reference to this system's matrix
      /// named \p mat_name. 
      inline SparseMatrix<double> & getMatrix (const std::string& mat_name) const;
      
      
      /// @returns \p true if this \p AnalysisDriver has a matrix associated with the
      /// given name, \p false otherwise.
      inline bool haveMatrix (const std::string& mat_name) const;
      
      
      /// returns the current load case number
      inline unsigned int getCurrentLoadCase() const;
      
      
      /// returns the current design variable
      inline const DesignData::DesignParameter& getCurrentDesignParameter() const;
      
      
      /// @returns the current solution vector for the system
      NumericVector<double>& getCurrentSolution(const unsigned int transient_order);
      
      /// @returns the kind of analysis being performed
      inline unsigned int getCurrentAnalysisKind() const;
      
      /// attaches analysis discipline for solution
      void attachAnalysisDiscipline(Discipline::AnalysisDisciplineBase* discipline);
      
      /// clears the data structure
      virtual void clear();
      
      /// this will loop over all load cases, and perform basic analysis for each
      void solveForLoadCases(const std::vector<unsigned int>& load_cases);
      
      /// this will loop over all load cases, and perform sensitivity analysis for 
      /// each
      void solveForLoadCaseSensitivity
      (const std::vector<DesignData::DesignParameter*>& design_params,
       const std::vector<unsigned int>& load_cases);
      
      /// this method performs the post process operations 
      virtual void postProcess(const std::vector<unsigned int>& load_cases, unsigned int disciplin_enum_ID);
      
    protected:
      
      
      /// this method initializes itself before the solutions are started for the driver
      virtual void initialize() = 0;
      
      /// writes a vector to the output stream
      void writeSolutionVectorToLog( FESystemDatabase::GenericDataInfoBase& data_info,
                                    const NumericVector<double>& sol, unsigned int disciplin_enum_ID);

      
      /// writes a vector to the output stream
      void writeComplexSolutionVectorToLog( FESystemDatabase::GenericDataInfoBase& data_info,
                                           const NumericVector<double>& real_sol, 
                                           const NumericVector<double>& img_sol, 
                                           unsigned int disciplin_enum_ID);
      
      
      /// apply dirichlet boundary condition to vector. The second arguement specifies if 
      /// the value of the BC should be applied to the DOF or if it should be zeroed. If true, 
      /// the value will be applied.
      void applyBoundaryConditionsToVectorForCurrentAnalysis(NumericVector<double>& vec, 
                                                             const bool if_apply_bc_vals,
                                                             unsigned int disciplin_enum_ID);
      
      
      /// this applies the dirichlet boundary condition the matrix. The matrix columns and rows 
      /// corresponding to the dofs with dirichlet BCs are zeroed, and the value specified 
      /// in the second arguement is set at the diagonal of all those dofs.
      void applyBoundaryConditionsToMatrixForCurrentAnalysis(SparseMatrix<double>& mat, 
                                                             double val,
                                                             unsigned int disciplin_enum_ID);
      
      
      /// localizes the solution vector. The first arguement is the vector from which
      /// the solution will be localized, and the second arguement will contain the 
      /// localized solution
      void localizeSolution(const DofMap& dof_map,
                            NumericVector<double>& source_vec, 
                            NumericVector<double>& localized_vec);
      
      /// this method solves the system for the current load case and analysis kind
      virtual void solveCurrentLoadCase() = 0;
      
      /// adds matrices and vectors to this analysis driver
      virtual void addMatricesAndVectors() = 0;
      
      /// this exception is thrown if an analysis discipline is attached without 
      /// clearing the object. 
      DeclException0(ExcInitBeforeClear);
      
      /// this exception is thrown when an invalid solver is specified in the input
      DeclException0(ExcInvalidSolver);
      
      DeclException2(ExcObjectNameDoesNotExist, std::string&, std::string&,
                     << arg1 << " Name: " << arg2 << "does not exist." );
      
      
      /// a number to identify this driver
      const unsigned int driver_ID;
      
      /// refernce to the FESystemController that owns this object
      FESystem::FESystemController& fesystem_controller;
      
      /// analysis discipline for this driver
      const unsigned int analysis_driver_type_enum_ID;
      
      /// solver that the analysis driver will use
      Solver::FESystemSolverBase* solver;
      
      /// solution for this analysis driver for the current analysis
      Solution::SolutionBase* solution_base;
      
      /// load case for the analysis
      unsigned int current_load_case;	
      
      /// design variable variable for which sensitivity analysis is being
      /// performed
      DesignData::DesignParameter *current_DV;
      
      /// this stores the current iterate for different (nonlinear/transient) analysis. 
      /// for transient solutions, the time derivative of the solution is stored in the 
      /// location in the vector equal to its derivative order.
      FESystemUtility::AutoPtrVector<NumericVector<double> > current_solution; 
      
      /// analysis being performed 
      unsigned int current_analysis_kind_enum_ID;
      
      /// Some systems need an arbitrary number of vectors.
      /// This map allows names to be associated with arbitrary
      /// vectors.  All the vectors in this map will be distributed
      /// in the same way as the solution vector.
      std::map<std::string, NumericVector<double>* > vectors;
      
      /// Some systems need an arbitrary number of matrices.
      std::map<std::string, SparseMatrix<double>* > matrices;

      /// map of analysis disciplines used in this driver
      std::map<unsigned int, Discipline::AnalysisDisciplineBase*> analysis_discipline_map;

    };
}



inline
void 
Driver::AnalysisDriver::attachSolution(Solution::SolutionBase* sol)
{
  Assert(sol != NULL, ExcEmptyObject());
  
  this->solution_base = sol;
}


inline 
FESystem::FESystemController&
Driver::AnalysisDriver::getFESystemController()
{
  return this->fesystem_controller;
}



inline 
unsigned int
Driver::AnalysisDriver::getID() const
{
  return this->driver_ID;
}



inline
unsigned int 
Driver::AnalysisDriver::getAnalysisDriverTypeEnumID()
{
  return this->analysis_driver_type_enum_ID;
}




//inline 
//Discipline::AnalysisDisciplineBase& 
//Driver::AnalysisDriver::getAnalysisDiscipline()
//{
//  Assert(this->analysis_discipline != NULL,
//         ExcEmptyObject());
//  return *(this->analysis_discipline);
//}




inline 
Solver::FESystemSolverBase& 
Driver::AnalysisDriver::getSolver()
{
  Assert(this->solver != NULL,
         ExcEmptyObject());
  return *(this->solver);
}


inline 
NumericVector<double> &
Driver::AnalysisDriver::getVector (const std::string& vec_name)
{
  // Make sure the vector exists
  std::map<std::string, NumericVector<double>* >::iterator pos = this->vectors.find(vec_name);
	
  Assert(pos != this->vectors.end(),
         Driver::AnalysisDriver::ExcObjectNameDoesNotExist("Vector", vec_name));
  
  return *(pos->second);
}





inline 
bool
Driver::AnalysisDriver::haveVector (const std::string& vec_name) const
{
  return (this->vectors.count(vec_name));
}




inline 
SparseMatrix<double> & 
Driver::AnalysisDriver::getMatrix(const std::string& mat_name) const
{
  // Make sure the matrix exists
  std::map<std::string, SparseMatrix<double>* >::const_iterator pos = this->matrices.find (mat_name);
	
  Assert(pos != this->matrices.end(),
         Driver::AnalysisDriver::ExcObjectNameDoesNotExist("Matrix", mat_name));
	
  return *(pos->second);
}




inline 
bool
Driver::AnalysisDriver::haveMatrix (const std::string& mat_name) const
{
  return (this->matrices.count(mat_name));
}




inline 
unsigned int 
Driver::AnalysisDriver::getCurrentLoadCase() const
{
  return this->current_load_case;
}




inline 
const DesignData::DesignParameter& 
Driver::AnalysisDriver::getCurrentDesignParameter() const
{
  return *(this->current_DV);
}




inline
unsigned int
Driver::AnalysisDriver::getCurrentAnalysisKind() const
{
  return this->current_analysis_kind_enum_ID;
}


#endif // __fesystem_analysis_driver_h__
