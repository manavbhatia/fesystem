// $Id: LinearAnalysisDriver.C,v 1.25.6.6 2008/08/21 00:50:52 manav Exp $

// C++ includes
#include <fstream>
#include <sstream>

// FESystem includes
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/LinearSolverInfo.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/DisciplineInfo.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"
#include "Solutions/SolutionBase.h"

// libMesh includes
#include "base/dof_map.h"
#include "geom/node.h"
#include "numerics/sparse_matrix.h"

Driver::LinearAnalysisDriver::LinearAnalysisDriver(const unsigned int ID,
                                                   FESystem::FESystemController& fesys_controller):
Driver::AnalysisDriver(ID, fesys_controller, LINEAR_ANALYSIS_DRIVER::num())
{  
}





Driver::LinearAnalysisDriver::~LinearAnalysisDriver()
{
	
}




void
Driver::LinearAnalysisDriver::addMatricesAndVectors()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  // add a system matrix, a solution vector and a rhs
  it->second->initMatrix(this->addMatrix("SystemMatrix"));
  it->second->initMatrix(this->addMatrix("SystemMatrixSensitivity"));
  it->second->initVector(this->addVector("Solution"));
  it->second->initVector(this->addVector("RHS"));
}



void
Driver::LinearAnalysisDriver::solveCurrentLoadCase()
{    
  
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  NumericVector<double>& solution = this->getVector("Solution");
  NumericVector<double>& rhs = this->getVector("RHS");
  SparseMatrix<double>& matrix = this->getMatrix("SystemMatrix");
  Solver::LinearSolver* linear_solver = dynamic_cast<Solver::LinearSolver*>(this->solver);
  
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
    
  // prepare the load vector for the analysis and then init the solver and then finally 
  // solve the system
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      matrix_map[Driver::SYSTEM_MATRIX::num()] = &matrix;
      vector_map[Driver::FORCE_VECTOR::num()] = &rhs;
      
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      
      matrix.close();
      rhs.close();

      this->applyBoundaryConditionsToVectorForCurrentAnalysis(rhs, true, it->second->getDisciplineEnumID());
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis(matrix, 1.0, it->second->getDisciplineEnumID());
      
      matrix.close();
      rhs.close();
    }
      break;
      
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      // for a sensitivity analysis, get the dK/dx and myltiply it with the solution vector 
      // then get dF/dx and calculate the RHS 
      
      
      // multiply this matrix with the solution vector 
      SparseMatrix<double>& matrix_sens = this->getMatrix("SystemMatrixSensitivity"); 
      SparseMatrix<double>& matrix = this->getMatrix("SystemMatrix"); 
      
      matrix_map[Driver::SYSTEM_MATRIX::num()] = &matrix;
      matrix_map[Driver::SYSTEM_MATRIX_SENSITIVITY::num()] = &matrix_sens;
      vector_map[Driver::FORCE_VECTOR_SENSITIVITY::num()] = &rhs;
      
      it->second->fillQuantity(real_map, matrix_map, vector_map); 
      
      
      // create a numeric vector and get the solution from the database
      std::auto_ptr<NumericVector<double> > 
      sol(NumericVector<double>::build().release()), 
      prod_vec(NumericVector<double>::build().release());
      
      // create a numeric vector and get the solution from the database
      FESystemDatabase::TimeIndependentDataInfo data_info;
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setLoadCase(this->current_load_case);
      
      this->fesystem_controller.global_data_storage->fillVector(data_info, *(sol.get()));
      
      // init second vector to the appropriate size
      prod_vec->init(sol->size());
      
      // perform the multiplication
      matrix_sens.multiply_vector(*(sol.get()), *(prod_vec.get()));
      
      // finally, perform the subtraction to get the RHS
      rhs.add(-1.0, *(prod_vec.get()));
      matrix.close();
      rhs.close();

      this->applyBoundaryConditionsToVectorForCurrentAnalysis
      (rhs, false, it->second->getDisciplineEnumID());
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis(matrix, 1.0, it->second->getDisciplineEnumID());

      matrix.close();
      rhs.close();

    }
      break;
			
    default:
      abort();
      break;
			
  }
  
  
  // now solve
  solution.zero();
  solution.close();
  linear_solver->setSystemMatrix(matrix);
  linear_solver->setPreconditionerMatrix(matrix);
  linear_solver->solve(rhs, solution);
  
  FESystemDatabase::TimeIndependentDataInfo data_info;
  data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
  data_info.setName("Solution");
  data_info.setLoadCase(this->current_load_case);
  if (this->getCurrentAnalysisKind() == Driver::SENSITIVITY_ANALYSIS::num())
    data_info.setDVID(this->getCurrentDesignParameter().getID());
  
  
  // save the solution vector to the database
  this->fesystem_controller.global_data_storage->storeVector(data_info,solution);
  this->writeSolutionVectorToLog(data_info, solution, it->second->getDisciplineEnumID());
}




void
Driver::LinearAnalysisDriver::initialize()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());

  // add the matrices and vectors
  this->addMatricesAndVectors();
}

