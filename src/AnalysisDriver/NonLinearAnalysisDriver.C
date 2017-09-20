// $Id: NonLinearAnalysisDriver.C,v 1.23.4.8 2008/08/21 20:43:10 manav Exp $

// C++ includes
#include <fstream>
#include <sstream>

// FESystem includes
#include "AnalysisDriver/NonLinearAnalysisDriver.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/NonlinearSolver.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/DisciplineInfo.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"
#include "Solutions/SolutionBase.h"


// libMesh includes
#include "base/dof_map.h"
#include "geom/node.h"
#include "numerics/petsc_vector.h"

// PETSc includes
#include "petscvec.h"

Driver::NonLinearAnalysisDriver::NonLinearAnalysisDriver
(const unsigned int ID,
 FESystem::FESystemController& fesys_controller):
Driver::AnalysisDriver(ID, fesys_controller, NONLINEAR_ANALYSIS_DRIVER::num()),
iteration_number(0),
n_iterations(0)
{
}





Driver::NonLinearAnalysisDriver::~NonLinearAnalysisDriver()
{
	
}



void
Driver::NonLinearAnalysisDriver::clear()
{
	this->iteration_number = 0;
  this->n_iterations = 0;
  
  Driver::AnalysisDriver::clear();
}



void
Driver::NonLinearAnalysisDriver::initialize()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // now resize the current solution vector
  Assert(this->current_solution.size() == 0,ExcInternalError());
  this->current_solution.resize(1);
  for (unsigned int i=0; i < this->current_solution.size(); i++)
    {
      this->current_solution.reset(i, NumericVector<double>::build().release());
      this->current_solution[i]->init(it->second->getAnalysisDofMap().n_dofs());
    }
  
  // add the matrices and vectors
  this->addMatricesAndVectors();
}



void Driver::NonLinearAnalysisDriver::addMatricesAndVectors()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  // add a system matrix, a solution vector and a rhs
  it->second->initMatrix(this->addMatrix("JacobianMatrix"));
  it->second->initMatrix(this->addMatrix("ScratchMatrix"));
  it->second->initVector(this->addVector("Solution"));
  it->second->initVector(this->addVector("Residual"));
  it->second->initVector(this->addVector("ScratchVector"));
}





void Driver::NonLinearAnalysisDriver::solveCurrentLoadCase()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      this->setInitialGuessForCurrentLoadCase();
      
      //  set the matrices and vectors for the solver before the solution
      Solver::NonlinearSolver* nonlinear_solver = 
      dynamic_cast<Solver::NonlinearSolver*>(this->solver);
      
      nonlinear_solver->attachMatrixAndVector(this->getMatrix("JacobianMatrix"),
                                              this->getVector("Residual"));
      
      nonlinear_solver->solve(this->getVector("Solution"));
    }
      break;
			
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      // for a sensitivity analysis, get the dK/dx and myltiply it with the solution vector
      // then get dF/dx and calculate the RHS
      // then get the system matrix, apply BC and return
      
      std::map<unsigned int, double> real_map;
      std::map<unsigned int, SparseMatrix<double>*> matrix_map;
      std::map<unsigned int, NumericVector<double>*> vector_map;
      
      // multiply this matrix with the solution vector
      SparseMatrix<double>& system_matrix = this->getMatrix("JacobianMatrix");
      SparseMatrix<double>& scratch_matrix = this->getMatrix("ScratchMatrix");
      NumericVector<double>& scratch_vec = this->getVector("ScratchVector");
      NumericVector<double>& rhs_vec = this->getVector("Residual");
      NumericVector<double>& sol_vec = this->getVector("Solution");
      
      system_matrix.zero();
      scratch_matrix.zero();
      scratch_vec.zero();
      rhs_vec.zero();
      sol_vec.zero();
      
      matrix_map[Driver::JACOBIAN_MATRIX::num()] = &system_matrix;
      matrix_map[Driver::SYSTEM_MATRIX_SENSITIVITY::num()] = &scratch_matrix;
      vector_map[Driver::FORCE_VECTOR_SENSITIVITY::num()] = &scratch_vec;
      
      // create a numeric vector and get the solution from the database
      FESystemDatabase::TimeIndependentDataInfo data_info;
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setLoadCase(this->current_load_case);
      
      this->fesystem_controller.global_data_storage->fillVector(data_info, sol_vec);
      
      // now localize the current solution vector
      this->localizeSolution(it->second->getAnalysisDofMap(), sol_vec, 
                             *(this->current_solution[0]));
      
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      
      // perform the multiplication
      scratch_matrix.multiply_vector(sol_vec, rhs_vec);
      
      // finally, perform the subtraction to get the RHS
      rhs_vec.add(-1.0, scratch_vec);
      rhs_vec.scale(-1.0);

      this->applyBoundaryConditionsToVectorForCurrentAnalysis
      (rhs_vec, false, it->second->getDisciplineEnumID());
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis
      (system_matrix, 1.0, it->second->getDisciplineEnumID());

      system_matrix.close();
      rhs_vec.close();
      
      // finally, solve the system
      Solver::LinearSolver* linear_solver = dynamic_cast<Solver::LinearSolver*>(this->solver);
      
      sol_vec.zero();	
      linear_solver->setSystemMatrix(system_matrix);
      linear_solver->setPreconditionerMatrix(system_matrix);
      linear_solver->solve(rhs_vec, sol_vec);
    }
      break;
      
    default:
      abort();
      break;      
  }
  
  
  
  FESystemDatabase::TimeIndependentDataInfo data_info;
  data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
  data_info.setName("Solution");
  data_info.setLoadCase(this->current_load_case);
  if (this->getCurrentAnalysisKind() == Driver::SENSITIVITY_ANALYSIS::num())
    data_info.setDVID(this->getCurrentDesignParameter().getID());
  
  // save the solution vector to the database
  NumericVector<double>& sol = this->getVector("Solution");
  this->fesystem_controller.global_data_storage->storeVector(data_info, sol);
  this->writeSolutionVectorToLog(data_info, sol, it->second->getDisciplineEnumID());
}




void Driver::NonLinearAnalysisDriver::setInitialGuessForCurrentLoadCase()
{  
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // get reference to the current solution vector from the solver
  NumericVector<double>& vector = this->getVector("Solution");
	
  // set the value of the vector to the initialization value/vector
  // specified in the input file
  
  // first check the method that has been specified for initialization
  std::string method = 
  this->fesystem_controller.analysis_case->getStringParameter("SOL_VECTOR_INIT_METHOD");
  
  if (method == "VALUE")
    {
      double init_value = 
      this->fesystem_controller.analysis_case->getRealParameter("SOL_VECTOR_INIT_VALUE");
      
      Vec sol_vec = dynamic_cast<PetscVector<double>&>(vector).vec();
      VecSet(sol_vec, init_value);
    }
  else if (method == "VECTOR")
    {
      // this needs to be implemented. This requires support for reading arbitrary vectors
      // and matrices from the disk.
      Assert(false, ExcInternalError());
      //       std::string vec_name = 
      // 	this->fesystem_controller.analysis_case->getStringParameter("SOL_VECTOR_INIT_VECTOR_FILE");
      //       unsigned int load_case = this->getCurrentLoadCase();
      //       this->fesystem_controller.global_data_storage->fillVector(, vector);
    }
  
  this->applyBoundaryConditionsToVectorForCurrentAnalysis(vector, true, it->second->getDisciplineEnumID());
}





void Driver::NonLinearAnalysisDriver::getResidualAtVector(NumericVector<double>& sol,
                                                          NumericVector<double>& res)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // set the current solution and calculate the residual
  Assert(this->current_solution.size() > 0, ExcInternalError());
  // now localize the current solution vector
  this->localizeSolution(it->second->getAnalysisDofMap(), sol, 
                         *(this->current_solution[0]));
  
  // if a basic analysis is being performed, then ask for the basic quantities, apply BC and return
  SparseMatrix<double>& matrix = this->getMatrix("ScratchMatrix");
  NumericVector<double>& scratch_vec = this->getVector("ScratchVector"); 
  matrix.zero();
  matrix.close();
  scratch_vec.zero();
  scratch_vec.close();
  
  res.zero();
  res.close();
  
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  matrix_map[Driver::SYSTEM_MATRIX::num()] = &matrix;
  vector_map[Driver::FORCE_VECTOR::num()] = &scratch_vec;
  
  // fill the matrices and vectors with the required quantities
  it->second->fillQuantity(real_map, matrix_map, vector_map);
  
  // perform multiplication K * X - F
  matrix.close();
  scratch_vec.close();
  matrix.multiply_vector(sol, res);
  res.add(-1.0, scratch_vec);
  res.close();
  
  this->applyBoundaryConditionsToVectorForCurrentAnalysis(res, false, it->second->getDisciplineEnumID());
  res.close();
  
  this->current_solution[0]->zero();
  
  // return control back to the solver
}



void Driver::NonLinearAnalysisDriver::getJacobianAtVector(NumericVector<double>& sol,
                                                          SparseMatrix<double>& jac)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // set the current solution and calculate the jacobian
  Assert(this->current_solution.size() > 0, ExcInternalError());
  // now localize the current solution vector
  this->localizeSolution(it->second->getAnalysisDofMap(), sol, 
                         *(this->current_solution[0]));
  
  // fill the system matrix with the jacobian matrix
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  jac.zero();
  jac.close();
  matrix_map[Driver::JACOBIAN_MATRIX::num()] = &jac;
  
  it->second->fillQuantity(real_map, matrix_map, vector_map);
  jac.close();
  
  this->applyBoundaryConditionsToMatrixForCurrentAnalysis
  (jac, 1.0, it->second->getDisciplineEnumID());
  jac.close();
  
  this->current_solution[0]->zero();
  // return control back to the solver
}


unsigned int Driver::NonLinearAnalysisDriver::currentNonlinearIterationNumber() const
{
  return (dynamic_cast<Solver::NonlinearSolver*>(this->solver)->getCurrentIterationNumber());
}
