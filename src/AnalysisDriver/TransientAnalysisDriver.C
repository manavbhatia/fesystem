// $Id: TransientAnalysisDriver.C,v 1.1.2.10 2008/08/21 00:50:53 manav Exp $

// C++ includes
#include <fstream>
#include <sstream>

// FESystem includes
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Solvers/TransientSolver.h"
#include "Solvers/LinearSolver.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/DisciplineInfo.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"
#include "Solutions/SolutionBase.h"
#include "Numerics/PetscSeqVector.h"
#include "Utilities/Log.h"


// libMesh includes
#include "base/dof_map.h"
#include "geom/node.h"
#include "numerics/petsc_vector.h"

// PETSc includes
#include "petscvec.h"

Driver::TransientAnalysisDriver::TransientAnalysisDriver
(const unsigned int ID,
 FESystem::FESystemController& fesys_controller,
 const unsigned int driver_enum_ID):
Driver::AnalysisDriver(ID, fesys_controller, driver_enum_ID),
current_analysis_time(0.0)
{
  
}





Driver::TransientAnalysisDriver::~TransientAnalysisDriver()
{
	
}



void
Driver::TransientAnalysisDriver::clear()
{
  this->current_analysis_time = 0.0;
  this->time_values.clear();
  
  Driver::AnalysisDriver::clear();
}



double
Driver::TransientAnalysisDriver::getCurrentAnalysisTime() const
{
  return this->current_analysis_time;
}


double
Driver::TransientAnalysisDriver::getSimulatedTime() const
{
  return dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getSimulatedTime();
}



unsigned int 
Driver::TransientAnalysisDriver::getCurrentIterationNumber() const
{
  return dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getCurrentIterationNumber();
}


unsigned int 
Driver::TransientAnalysisDriver::getSimulatedIterationNumber() const
{
  return dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getSimulatedIterationNumber();
}



void
Driver::TransientAnalysisDriver::setNewIteration(const double time,
                                                 const unsigned int iter_num,
                                                 const std::vector<NumericVector<double>*>& sol_vec)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // make sure that the numbers make some sense
  Assert(time == this->getSimulatedTime() && iter_num == this->getSimulatedIterationNumber(), 
         ExcInternalError());
  
  // depending upon the kind of analysis, check the iteration time
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
      this->time_values.push_back(time);
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
      Assert(time == this->time_values[iter_num-1], ExcInternalError());
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  // make sure that the number of solutions in the vector is same as the order of the
  // system
  unsigned int order = it->second->getTransientSystemOrder();
  
  Assert(sol_vec.size() == (order+1),
         ExcInternalError());
  
  FESystemDatabase::TimeDependentDataInfo data_info;
  
  // now store the solution vector
  for (unsigned int i=0; i <= order; i++)
    {
      data_info.clear();
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(i);
      if (this->current_analysis_kind_enum_ID == SENSITIVITY_ANALYSIS::num())
        data_info.setDVID(this->getCurrentDesignParameter().getID());	
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo(iter_num, time);
      
      this->fesystem_controller.global_data_storage->storeVector(data_info, *(sol_vec[i]));
      this->writeSolutionVectorToLog(data_info, *(sol_vec[i]), it->second->getDisciplineEnumID());
    }
}






void Driver::TransientAnalysisDriver::setInitialConditionForCurrentLoadCase
(std::vector<NumericVector<double>*>& vectors)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  unsigned int order = it->second->getTransientSystemOrder();
  
  // make sure that the right number of vectors are given 
  Assert(order > 0 , ExcInternalError());
  Assert(vectors.size() == order,	ExcInternalError());
  
  // get references to the necessary data strucutres
  LoadDatabase& load_database = *(this->fesystem_controller.load_database.get());
  const MeshDS::FEMeshData& mesh_data = it->second->getAnalysisMeshData();
  
  // set the value of the vector to the initialization value/vector
  // specified in the input file
  
  //   // first check the method that has been specified for initialization
  //   std::string method = 
  //     this->fesystem_controller.analysis_case->getStringParameter("SOL_VECTOR_INIT_METHOD");
  
  //   if (method == "VALUE")
  //     {
  //       double init_value = 
  // 	this->fesystem_controller.analysis_case->getRealParameter("SOL_VECTOR_INIT_VALUE");
  
  //       Vec sol_vec = dynamic_cast<PetscVector<double>&>(vector).vec();
  //       VecSet(sol_vec, init_value);
  //     }
  //   else if (method == "VECTOR")
  //     {
  //       // this needs to be implemented.
  //       Assert(false, ExcInternalError());
  //       //     std::string vec_name = 
  //       //     this->fesystem_controller.analysis_case->getStringParameter("SOL_VECTOR_INIT_VECTOR_FILE");
  //       //     unsigned int load_case = this->getCurrentLoadCase();
  //       //     this->fesystem_controller.global_data_storage->fillVector(load_case, vec_name, vector);
  //     }
  
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      // the Initial conditions for now, are being hardcoded as zero, but will be changesd later
      for (unsigned int i=0; i < order; i++)
        {
          vectors[i]->zero();
          vectors[i]->perform_intermediate_vector_assembly();
          if (i==0)
            this->applyBoundaryConditionsToVectorForCurrentAnalysis(*(vectors[i]), true, it->second->getDisciplineEnumID());
          else
            this->applyBoundaryConditionsToVectorForCurrentAnalysis(*(vectors[i]), false, it->second->getDisciplineEnumID());
          
          vectors[i]->close();
        }
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
}




Driver::LinearTransientAnalysisDriver::LinearTransientAnalysisDriver
(const unsigned int ID, FESystem::FESystemController& fesys_controller):
Driver::TransientAnalysisDriver(ID, fesys_controller, 
                                LINEAR_TRANSIENT_ANALYSIS_DRIVER::num())
{
  
}



Driver::LinearTransientAnalysisDriver::~LinearTransientAnalysisDriver()
{
  
}


void
Driver::LinearTransientAnalysisDriver::initialize()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  // add the matrices and vectors
  this->addMatricesAndVectors();
}



void Driver::LinearTransientAnalysisDriver::addMatricesAndVectors()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // add a system matrix, a solution vector and a rhs
  switch (it->second->getTransientSystemOrder())
  {
    case 2:
    {
      if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
        {
          SparseMatrix<double>& mat = this->addMatrix(Driver::TRANSIENT_C2_MATRIX::name());
          it->second->initMatrix(mat);
        }
      
    }
      
      
    case 1:
    {
      if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
        {
          SparseMatrix<double>& mat = this->addMatrix(Driver::TRANSIENT_C1_MATRIX::name());
          it->second->initMatrix(mat);
        }
      
      if (it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
        {
          SparseMatrix<double>& mat = this->addMatrix(Driver::SYSTEM_MATRIX::name());
          it->second->initMatrix(mat);
        }
    }
      break;
      
    default:
      Assert(false, ExcInternalError()); 
  }
}






void Driver::LinearTransientAnalysisDriver::solveCurrentLoadCase()
{
  Assert(this->time_values.size() == 0, ExcInternalError());
  
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  unsigned int order = it->second->getTransientSystemOrder();
  Assert(order > 0, ExcInternalError());
  Assert(this->solver != NULL, ExcInternalError());
  
  Solver::LinearTransientSolverBase* linear_transient_solver = 
  dynamic_cast<Solver::LinearTransientSolverBase*>(this->solver);
  
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      FESystemUtility::AutoPtrVector<NumericVector<double> > vecs(order);
      
      std::auto_ptr<NumericVector<double> > vec_ptr;
      for (unsigned int i=0; i < order; i++)
        {
          // create a vector for the initial conditions, and init it to the 
          // dof distribution of the discipline
          vec_ptr.reset(NumericVector<double>::build().release());
          it->second->initVector(*vec_ptr);
          vecs.reset(i, vec_ptr.release());
        }
      
      this->setInitialConditionForCurrentLoadCase(vecs.getReference());
      linear_transient_solver->setInitialCondition(vecs.getReference());
      
      // now clear the vectors, since they are not needed anymore
      vecs.clear();
      
      
      // create a vector of system matrices to be given to the solver
      std::vector<SparseMatrix<double>*> mat_vec(order+1);
      
      switch (order)
      {
        case 2:
        {
          if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
            mat_vec[2] = &(this->getMatrix(Driver::TRANSIENT_C2_MATRIX::name()));
          else 
            mat_vec[2] = NULL;
        }
          
          
        case 1:
        {
          if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
            mat_vec[1] = &(this->getMatrix(Driver::TRANSIENT_C1_MATRIX::name()));
          else 
            mat_vec[1] = NULL;
          
          if (it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
            mat_vec[0] = &(this->getMatrix(Driver::SYSTEM_MATRIX::name()));
          else
            mat_vec[0] = NULL;
        }
          break;
          
        default:
          Assert(false, ExcInternalError()); 
      }
      
      linear_transient_solver->attachCoefficientMatrices(mat_vec);
      
      linear_transient_solver->solve();
    }
      break;
      
    default:
      abort();
      break;      
  }
  
  // once the solution is done, store the time values
  if (this->current_analysis_kind_enum_ID == ANALYSIS::num())
    {
      FESystemDatabase::TimeIndependentDataInfo data_info;
      
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("TimeValues");
      data_info.setLoadCase(this->current_load_case);
      
      FESystemNumerics::PetscSeqVector<double> vec(this->time_values.size());
      for (unsigned int i=0; i < this->time_values.size(); i++)
        vec.set(i, this->time_values[i]);
      
      this->fesystem_controller.global_data_storage->storeVector(data_info, vec);
    }
}




void 
Driver::LinearTransientAnalysisDriver::getCoefficientMatrices
(const double time, std::vector<SparseMatrix<double>* >& matrices)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator discipline_it;
  discipline_it = this->analysis_discipline_map.begin();

  this->current_analysis_time = time;
  
  unsigned int order = discipline_it->second->getTransientSystemOrder();
  Assert(matrices.size() == (order+1), ExcInternalError());
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (order)
  {
    case 2:
    {
      if (discipline_it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
        {
          matrices[2]->zero(); matrices[2]->close();
          matrix_map[Driver::TRANSIENT_C2_MATRIX::num()] = matrices[2];
        }
    }
      
      
    case 1:
    {
      if (discipline_it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
        {
          matrices[1]->zero(); matrices[1]->close();
          matrix_map[Driver::TRANSIENT_C1_MATRIX::num()] = matrices[1];	
        }
      
      if (discipline_it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
        {
          matrices[0]->zero(); matrices[0]->close();
          matrix_map[Driver::SYSTEM_MATRIX::num()] = matrices[0];
        }
    }
      break;
	    
    default:
      Assert(false, ExcInternalError()); 
  }
  
  // fill the matrices and vectors with the required quantities
  discipline_it->second->fillQuantity(real_map, matrix_map, vector_map);
  
  std::map<unsigned int, SparseMatrix<double>*>::iterator it, end;
  it = matrix_map.begin();
  end = matrix_map.end();
  
  double time_step = dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getTimeStepSize();
  unsigned int i = 0;
  
  for (; it != end; it++)
    {
      it->second->close();
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis(*(it->second), pow(time_step,i), discipline_it->second->getDisciplineEnumID());
      it->second->close();
      i++;
    }

  // return control back to the solver
}



void 
Driver::LinearTransientAnalysisDriver::getForceVector
(const double time, NumericVector<double>& vec)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  this->current_analysis_time = time;
  
  unsigned int order = it->second->getTransientSystemOrder();
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      vec.zero(); vec.close();
      vector_map[Driver::FORCE_VECTOR::num()] = &vec;
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      vec.close();
      this->applyBoundaryConditionsToVectorForCurrentAnalysis(vec, true, it->second->getDisciplineEnumID());
      vec.close();
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      if (!this->haveMatrix("ScratchMatrix"))
        this->addMatrix("ScratchMatrix");
      if (!this->haveVector("ScratchVector1"))
        this->addVector("ScratchVector1");
      if (!this->haveVector("ScratchVector2"))
        this->addVector("ScratchVector2");
        
      SparseMatrix<double>& scratch_matrix = this->getMatrix("ScratchMatrix");
      NumericVector<double>& scratch_vec1 = this->getVector("ScratchVector1");
      NumericVector<double>& scratch_vec2 = this->getVector("ScratchVector2");
      
      vec.zero(); vec.perform_intermediate_vector_assembly();
      vector_map[Driver::FORCE_VECTOR_SENSITIVITY::num()] = &vec;
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      vector_map.clear();
      vec.perform_intermediate_vector_assembly();
      
      FESystemDatabase::TimeDependentDataInfo data_info;
      
      unsigned int qty_enum_ID = 0;
      for (unsigned int i=0; i<= order; i++)
        {
          switch (i)
          {
            case 0:
              qty_enum_ID = Driver::SYSTEM_MATRIX_SENSITIVITY::num();
              break;
              
            case 1:
              qty_enum_ID = Driver::TRANSIENT_C1_MATRIX_SENSITIVITY::num();
              break;
              
            case 2:
              qty_enum_ID = Driver::TRANSIENT_C2_MATRIX_SENSITIVITY::num();
              break;		
              
            default:
              Assert(false, ExcInternalError());
          }
          
          // if the matrix for this order exists, then process it, else go to the 
          // next order
          if (it->second->disciplineHasMatrix(qty_enum_ID))
            {
              // get the matrix and the solution
              scratch_matrix.zero(); scratch_matrix.perform_intermediate_matrix_assembly();
              matrix_map.clear();
              matrix_map[qty_enum_ID] = &scratch_matrix;
              it->second->fillQuantity(real_map, matrix_map, vector_map);
              scratch_matrix.close();
              
              // create the data info object for this solution
              data_info.clear();
              data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
              data_info.setName("Solution");
              data_info.setOrder(i);
              data_info.setLoadCase(this->current_load_case); 
              // this needs to be implemted, since iter_num is not available here
              Assert(false, ExcInternalError());
              data_info.setTransientIterationInfo(this->getCurrentIterationNumber(), time); 
              
              this->fesystem_controller.global_data_storage->fillVector(data_info, scratch_vec1);
              scratch_vec1.close();
              scratch_vec2.zero(); 
              scratch_vec2.close();
              
              scratch_matrix.multiply_vector(scratch_vec1, scratch_vec2);
              vec.add(-1.0, scratch_vec2);
              vec.perform_intermediate_vector_assembly();
            }
        }
      
      this->applyBoundaryConditionsToVectorForCurrentAnalysis(vec, false, it->second->getDisciplineEnumID());
      vec.close();
    }
      break;
      
    default:
      Assert(false, ExcInternalError());
  }
  
  // return control back to the solver
}




Driver::NonlinearTransientAnalysisDriver::NonlinearTransientAnalysisDriver
(const unsigned int ID, FESystem::FESystemController& fesys_controller):
Driver::TransientAnalysisDriver(ID, fesys_controller, 
                                NONLINEAR_TRANSIENT_ANALYSIS_DRIVER::num())
{
  
}



Driver::NonlinearTransientAnalysisDriver::~NonlinearTransientAnalysisDriver()
{
  
}




void
Driver::NonlinearTransientAnalysisDriver::initialize()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  // now resize the current solution vector
  Assert(this->current_solution.size() == 0,ExcInternalError());
  this->current_solution.resize(it->second->getTransientSystemOrder()+1);
  for (unsigned int i=0; i < this->current_solution.size(); i++)
    {
      this->current_solution.reset(i, NumericVector<double>::build().release());
      this->current_solution[i]->init(it->second->getAnalysisDofMap().n_dofs());
    }
  
  // add the matrices and vectors
  this->addMatricesAndVectors();
}



void Driver::NonlinearTransientAnalysisDriver::addMatricesAndVectors()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // add the jacobian matrix
  SparseMatrix<double>& jac = this->addMatrix(Driver::JACOBIAN_MATRIX::name());
  it->second->initMatrix(jac);
  SparseMatrix<double>& scr_mat = this->addMatrix("ScratchMatrix");
  it->second->initMatrix(scr_mat);
  NumericVector<double>& scr_vec = this->addVector("ScratchVector1");
  it->second->initVector(scr_vec);
  
  
  // add a system matrix, a solution vector and a rhs
  switch (it->second->getTransientSystemOrder())
  {
    case 2:
    {
      if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
        {
          SparseMatrix<double>& mat = this->addMatrix(Driver::TRANSIENT_C2_MATRIX::name());
          it->second->initMatrix(mat);
        }
      
    }
      
      
    case 1:
    {
      if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
        {
          SparseMatrix<double>& mat = this->addMatrix(Driver::TRANSIENT_C1_MATRIX::name());
          it->second->initMatrix(mat);
        }
      
      if (it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
        {
          SparseMatrix<double>& mat = this->addMatrix(Driver::SYSTEM_MATRIX::name());
          it->second->initMatrix(mat);
        }
    }
      break;
      
    default:
      Assert(false, ExcInternalError()); 
  }
}






void Driver::NonlinearTransientAnalysisDriver::solveCurrentLoadCase()
{
  
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  unsigned int order = it->second->getTransientSystemOrder();
  Assert(order > 0, ExcInternalError());
  // make sure that the solver is the right kind
  Assert(this->solver != NULL, ExcInternalError());
    
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      Assert(this->time_values.size() == 0, ExcInternalError());

      Assert(this->solver->getSolverClassEnumID() == 
             Solver::NONLINEAR_TRANSIENT_SOLVER::num(), ExcInternalError());

      Solver::NonlinearTransientSolverBase* nonlinear_transient_solver = 
      dynamic_cast<Solver::NonlinearTransientSolverBase*>(this->solver);

      
      FESystemUtility::AutoPtrVector<NumericVector<double> > vecs(order);
      
      std::auto_ptr<NumericVector<double> > vec_ptr;
      for (unsigned int i=0; i < order; i++)
        {
          // create a vector for the initial conditions, and init it to the 
          // dof distribution of the discipline
          vec_ptr.reset(NumericVector<double>::build().release());
          it->second->initVector(*vec_ptr);
          vecs.reset(i, vec_ptr.release());
        }
      
      this->setInitialConditionForCurrentLoadCase(vecs.getReference());
      nonlinear_transient_solver->setInitialCondition(vecs.getReference());
      
      SparseMatrix<double>& jac_matrix = this->getMatrix(Driver::JACOBIAN_MATRIX::name());
      jac_matrix.zero();
      nonlinear_transient_solver->attachJacobianMatrix(jac_matrix);
      
      // now clear the vectors, since they are not needed anymore
      vecs.clear();
      
      
      // create a vector of system matrices to be given to the solver
      std::vector<SparseMatrix<double>*> mat_vec(order+1);
      
      switch (order)
      {
        case 2:
        {
          if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
            {
              mat_vec[2] = &(this->getMatrix(Driver::TRANSIENT_C2_MATRIX::name()));
              mat_vec[2]->zero();
            }
          else 
            mat_vec[2] = NULL;
        }
          
          
        case 1:
        {
          if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
            {
              mat_vec[1] = &(this->getMatrix(Driver::TRANSIENT_C1_MATRIX::name()));
              mat_vec[1]->zero();
            }
          else 
            mat_vec[1] = NULL;
          
          if (it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
            {
              mat_vec[0] = &(this->getMatrix(Driver::SYSTEM_MATRIX::name()));
              mat_vec[0]->zero();
            }
          else
            mat_vec[0] = NULL;
        }
          break;
          
        default:
          Assert(false, ExcInternalError()); 
      }
      
      nonlinear_transient_solver->attachCoefficientMatrices(mat_vec);

      // now that the matrices have been set and initialized, solve the system
      nonlinear_transient_solver->solve();
      
      // after the solution, store the time values at which solution was obtained
      if (FESystem::local_processor == 0)
        {
          FESystemDatabase::TimeIndependentDataInfo data_info;
          
          data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
          data_info.setName("TimeValues");
          data_info.setLoadCase(this->current_load_case);
          
          FESystemNumerics::PetscSeqVector<double> vec(this->time_values.size());
          for (unsigned int i=0; i < this->time_values.size(); i++)
            vec.set(i, this->time_values[i]);
          
          this->fesystem_controller.global_data_storage->storeVector(data_info, vec);
        }
      
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      // the time values should already be available. If not, load them from disc
      if (this->time_values.size() == 0)
        {
          FESystemDatabase::TimeIndependentDataInfo data_info;
          
          data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
          data_info.setName("TimeValues");
          data_info.setLoadCase(this->current_load_case);
          
          FESystemNumerics::PetscSeqVector<double> vec;
          this->fesystem_controller.global_data_storage->fillVector(data_info, vec);
          this->time_values.resize(vec.size());
          for (unsigned int i=0; i < vec.size(); i++)
            this->time_values[i] = vec.el(i);
        }
      
      
      // get a pointer to the solver, and make sure it is of the right kind
      Assert(this->solver->getSolverClassEnumID() == 
             Solver::LINEAR_TRANSIENT_SOLVER::num(), ExcInternalError());
      
      Solver::LinearTransientSolverBase* linear_transient_solver = 
      dynamic_cast<Solver::LinearTransientSolverBase*>(this->solver);
      
      
      FESystemUtility::AutoPtrVector<NumericVector<double> > vecs(order);
      
      std::auto_ptr<NumericVector<double> > vec_ptr;
      for (unsigned int i=0; i < order; i++)
        {
          // create a vector for the initial conditions, and init it to the 
          // dof distribution of the discipline
          vec_ptr.reset(NumericVector<double>::build().release());
          it->second->initVector(*vec_ptr);
          vecs.reset(i, vec_ptr.release());
        }
      
      // all sensitivity analysis solutions are assumed to have a zero 
      // initial condition
      this->setInitialConditionForCurrentLoadCase(vecs.getReference());
      linear_transient_solver->setInitialCondition(vecs.getReference());

      // now clear the vectors, since they are not needed anymore
      vecs.clear();
      
      
      // create a vector of system matrices to be given to the solver
      std::vector<SparseMatrix<double>*> mat_vec(order+1);
      
      // create the matrices and attach them to the solver
      switch (order)
      {
        case 2:
        {
          if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
            {
              mat_vec[2] = &(this->getMatrix(Driver::TRANSIENT_C2_MATRIX::name()));
              mat_vec[2]->zero();
            }
          else 
            mat_vec[2] = NULL;
        }
          
          
        case 1:
        {
          if (it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
            {
              mat_vec[1] = &(this->getMatrix(Driver::TRANSIENT_C1_MATRIX::name()));
              mat_vec[1]->zero();
            }
          else 
            mat_vec[1] = NULL;
          
          if (it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
            {
              mat_vec[0] = &(this->getMatrix(Driver::SYSTEM_MATRIX::name()));
              mat_vec[0]->zero();
            }
          else
            mat_vec[0] = NULL;
        }
          break;
          
        default:
          Assert(false, ExcInternalError()); 
      }
      
      linear_transient_solver->attachCoefficientMatrices(mat_vec);
      
      // now that the matrices have been set and initialized, solve the system
      linear_transient_solver->solve();
    }
      break;
      
    default:
      abort();
      break;      
  }
  
}




void 
Driver::NonlinearTransientAnalysisDriver::getCoefficientMatrices
(const double time,
 std::vector<NumericVector<double>* >& vecs,                                                                 
 std::vector<SparseMatrix<double>* >& matrices)
{  
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator discipline_it;
  discipline_it = this->analysis_discipline_map.begin();

  this->current_analysis_time = time;
  
  unsigned int order = discipline_it->second->getTransientSystemOrder();
  
  Assert(vecs.size() == (order+1), ExcInternalError());
  Assert(matrices.size() == (order+1), ExcInternalError());
  Assert(this->current_solution.size() == (order+1),ExcInternalError());
  
  
  for (unsigned int i=0; i <= order; i++)
    {
      Assert(vecs[i] != NULL, ExcInternalError());
      this->current_solution[i]->zero();
      vecs[i]->close();
      this->localizeSolution(discipline_it->second->getAnalysisDofMap(), *(vecs[i]),
                             *(this->current_solution[i]));
      this->current_solution[i]->close();
    }
  
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (order)
  {
    case 2:
    {
      if (matrices[2] != NULL)
        {
          Assert(discipline_it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()),
                 ExcInternalError());
          matrices[2]->zero(); matrices[2]->perform_intermediate_matrix_assembly();
          matrix_map[Driver::TRANSIENT_C2_MATRIX::num()] = matrices[2];
        }
    }
      
      
    case 1:
    {
      if (matrices[1] != NULL)
        {
          Assert(discipline_it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()),
                 ExcInternalError());
          matrices[1]->zero(); matrices[1]->perform_intermediate_matrix_assembly();
          matrix_map[Driver::TRANSIENT_C1_MATRIX::num()] = matrices[1];
        }
      
      
      if (matrices[0] != NULL)
        {
          Assert(discipline_it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()),
                 ExcInternalError());
          matrices[0]->zero(); matrices[0]->perform_intermediate_matrix_assembly();
          matrix_map[Driver::SYSTEM_MATRIX::num()] = matrices[0];
        }
    }
      break;
	    
    default:
      Assert(false, ExcInternalError()); 
  }
  
  // fill the matrices and vectors with the required quantities
  discipline_it->second->fillQuantity(real_map, matrix_map, vector_map);
  
  std::map<unsigned int, SparseMatrix<double>*>::iterator it, end;
  it = matrix_map.begin();
  end = matrix_map.end();
  
  double time_step = dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getTimeStepSize();
  unsigned int i = 0;

  for (; it != end; it++)
    {
      it->second->perform_intermediate_matrix_assembly();
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis(*(it->second), pow(time_step,i), discipline_it->second->getDisciplineEnumID());
      it->second->close();
      i++;
    }
  
  // clear the solution vectors
  for (unsigned int i=0; i <= order; i++)
    this->current_solution[i]->zero();
  
  // return control back to the solver
}


void 
Driver::NonlinearTransientAnalysisDriver::getForceVector(const double time,
                                                         std::vector<NumericVector<double>* >& vecs,                                                                 
                                                         NumericVector<double>& vec)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  this->current_analysis_time = time;
  
  unsigned int order = it->second->getTransientSystemOrder();
  Assert(vecs.size() == (order+1), ExcInternalError());
  Assert(this->current_solution.size() == (order+1),ExcInternalError());
  
  for (unsigned int i=0; i <= order; i++)
    {
      Assert(vecs[i] != NULL, ExcInternalError());
      this->current_solution[i]->zero();
      vecs[i]->close();
      this->localizeSolution(it->second->getAnalysisDofMap(), *(vecs[i]),
                              *(this->current_solution[i]));
      this->current_solution[i]->close();
    }
  vec.zero();
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      vector_map[Driver::FORCE_VECTOR::num()] = &vec;
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      this->applyBoundaryConditionsToVectorForCurrentAnalysis(vec, true, it->second->getDisciplineEnumID());
      vec.close();
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    default:
      Assert(false, ExcInternalError());
  }  
  
  // clear the solution vectors
  for (unsigned int i=0; i <= order; i++)
    this->current_solution[i]->zero();
  
  // return control back to the solver
}



void 
Driver::NonlinearTransientAnalysisDriver::getResidualVector(const double time,
                                                            std::vector<NumericVector<double>* >& vecs,                                                                 
                                                            NumericVector<double>& vec)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  this->current_analysis_time = time;
  
  unsigned int order = it->second->getTransientSystemOrder();
  Assert(vecs.size() == (order+1), ExcInternalError());
  Assert(this->current_solution.size() == (order+1),ExcInternalError());
  
  for (unsigned int i=0; i <= order; i++)
    {
      Assert(vecs[i] != NULL, ExcInternalError());
      this->current_solution[i]->zero();
      vecs[i]->close();
      this->localizeSolution(it->second->getAnalysisDofMap(), *(vecs[i]),
                              *(this->current_solution[i]));
      this->current_solution[i]->close();
    }
  vec.zero();
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      SparseMatrix<double>& scratch_matrix = this->getMatrix("ScratchMatrix");
      NumericVector<double>& scratch_vec1 = this->getVector("ScratchVector1");

      double time_step = dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getTimeStepSize();

      unsigned int qty_enum_ID = 0;
      for (unsigned int i=0; i<= order; i++)
        {
          switch (i)
          {
            case 0:
              qty_enum_ID = Driver::SYSTEM_MATRIX::num();
              break;
              
            case 1:
              qty_enum_ID = Driver::TRANSIENT_C1_MATRIX::num();
              break;
              
            case 2:
              qty_enum_ID = Driver::TRANSIENT_C2_MATRIX::num();
              break;		
              
            default:
              Assert(false, ExcInternalError());
          }
          
          // if the matrix for this order exists, then process it, else go to the 
          // next order
          if (it->second->disciplineHasMatrix(qty_enum_ID))
            {
              // get the matrix and the solution
              matrix_map.clear();
              matrix_map[qty_enum_ID] = &scratch_matrix;
              scratch_matrix.zero();
              scratch_matrix.perform_intermediate_matrix_assembly();
              it->second->fillQuantity(real_map, matrix_map, vector_map);

              // apply the boundary conditions to this matrix
              scratch_matrix.close();
              this->applyBoundaryConditionsToMatrixForCurrentAnalysis(scratch_matrix,
                                                                      pow(time_step, i), it->second->getDisciplineEnumID());
              
              scratch_matrix.close();
              scratch_vec1.zero(); 
              scratch_vec1.close();
              
              scratch_matrix.multiply_vector(*(vecs[i]), scratch_vec1);
              scratch_vec1.close();
              
              vec.add(1.0, scratch_vec1);
              vec.perform_intermediate_vector_assembly();
            }
        }
      
      scratch_vec1.zero();
      scratch_vec1.close();
      matrix_map.clear();
      vector_map.clear();
      vector_map[Driver::FORCE_VECTOR::num()] = &scratch_vec1;
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      scratch_vec1.close();
      this->applyBoundaryConditionsToVectorForCurrentAnalysis(scratch_vec1, true, it->second->getDisciplineEnumID());
      scratch_vec1.close();
      vec.add(-1.0, scratch_vec1);
      vec.close();
      
      // this has been temporarily removed to check the approach to application of 
      // dirichlet BCs. The BCs are applied to the coefficient matrices instead of this
      //       this->applyBoundaryConditionsToVectorForCurrentAnalysis(vec, false);
      
    }
      break;
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    default:
      Assert(false, ExcInternalError());
  }  
  
  // clear the solution vectors
  for (unsigned int i=0; i <= order; i++)
    this->current_solution[i]->zero();
  
  // return control back to the solver
}



void 
Driver::NonlinearTransientAnalysisDriver::getJacobianMatrix(const double time,
                                                            std::vector<NumericVector<double>* >& vecs,                                                                 
                                                            SparseMatrix<double>& jac)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  this->current_analysis_time = time;
  
  unsigned int order = it->second->getTransientSystemOrder();
  Assert(vecs.size() == (order+1), ExcInternalError());
  Assert(this->current_solution.size() == (order+1),ExcInternalError());
  
  for (unsigned int i=0; i <= order; i++)
    {
      Assert(vecs[i] != NULL, ExcInternalError());
      this->current_solution[i]->zero();
      vecs[i]->close();
      this->localizeSolution(it->second->getAnalysisDofMap(), *(vecs[i]),
                              *(this->current_solution[i]));
      this->current_solution[i]->close();
    }
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      jac.zero(); jac.perform_intermediate_matrix_assembly();
      matrix_map[Driver::JACOBIAN_MATRIX::num()] = &jac;
      
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      jac.perform_intermediate_matrix_assembly();
      
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis(jac, 1.0, it->second->getDisciplineEnumID());
      jac.close();
    }
      break;
      
      // jacobian for a sensitivity analysis does not make sense, since the 
      // solution is already known
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    default:
      Assert(false, ExcInternalError());
  }
  
  // clear the solution vectors
  for (unsigned int i=0; i <= order; i++)
    this->current_solution[i]->zero();
  
  // return control back to the solver
}


void 
Driver::NonlinearTransientAnalysisDriver::getCoefficientMatrices
(const double time, std::vector<SparseMatrix<double>* >& matrices)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator discipline_it;
  discipline_it = this->analysis_discipline_map.begin();
  

  // this routine is to be used only for sensitivity analysis routines, which are linear
  Assert(this->current_analysis_kind_enum_ID == SENSITIVITY_ANALYSIS::num(),
         ExcInternalError());
  
  // set the analysis time, and make sure that it matches with the time in the time vector
  // for which the analysis was performed
  this->current_analysis_time = time;
  Assert(this->getCurrentAnalysisTime() ==
         this->time_values[this->getCurrentIterationNumber()-1],
         ExcInternalError());
  
  unsigned int order = discipline_it->second->getTransientSystemOrder();
  Assert(matrices.size() == (order+1), ExcInternalError());
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  switch (order)
  {
    case 2:
    {
      if (discipline_it->second->disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
        {
          matrices[2]->zero(); matrices[2]->close();
          matrix_map[Driver::TRANSIENT_C2_MATRIX::num()] = matrices[2];
        }
    }
      
      
    case 1:
    {
      if (discipline_it->second->disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
        {
          matrices[1]->zero(); matrices[1]->close();
          matrix_map[Driver::TRANSIENT_C1_MATRIX::num()] = matrices[1];	
        }
      
      if (discipline_it->second->disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
        {
          matrices[0]->zero(); matrices[0]->close();
          // ask for the jacobian matrix, since that is needed for nonlinear problems
          matrix_map[Driver::JACOBIAN_MATRIX::num()] = matrices[0];
        }
    }
      break;
	    
    default:
      Assert(false, ExcInternalError()); 
  }
  
  // set the solution vectors for the current iteration
  Assert(this->current_solution.size() == (order+1), 
         ExcInternalError());

  FESystemUtility::AutoPtrVector<NumericVector<double> > vecs(order+1);
  
  std::auto_ptr<NumericVector<double> > vec_ptr;
  for (unsigned int i=0; i <= order; i++)
    {
      // create a vector for the initial conditions, and init it to the 
      // dof distribution of the discipline
      vec_ptr.reset(NumericVector<double>::build().release());
      discipline_it->second->initVector(*vec_ptr);
      vecs.reset(i, vec_ptr.release());
    }

  FESystemDatabase::TimeDependentDataInfo data_info;

  // now get the solution vectors from the database
  for (unsigned int i=0; i <= order; i++)
    {
      data_info.clear();
      data_info.setDisciplineEnumID(discipline_it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(i);
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo(this->getCurrentIterationNumber(), 
                                          this->getCurrentAnalysisTime());
      
      this->fesystem_controller.global_data_storage->fillVector(data_info, *(vecs[i]));
      vecs[i]->close();
      this->current_solution[i]->zero();
      this->localizeSolution(discipline_it->second->getAnalysisDofMap(),*(vecs[i]),
                             *(this->current_solution[i]));
    }
  
  
  // fill the matrices and vectors with the required quantities
  discipline_it->second->fillQuantity(real_map, matrix_map, vector_map);
  
  
  std::map<unsigned int, SparseMatrix<double>*>::iterator it, end;
  it = matrix_map.begin();
  end = matrix_map.end();
  
  double time_step = dynamic_cast<Solver::TransientSolverBase*>(this->solver)->getTimeStepSize();
  unsigned int i = 0;

  for (; it != end; it++)
    {
      it->second->close();
      this->applyBoundaryConditionsToMatrixForCurrentAnalysis(*(it->second), pow(time_step,i), discipline_it->second->getDisciplineEnumID());
      it->second->close();
      i++;
    }
  
  // clear the solution vectors before returning
  for (unsigned int i=0; i <= order; i++)
    this->current_solution[i]->zero();
  
  // return control back to the solver
}


void 
Driver::NonlinearTransientAnalysisDriver::getForceVector
(const double time, NumericVector<double>& vec)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  // this routine is to be used only for sensitivity analysis routines, which are linear
  Assert(this->current_analysis_kind_enum_ID == SENSITIVITY_ANALYSIS::num(),
         ExcInternalError());
  
  // set the analysis time, and make sure that it matches with the time in the time vector
  // for which the analysis was performed
  this->current_analysis_time = time;
  Assert(this->getCurrentAnalysisTime() ==
         this->time_values[this->getCurrentIterationNumber()-1],
         ExcInternalError());
  
  unsigned int order = it->second->getTransientSystemOrder();
    
  // set the solution vectors for the current iteration
  Assert(this->current_solution.size() == (order+1), ExcInternalError());
  FESystemUtility::AutoPtrVector<NumericVector<double> > vecs(order+1);
  
  std::auto_ptr<NumericVector<double> > vec_ptr;
  for (unsigned int i=0; i <= order; i++)
    {
      // create a vector for the initial conditions, and init it to the 
      // dof distribution of the discipline
      vec_ptr.reset(NumericVector<double>::build().release());
      it->second->initVector(*vec_ptr);
      vecs.reset(i, vec_ptr.release());
    }
  
  FESystemDatabase::TimeDependentDataInfo data_info;
  
  // now get the solution vectors from the database
  for (unsigned int i=0; i <= order; i++)
    {
      data_info.clear();
      data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      data_info.setName("Solution");
      data_info.setOrder(i);
      data_info.setLoadCase(this->current_load_case);
      data_info.setTransientIterationInfo(this->getCurrentIterationNumber(), 
                                          this->getCurrentAnalysisTime());
      
      this->fesystem_controller.global_data_storage->fillVector(data_info, *(vecs[i]));
      vecs[i]->close();
      this->current_solution[i]->zero();
      this->localizeSolution(it->second->getAnalysisDofMap(),*(vecs[i]),
                             *(this->current_solution[i]));
    }
  
  // create a map of the matrices to be requested from the system
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  if (!this->haveMatrix("ScratchMatrix"))
    this->addMatrix("ScratchMatrix");
  if (!this->haveVector("ScratchVector1"))
    this->addVector("ScratchVector1");
  
  SparseMatrix<double>& scratch_matrix = this->getMatrix("ScratchMatrix");
  NumericVector<double>& scratch_vec1 = this->getVector("ScratchVector1");
  
  vec.zero(); vec.perform_intermediate_vector_assembly();
  vector_map[Driver::FORCE_VECTOR_SENSITIVITY::num()] = &vec;
  it->second->fillQuantity(real_map, matrix_map, vector_map);
  vector_map.clear();
  vec.perform_intermediate_vector_assembly();
  
  unsigned int qty_enum_ID = 0, qty_sensitivity_enum_ID=0;
  for (unsigned int i=0; i<= order; i++)
    {
      switch (i)
      {
        case 0:
        {
          qty_enum_ID = Driver::SYSTEM_MATRIX::num();
          qty_sensitivity_enum_ID = Driver::SYSTEM_MATRIX_SENSITIVITY::num();
        }
          break;
          
        case 1:
        {
          qty_enum_ID = Driver::TRANSIENT_C1_MATRIX::num();
          qty_sensitivity_enum_ID = Driver::TRANSIENT_C1_MATRIX_SENSITIVITY::num();
        }
          break;
          
        case 2:
        {
          qty_enum_ID = Driver::TRANSIENT_C2_MATRIX::num();
          qty_sensitivity_enum_ID = Driver::TRANSIENT_C2_MATRIX_SENSITIVITY::num();
        }
          break;		
          
        default:
          Assert(false, ExcInternalError());
      }
      
      // if the matrix for this order exists, then process it, else go to the 
      // next order
      if (it->second->disciplineHasMatrix(qty_enum_ID))
        {
          // get the matrix and the solution
          scratch_matrix.zero(); scratch_matrix.perform_intermediate_matrix_assembly();
          matrix_map.clear();
          matrix_map[qty_sensitivity_enum_ID] = &scratch_matrix;
          it->second->fillQuantity(real_map, matrix_map, vector_map);
          scratch_matrix.close();
          
          scratch_vec1.zero(); 
          
          scratch_matrix.multiply_vector(*(this->current_solution[i]), scratch_vec1);
          vec.add(-1.0, scratch_vec1);
          vec.perform_intermediate_vector_assembly();
        }
    }
  
  this->applyBoundaryConditionsToVectorForCurrentAnalysis(vec, false, it->second->getDisciplineEnumID());
  vec.close();
  
  // clear the solution vectors before returning
  for (unsigned int i=0; i <= order; i++)
    this->current_solution[i]->zero();
  
  // return control back to the solver
}
