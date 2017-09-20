// $Id: EigenProblemAnalysisDriver.C,v 1.13.4.8 2008/08/21 20:44:25 manav Exp $

// C++ includes


// FESystem includes
#include "AnalysisDriver/EigenProblemAnalysisDriver.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Solvers/EigenSolver.h"
#include "Solvers/EigenSolverInfo.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/DisciplineInfo.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"
#include "Solutions/SolutionBase.h"
#include "Utilities/Log.h"

// libMesh includes
#include "base/dof_map.h"
#include "geom/node.h"
#include "numerics/sparse_matrix.h"
#include "numerics/petsc_matrix.h"

int slepc_eig_solve(Mat A, Mat B);


Driver::EigenProblemAnalysisDriver::EigenProblemAnalysisDriver
(const unsigned int ID,
 FESystem::FESystemController& fesys_controller):
Driver::AnalysisDriver(ID, fesys_controller, EIGENPROBLEM_ANALYSIS_DRIVER::num())
{
}





Driver::EigenProblemAnalysisDriver::~EigenProblemAnalysisDriver()
{
  
}



void
Driver::EigenProblemAnalysisDriver::initialize()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  // add the matrices and vectors
  this->addMatricesAndVectors();
}



void
Driver::EigenProblemAnalysisDriver::addMatricesAndVectors()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  // before adding the matrices, the driver should check with the discipline
  // about the kind of problem. If the proble is generlized eigenvalue problem, 
  // then we will need two matrices. Else, only one.
  
  Assert(it->second != NULL, ExcInternalError())
  unsigned int kind_enum_ID = it->second->getEigenProblemKindEnumID();
  
  // add the matrices to the driver
  it->second->initMatrix(this->addMatrix("A_matrix"));
  
  if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
      kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
    it->second->initMatrix(this->addMatrix("B_matrix"));
  
  it->second->initVector(this->addVector("EigenVector_Real"));
  it->second->initVector(this->addVector("EigenVector_Img"));
  it->second->initVector(this->addVector("ScratchVec"));
}





void
Driver::EigenProblemAnalysisDriver::solveCurrentLoadCase()
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  SparseMatrix<double>& A_matrix = this->getMatrix("A_matrix");
  SparseMatrix<double>* B_matrix = NULL;
  
  unsigned int kind_enum_ID = it->second->getEigenProblemKindEnumID();
  if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
      kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
    B_matrix = &(this->getMatrix("B_matrix"));
  
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  real_map[Driver::MODEL_MASS::num()] = 0.0;
  // prepare the load vector for the analysis and then init the solver and then finally 
  // solve the system
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      Solver::EigenSolverBase* eigen_solver = dynamic_cast<Solver::EigenSolverBase*>(this->solver);
      
      matrix_map[Driver::EIGENPROBLEM_A_MATRIX::num()] = &A_matrix;
      
      if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
          kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
        matrix_map[Driver::EIGENPROBLEM_B_MATRIX::num()] = B_matrix;
      
      it->second->fillQuantity(real_map, matrix_map, vector_map);
      
      std::ostringstream oss;
      oss << std::setprecision(14) << real_map[Driver::MODEL_MASS::num()] << std::endl;
      FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, oss.str());
      
      
      this->applyBoundaryConditionsForCurrentAnalysis(A_matrix, B_matrix);
      
      A_matrix.close();
      if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
          kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
        B_matrix->close();
      
      //	// store the matrices. This is only temporary, and will be removed.
      //	FESystemDatabase::TimeIndependentDataInfo data_info;
      //	data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      //	data_info.setName("A_mat");
      //	data_info.setLoadCase(this->current_load_case);
      //        this->fesystem_controller.global_data_storage->storeMatrix(data_info, A_matrix);
      //
      //	data_info.clear();
      //	data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
      //	data_info.setName("B_mat");
      //	data_info.setLoadCase(this->current_load_case);
      //        this->fesystem_controller.global_data_storage->storeMatrix(data_info, *B_matrix);
      //        
      
      //  eigen_solver->solve(&A_matrix, B_matrix);
      eigen_solver->solve(&A_matrix, B_matrix);
      
      
      //  std::auto_ptr<SparseMatrix<double> > A_matrix_sub(SparseMatrix<double>::build().release() );
      //  std::auto_ptr<SparseMatrix<double> > B_matrix_sub(SparseMatrix<double>::build().release() );
      //
      //  this->getSubmatrixAfterBoundaryCondition(A_matrix, B_matrix, *(A_matrix_sub.get()), B_matrix_sub.get());
      
      
      FESystemDatabase::TimeIndependentEigenDataInfo eigen_vec_info_real, eigen_vec_info_img;
      double value_r, value_i;
      NumericVector<double>& real_vec = this->getVector("EigenVector_Real");
      NumericVector<double>& img_vec = this->getVector("EigenVector_Img");
      std::string load_case_string = "Load Case : ";
      
      std::ostringstream lc_oss;
      lc_oss << this->getCurrentLoadCase();
      load_case_string += lc_oss.str();
      
      unsigned int n_converged = eigen_solver->getNConvergedEigenPairs();
      DenseVector<double> eig_value_real(n_converged), eig_value_img(n_converged);
      
      // also, save the eigen vectors
      FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, n_converged);
      for (unsigned int i=0; i < n_converged; i++)
        {
          std::ostringstream oss;
          oss << i;
          eigen_vec_info_real.clear();
          eigen_vec_info_real.setDisciplineEnumID(it->second->getDisciplineEnumID());
          eigen_vec_info_real.setName("EigenVector_Real");
          eigen_vec_info_real.setLoadCase(this->getCurrentLoadCase());
          eigen_vec_info_real.setModeInfo(i+1);
          
          eigen_vec_info_img.clear();
          eigen_vec_info_img.setDisciplineEnumID(it->second->getDisciplineEnumID());
          eigen_vec_info_img.setName("EigenVector_Img");
          eigen_vec_info_img.setLoadCase(this->getCurrentLoadCase());
          eigen_vec_info_img.setModeInfo(i+1);
          
          eigen_solver->getEigenPair(i, &value_r, &value_i, real_vec, img_vec);
          //          eigen_solver->getEigenValue(i, &value_r, &value_i);
          FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, value_r);//<< " + i " << value_i << std::endl;
          //            << "  residual: " << eigen_solver->getResidualForEigenPair(i) << std::endl;
          
          eig_value_real(i) = value_r;
          eig_value_img(i) = value_i;
          
          
          // now write the solutions to the log file
          if (FESystem::local_processor == 0)
            {
              this->fesystem_controller.log_file->write(load_case_string);
              this->fesystem_controller.log_file->write("EigenSolution");
              std::string value;
              std::ostringstream oss_value;
              value = "EigenValue Number  ";
              value += oss.str();
              value += "   =   ";
              oss_value.setf(std::ios::showpoint);
              oss_value << std::setprecision(14) << value_r;
              if (value_i != 0.0)
                {
                  oss_value << "  + i  ";
                  oss_value << std::setprecision(14) << value_i;
                }
              value += oss_value.str();
              this->fesystem_controller.log_file->write(value);
              this->fesystem_controller.log_file->write("Eigen Vector: Real Part");
            }
              

          this->fesystem_controller.global_data_storage->storeVector(eigen_vec_info_real,real_vec);
          this->writeSolutionVectorToLog(eigen_vec_info_real,
                                         real_vec, it->second->getDisciplineEnumID());
          
          if (value_i != 0.0)
            {
              if (FESystem::local_processor == 0)
                this->fesystem_controller.log_file->write("Eigen Vector: Imaginary Part");

              this->fesystem_controller.global_data_storage->storeVector(eigen_vec_info_img, img_vec);
              this->writeSolutionVectorToLog(eigen_vec_info_img,
                                             img_vec, it->second->getDisciplineEnumID());
            }

          // now register the eigen vector and store it in the database
          FESystemDatabase::TimeIndependentDataInfo data_info;
          data_info.clear();
          data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
          data_info.setName("EigenValuesReal");
          data_info.setLoadCase(this->current_load_case);
          this->fesystem_controller.global_data_storage->storeVector(data_info, eig_value_real);

          data_info.clear();
          data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
          data_info.setName("EigenValuesImg");
          data_info.setLoadCase(this->current_load_case);
          this->fesystem_controller.global_data_storage->storeVector(data_info, eig_value_img);
          
          
        }
      
      
    }
      break;
      
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      // this implementation is only for Hermitian and generalized hermitian cases
      // where all the eigenvalues are real. Cases with complex eigenvalues will 
      // have to be treated specially. 
      
      // it is also assumed that the eigenvectors are normalized with respect to the B 
      // matrix
      
      switch(kind_enum_ID)
      {
        case HERMITIAN_EIGENPROBLEM_ENUM_ID: 
        case GENERALIZED_HERMITIAN_EIGENPROBLEM_ENUM_ID: 
        {
          // get the vector eigenvalues for this load case. 
          DenseVector<double> eig_value_real, eig_value_real_sens;
          FESystemDatabase::TimeIndependentDataInfo data_info;
          data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
          data_info.setName("EigenValuesReal");
          data_info.setLoadCase(this->current_load_case);
          this->fesystem_controller.global_data_storage->fillVector(data_info, eig_value_real);
          eig_value_real_sens.resize(eig_value_real.size());
          
          // get the sensitivity of the system matrices for the current DV
          matrix_map[Driver::EIGENPROBLEM_A_MATRIX_SENSITIVITY::num()] = &A_matrix;
          
          if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
              kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
            {
              matrix_map[Driver::EIGENPROBLEM_B_MATRIX_SENSITIVITY::num()] = B_matrix;
              B_matrix->zero();
            }
          
          A_matrix.zero();
          it->second->fillQuantity(real_map, matrix_map, vector_map);
          
          this->applyBoundaryConditionsForCurrentAnalysis(A_matrix, B_matrix);
          
          A_matrix.close();
          if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
              kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
            B_matrix->close();
          
          NumericVector<double>& real_vec = this->getVector("EigenVector_Real");
          NumericVector<double>& scratch_vec = this->getVector("ScratchVec");
          FESystemDatabase::TimeIndependentEigenDataInfo eigen_vec_info_real;
          double A_factor, B_factor, eig_sens;
          // now, iterate over all the eigenvalues and calculate their sensitivity
          for (unsigned int i=0; i < eig_value_real.size(); i++)
            {
              // get the eigenvector, and calculate the sensitivities
              std::ostringstream oss;
              oss << i;
              eigen_vec_info_real.clear();
              eigen_vec_info_real.setDisciplineEnumID(it->second->getDisciplineEnumID());
              eigen_vec_info_real.setName("EigenVector_Real");
              eigen_vec_info_real.setLoadCase(this->getCurrentLoadCase());
              eigen_vec_info_real.setModeInfo(i+1);
              this->fesystem_controller.global_data_storage->fillVector(eigen_vec_info_real,real_vec);
              
              real_vec.close();
              
              // calculate the individual terms
              A_matrix.multiply_vector(real_vec, scratch_vec);
              A_factor = real_vec.dot(scratch_vec);
              if (kind_enum_ID == Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::num() ||
                  kind_enum_ID == Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::num())
                {
                  B_matrix->multiply_vector(real_vec, scratch_vec);
                  B_factor = real_vec.dot(scratch_vec);
                }
              else
                B_factor = 1.0;
              
              eig_sens = A_factor - eig_value_real(i) * B_factor;
              eig_value_real_sens(i) = eig_sens;
            }
          
          // now register the eigen vector and store it in the database
          data_info.clear();
          data_info.setDisciplineEnumID(it->second->getDisciplineEnumID());
          data_info.setName("EigenValuesReal");
          data_info.setLoadCase(this->current_load_case);
          data_info.setDVID(this->getCurrentDesignParameter().getID());
          this->fesystem_controller.global_data_storage->storeVector(data_info, 
                                                                     eig_value_real_sens);
        }
          break;
          
        case NON_HERMITIAN_EIGENPROBLEM_ENUM_ID: 
        case GENERALIZED_NON_HERMITIAN_EIGENPROBLEM_ENUM_ID: 
        default:
          Assert(false, ExcInternalError());
      }
      
    }
      break;
			
    default:
      abort();
      break;
			
  }
  
}



void 
Driver::EigenProblemAnalysisDriver::applyBoundaryConditionsForCurrentAnalysis
(SparseMatrix<double>& A_matrix,
 SparseMatrix<double>* B_matrix)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();

  
  LoadDatabase& load_database = *(this->fesystem_controller.load_database.get());
  const MeshDS::FEMeshData& mesh_data = it->second->getAnalysisMeshData();
  
  
  // get the boundary conditions
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
  load_info(it->second->getBoundaryConditionLoadInfo().release());
  
  FESystemUtility::AutoPtrVector<Loads::DirichletBoundaryConditionCombination> 
  bcond(load_database.getAllLoadCombinations
        <Loads::DirichletBoundaryConditionDataInfo, Loads::DirichletBoundaryConditionCombination>
        (*load_info).release());
  
  // to apply the boundary condition, get the global dof_number for each
  // constrained dof and use the penalty method to apply the boundary condition
  std::vector<Loads::DirichletBoundaryConditionCombination*>::const_iterator bc_it, bc_end;
  bc_it = bcond.get()->begin();
  bc_end = bcond.get()->end();
  
  const double penalty = 1.0e16;
  std::vector<unsigned int> bc_rows;
  unsigned int dof, constrained_node, constrained_dof;
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      for (; bc_it != bc_end ; bc_it++)
        {
          constrained_node =  (**bc_it).getNodeID();
          constrained_dof = (**bc_it).getDofNumber();
          dof = 
          (mesh_data.getNodeFromForeignID(constrained_node))->dof_number(0,(constrained_dof-1),0);
          
          //matrix.add(dof,dof,penalty);
          bc_rows.push_back(dof);
        } 
      
      Mat A_mat;
      A_mat = dynamic_cast<PetscMatrix<double>& >(A_matrix).mat();
      MatSetOption(A_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      
      A_matrix.perform_intermediate_matrix_assembly();
      // now, diagonalize the rows and columns for the given dofs.
      // to zero the columns, iterate over all rows and if the element for the 
      // specified dof is not zero, zero it
      for (unsigned int i=0; i < A_matrix.m(); i++)
        for (unsigned int j=0; j < bc_rows.size(); j++)
          A_matrix.set(i, bc_rows[j], 0.0);
      A_matrix.zero_rows(bc_rows, penalty);
      
      MatSetOption(A_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
      A_matrix.close();
      
      // repeat the same process for the B_matrix if it has been specified
      if (B_matrix != NULL)
        {
          Mat B_mat;
          B_mat = dynamic_cast<PetscMatrix<double>* >(B_matrix)->mat();
          B_matrix->perform_intermediate_matrix_assembly();
          MatSetOption(B_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
          // now, diagonalize the rows and columns for the given dofs.
          // the diagonal matrix here is set to a very small value, so that the
          // eigenvalue for this will be a very high number.
          // to zero the columns, iterate over all rows and if the element for the 
          // specified dof is not zero, zero it
          for (unsigned int i=0; i < B_matrix->m(); i++)
            for (unsigned int j=0; j < bc_rows.size(); j++)
              B_matrix->set(i, bc_rows[j], 0.0);
          B_matrix->zero_rows(bc_rows, penalty);
          
          MatSetOption(B_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
          B_matrix->close();
        }
    }
      break;
			
    default:
      abort();
      break;
  }
}


void 
Driver::EigenProblemAnalysisDriver::getSubmatrixAfterBoundaryCondition
(SparseMatrix<double>& A_matrix,
 SparseMatrix<double>* B_matrix,
 SparseMatrix<double>& A_matrix_sub,
 SparseMatrix<double>* B_matrix_sub)
{
  // this is a uni-discipline driver. Make sure that the number of disciplines here is one
  Assert(this->analysis_discipline_map.size() == 1, ExcInternalError());
  
  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.begin();
  
  LoadDatabase& load_database = *(this->fesystem_controller.load_database.get());
  const MeshDS::FEMeshData& mesh_data = it->second->getAnalysisMeshData();
  
  
  // get the boundary conditions
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
  load_info(it->second->getBoundaryConditionLoadInfo().release());
  
  FESystemUtility::AutoPtrVector<Loads::DirichletBoundaryConditionCombination> 
  bcond(load_database.getAllLoadCombinations
        <Loads::DirichletBoundaryConditionDataInfo, Loads::DirichletBoundaryConditionCombination>
        (*load_info).release());
  
  // to apply the boundary condition, get the global dof_number for each
  // constrained dof and use the penalty method to apply the boundary condition
  std::vector<Loads::DirichletBoundaryConditionCombination*>::const_iterator bc_it, bc_end;
  bc_it = bcond.get()->begin();
  bc_end = bcond.get()->end();
  
  std::set<unsigned int> bc_rows;
  unsigned int dof, constrained_node, constrained_dof;
  
  
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      for (; bc_it != bc_end ; bc_it++)
        {
          constrained_node =  (**bc_it).getNodeID();
          constrained_dof = (**bc_it).getDofNumber();
          dof = 
          (mesh_data.getNodeFromForeignID(constrained_node))->dof_number(0,(constrained_dof-1),0);
          
          bc_rows.insert(dof);
        }
      
      std::vector<unsigned int> non_bc_rows(A_matrix.m() - bc_rows.size());
      
      // now create the vector of dofs that are not bcd
      unsigned int n_elem = 0;
      for (unsigned int i=0; i < A_matrix.m(); i++)
        if (bc_rows.count(i) == 0)
          {
            non_bc_rows[n_elem] = i;
            n_elem++;
          }
      
      // now, get the submatrices
      A_matrix.create_submatrix(A_matrix_sub, non_bc_rows, non_bc_rows);
      B_matrix->create_submatrix(*B_matrix_sub, non_bc_rows, non_bc_rows);
      
      Mat mat1 = dynamic_cast<PetscMatrix<double>&>(A_matrix).mat();
      PetscErrorCode ierr; PetscTruth symm;
      ierr = MatIsSymmetric(mat1, 1.0e-10, &symm);
      CHKERRABORT(FESystem::COMM_WORLD, ierr);
      //std::cout << "A " << symm << std::endl;
      
      mat1 = dynamic_cast<PetscMatrix<double>*>(B_matrix)->mat();
      ierr = MatIsSymmetric(mat1, 1.0e-10, &symm);
      CHKERRABORT(FESystem::COMM_WORLD, ierr);
      //std::cout << "B " << symm << std::endl;
      
      mat1 = dynamic_cast<PetscMatrix<double>&>(A_matrix_sub).mat();
      ierr = MatIsSymmetric(mat1, 1.0e-10, &symm);
      CHKERRABORT(FESystem::COMM_WORLD, ierr);
      //std::cout << "A_sub " << symm << std::endl;
      
      mat1 = dynamic_cast<PetscMatrix<double>*>(B_matrix_sub)->mat();
      ierr = MatIsSymmetric(mat1, 1.0e-10, &symm);
      CHKERRABORT(FESystem::COMM_WORLD, ierr);
      //std::cout << "B_sub " << symm << std::endl;
    }
      break;
			
    default:
      abort();
      break;
  }
}



