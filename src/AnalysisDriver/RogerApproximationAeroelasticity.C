/*
 *  RogerApproximationAeroelasticity.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/20/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// C++ includes
#include <string>
#include <sstream>
#include <map>
#include<iostream>
#include<iomanip>
#include<fstream>

// FESystem includes
#include "AnalysisDriver/RogerApproximationAeroelasticity.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/AerodynamicDisciplineBase.h"
#include "Solutions/AeroelasticitySolution.h"
#include "Numerics/MultiBlockSparseMatrix.h"
#include "Loads/FlightCondition.h"
#include "Solvers/EigenSolver.h"
#include "Loads/LoadDatabase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"
#include "Utilities/Log.h"
#include "Utilities/AutoptrVector.h"

// libMesh includes
#include "utils/utility.h"
#include "numerics/petsc_matrix.h"

Driver::RogerApproximationAeroelasticityDriver::RogerApproximationAeroelasticityDriver
(const unsigned int ID, FESystem::FESystemController& fesys_controller):
Driver::AeroelasticityAnalysisDriverBase(ID, fesys_controller, Driver::ROGER_APPROXIMATION_AEROELASTICITY_DRIVER::num())
{  
}


Driver::RogerApproximationAeroelasticityDriver::~RogerApproximationAeroelasticityDriver()
{
	
}



void
Driver::RogerApproximationAeroelasticityDriver::initialize()
{
  this->addMatricesAndVectors();
}



void
Driver::RogerApproximationAeroelasticityDriver::addMatricesAndVectors()
{
  Discipline::StructuralAnalysis& structure = 
  dynamic_cast<Discipline::StructuralAnalysis&>(this->getStructuralDiscipline());
  Discipline::AerodynamicDisciplineBase& aero = 
  dynamic_cast<Discipline::AerodynamicDisciplineBase&>(this->getAerodynamicDiscipline());
  
  
  // add a system matrix, a solution vector and a rhs
  if (structure.disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
    structure.initMatrix(this->addMatrix("StructuralStiffness"));
  if (structure.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
    structure.initMatrix(this->addMatrix("StructuralDamping"));
  if (structure.disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
    structure.initMatrix(this->addMatrix("StructuralMass"));
  
  
  // check with the aerodynamic disciplines and add the right ones
  // the dof map of the matrix is used from the structures, since 
  // aero discipline does not add any
  if (aero.disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
    structure.initMatrix(this->addMatrix("AerodynamicStiffness"));
  if (aero.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
    structure.initMatrix(this->addMatrix("AerodynamicDamping"));
  if (aero.disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
    structure.initMatrix(this->addMatrix("AerodynamicMass"));
  
  // get a reference to the solution 
  const Solution::AeroelasticitySolution& aeroelasticity_sol = 
  dynamic_cast<const Solution::AeroelasticitySolution&> (aero.getSolution());
  
  for (unsigned int i=0; i < aeroelasticity_sol.numberOfAerodynamicLagTerms(); i++)
    {
      std::ostringstream lag_num;
      lag_num << i;
      std::string name = "AerodynamicLagMatrix";
      name += lag_num.str();
      aero.initMatrix(this->addMatrix(name));
    }
}





void
Driver::RogerApproximationAeroelasticityDriver::solveCurrentLoadCase()
{    
  
  Discipline::StructuralAnalysis& structure = 
  dynamic_cast<Discipline::StructuralAnalysis&>(this->getStructuralDiscipline());
  Discipline::AerodynamicDisciplineBase& aero = 
  dynamic_cast<Discipline::AerodynamicDisciplineBase&>(this->getAerodynamicDiscipline());
  const Solution::AeroelasticitySolution& aeroelasticity_sol = 
  dynamic_cast<const Solution::AeroelasticitySolution&> (aero.getSolution());
    
  std::map<unsigned int, double> real_map;
  std::map<unsigned int, SparseMatrix<double>*> aero_matrix_map;
  std::map<unsigned int, SparseMatrix<double>*> structural_matrix_map;
  std::map<unsigned int, NumericVector<double>*> vector_map;
  
  // prepare the load vector for the analysis and then init the solver and then finally 
  // solve the system
  switch (this->current_analysis_kind_enum_ID)
  {
    case ANALYSIS_ENUM_ID:
    {
      // get the structural matrices
      if (structure.disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
        structural_matrix_map[Driver::SYSTEM_MATRIX::num()] = &this->getMatrix("StructuralStiffness");
      if (structure.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
        structural_matrix_map[Driver::TRANSIENT_C1_MATRIX::num()] = &this->getMatrix("StructuralDamping");
      if (structure.disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
        structural_matrix_map[Driver::TRANSIENT_C2_MATRIX::num()] = &this->getMatrix("StructuralMass");
      
      structure.fillQuantity(real_map, structural_matrix_map, vector_map);
      
      // get the aerodynamic matrices
      if (aero.disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
        aero_matrix_map[Driver::SYSTEM_MATRIX::num()] = &this->getMatrix("AerodynamicStiffness");
      if (aero.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
        aero_matrix_map[Driver::TRANSIENT_C1_MATRIX::num()] = &this->getMatrix("AerodynamicDamping");
      if (aero.disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
        aero_matrix_map[Driver::TRANSIENT_C2_MATRIX::num()] = &this->getMatrix("AerodynamicMass");
      
      // create the matrices for computation of eigenvalues
      std::auto_ptr<SparseMatrix<double> > mass_mat(new PetscMatrix<double> ()),
      damping_mat(new PetscMatrix<double> ()), stiffness_mat(new PetscMatrix<double>()), 
      mass_mat_bc(new PetscMatrix<double> ()), damping_mat_bc(new PetscMatrix<double> ()), 
      stiffness_mat_bc(new PetscMatrix<double>()), identity_mat(new PetscMatrix<double> ()),
      zero_matrix(new PetscMatrix<double>()); 
//      roger_identity_matrix(new PetscMatrix<double>());
      
      //std::auto_ptr<NumericVector<double> > mass_diagonalized(new PetscVector<double>());
      
      // now iterate over the dynamic pressures for which solution has to be obtained
      unsigned int n_lags_terms = aeroelasticity_sol.numberOfAerodynamicLagTerms();
      // get the current load case, and use the aeroelasticity load set to get the dynamic pressures
      const std::vector<Loads::FlightCondition>& flt_conds = 
      aeroelasticity_sol.getFlightConditionVector();
      
      // now, create the combined Roger Approximation matrix
      std::auto_ptr<FESystemNumerics::MultiBlockSparseMatrix<double> > 
      roger_matrix(new FESystemNumerics::PetscMultiBlockSparseMatrix<double>()),
      roger_lhs_matrix(new FESystemNumerics::PetscMultiBlockSparseMatrix<double>());
      bool updated_sparsity = false;
      double dyn_press = 0.0;
      
      std::vector<SparseMatrix<double>*> matrices;
      std::vector<double> eig_val_vec;

      std::vector<int> local_to_global_map;
      IS row_matrix_index_set, col_matrix_index_set;
      this->createSubMatrixIndices(structure, &row_matrix_index_set,
                                   &col_matrix_index_set, local_to_global_map);
      
      
      for (unsigned int ii=0; ii<flt_conds.size(); ii++)
        {
          //set the current flight condition number and ask for the matrix
          this->current_flight_condition = &(flt_conds[ii]);

          // zero the matrices in the aero matrix map
          std::map<unsigned int, SparseMatrix<double>*>::iterator it, end;
          it = aero_matrix_map.begin();
          end = aero_matrix_map.end();
          for ( ; it != end; it++)
            {
              it->second->zero();
              it->second->close();
            }
          
          aero.fillQuantity(real_map, aero_matrix_map, vector_map);
          
          this->current_flight_condition = NULL;
          dyn_press = flt_conds[ii].getDynamicPressure();
          
          // now prepare the matrices for solution
          // first create the composite mass, damping and stiffness matrices; 
          // mass matrix: M = M_str - q_D * M_A
          if (!mass_mat->initialized())
            mass_mat->duplicate_matrix(this->getMatrix("StructuralMass") , true);
          else 
            {
              mass_mat->zero();
              mass_mat->close();
              mass_mat->copy_matrix(this->getMatrix("StructuralMass"));
            }
          mass_mat->close();
          if (aero.disciplineHasMatrix(Driver::TRANSIENT_C2_MATRIX::num()))
            mass_mat->add(-dyn_press, this->getMatrix("AerodynamicMass"));
          mass_mat->close();
          
          this->getSubmatrixAfterBoundaryCondition(*mass_mat, *mass_mat_bc,
                                                   row_matrix_index_set,
                                                   col_matrix_index_set);
          if (!zero_matrix->initialized())
            {
              zero_matrix->duplicate_matrix(*mass_mat_bc, false);
              zero_matrix->zero();
              zero_matrix->close();
            }
          
          
//          // the mass matrix is diagonalized and then inverted to save CPU time. 
//          if (!mass_diagonalized->initialized())
//          {
//            mass_diagonalized->init(local_to_global_map.size());
//            mass_diagonalized->zero();
//            mass_diagonalized->close(); 
//          }
//          
//          mass_mat_bc->diagonalize(*mass_diagonalized);
//          mass_diagonalized->invert_values();
//          mass_diagonalized->close(); 
          
          // damping matrix: C = C_str - q_D * C_A 
          if (structure.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
            {
              if (!damping_mat->initialized())
                damping_mat->duplicate_matrix(this->getMatrix("StructuralDamping") , true); 
              else 
                {
                  damping_mat->zero();
                  damping_mat->close();
                  damping_mat->copy_matrix(this->getMatrix("StructuralDamping")); 
                }
              damping_mat->close();

              if (aero.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
                damping_mat->add(-dyn_press, this->getMatrix("AerodynamicDamping"));
              damping_mat->close();
            }
          else if (aero.disciplineHasMatrix(Driver::TRANSIENT_C1_MATRIX::num()))
            {
              if (!damping_mat->initialized())
                damping_mat->duplicate_matrix(this->getMatrix("AerodynamicDamping"), true);
              else
                {
                  damping_mat->zero();
                  damping_mat->close();
                  damping_mat->copy_matrix(this->getMatrix("AerodynamicDamping"));
                }
              
              damping_mat->close();
              damping_mat->scale(-dyn_press);
              damping_mat->close();
            }          


          this->getSubmatrixAfterBoundaryCondition(*damping_mat, *damping_mat_bc, 
                                                   row_matrix_index_set,
                                                   col_matrix_index_set);

          //damping_mat_bc->leftDiagonalScale(*mass_diagonalized);
          damping_mat_bc->close();
          damping_mat_bc->scale(-1.0);
          
          // stiffness matrix: K = K_str - q_D * K_A
          if (!stiffness_mat->initialized())
            stiffness_mat->duplicate_matrix(this->getMatrix("StructuralStiffness") , true);
          else
            {
              stiffness_mat->zero();
              stiffness_mat->close();
              stiffness_mat->copy_matrix(this->getMatrix("StructuralStiffness"));
            }
          
          stiffness_mat->close();
          if (aero.disciplineHasMatrix(Driver::SYSTEM_MATRIX::num()))
            stiffness_mat->add(-dyn_press, this->getMatrix("AerodynamicStiffness"));
          stiffness_mat->close();

         
          this->getSubmatrixAfterBoundaryCondition(*stiffness_mat, *stiffness_mat_bc,
                                                   row_matrix_index_set,
                                                   col_matrix_index_set);

          //stiffness_mat_bc->leftDiagonalScale(*mass_diagonalized);
          stiffness_mat_bc->close();
          stiffness_mat_bc->scale(-1.0);
                    
          if (!identity_mat->initialized())
            {
              identity_mat->duplicate_matrix(*stiffness_mat_bc, false);
              identity_mat->zero();
              identity_mat->shift(1.0);
              identity_mat->close();
            }

          // if the sparsity has not been defined for the matrix, do so
          if (!updated_sparsity)
            {
              // set the dimensions
              roger_matrix->setMatrixDimension(2+n_lags_terms, 2+n_lags_terms);
              // define the matrix
              roger_matrix->setSubMatrix(0,0, *zero_matrix);
              roger_matrix->setSubMatrix(0,1, *identity_mat);
              roger_matrix->setSubMatrix(1,0, *stiffness_mat_bc);
              roger_matrix->setSubMatrix(1,1, *damping_mat_bc);

              // the lhs matrix
              roger_lhs_matrix->setMatrixDimension(2+n_lags_terms, 2+n_lags_terms);
              roger_lhs_matrix->setSubMatrix(0,0, *identity_mat);
              roger_lhs_matrix->setSubMatrix(0,1, *zero_matrix);
              roger_lhs_matrix->setSubMatrix(1,0, *zero_matrix);
              roger_lhs_matrix->setSubMatrix(1,1, *mass_mat_bc);
              
              updated_sparsity = true;
            }
          
          // now, update the values of the matrix before finally solving the system
          roger_matrix->update();
          roger_lhs_matrix->update();
          SparseMatrix<double>& rm_rhs_mat = roger_matrix->getSparseMatrixObject();
          SparseMatrix<double>& rm_lhs_mat = roger_lhs_matrix->getSparseMatrixObject();
          
//          if (!roger_identity_matrix->initialized())
//            {
//              roger_identity_matrix->duplicate_matrix(rm, false);
//              roger_identity_matrix->zero();
//              roger_identity_matrix->shift(1.0);
//              roger_identity_matrix->close();
//            }
//          rm.print_matlab("mat.m");
          
          // now solve the eigenvalue problem
          Solver::EigenSolverBase* eigen_solver =dynamic_cast<Solver::EigenSolverBase*>(this->solver);
          eigen_solver->solve(&rm_rhs_mat, &rm_lhs_mat);
          

          // after the solution, get the eigenvectors and store them one by one
          {
            FESystemDatabase::TimeIndependentEigenDataInfo eigen_vec_info_real, eigen_vec_info_img;
            double value_r, value_i;
            std::string load_case_string = "Load Case : ";
            
            std::ostringstream lc_oss;
            lc_oss << this->getCurrentLoadCase();
            load_case_string += lc_oss.str();
            
            unsigned int n_converged = eigen_solver->getNConvergedEigenPairs();
            DenseVector<double> eig_value_real(n_converged), eig_value_img(n_converged);
            // prepare the vectors to store the eigen values
            std::auto_ptr<NumericVector<double> > 
            real_vec(NumericVector<double>::build().release()),
            img_vec(NumericVector<double>::build().release()),
            str_real_vec(NumericVector<double>::build().release()),
            str_img_vec(NumericVector<double>::build().release());
            
            // ask the roger matrix to initialize these vectors
            roger_matrix->initRightVector(*real_vec);
            roger_matrix->initRightVector(*img_vec);
            structure.initVector(*str_real_vec);
            structure.initVector(*str_img_vec);
            std::auto_ptr<NumericVector<double> > tmp_vec(new PetscVector<double>());
            tmp_vec->init(rm_rhs_mat.m() / 2);
            
            // create the data structure to store the eigenvectors and get them from the solver
            FESystemUtility::AutoPtrVector<NumericVector<double> > numeric_vectors(n_converged);
            for (unsigned int i=0; i<n_converged; i++)
            {
              numeric_vectors.reset(i, new PetscVector<double>());
              numeric_vectors[i]->init(rm_rhs_mat.m());
            }
            
            eigen_solver->getInvariantSubspace(numeric_vectors.getReference());
            unsigned int val_index = 0;
            // also, save the eigen vectors
            FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, n_converged);
            for (unsigned int i=0; i < n_converged; i++)
              {
                std::ostringstream oss;
                oss << i;
                std::ostringstream oss_eigvec_name_real;
                std::ostringstream oss_eigvec_name_img;
                
                oss_eigvec_name_real << "RogerEigenVector_Real_FltCond_" << ii;
                oss_eigvec_name_img << "RogerEigenVector_Img_FltCond_" << ii;
                
                eigen_vec_info_real.clear();
                eigen_vec_info_real.setDisciplineEnumID(structure.getDisciplineEnumID());
                eigen_vec_info_real.setName(oss_eigvec_name_real.str());
                eigen_vec_info_real.setLoadCase(this->getCurrentLoadCase());
                eigen_vec_info_real.setModeInfo(i+1);
                
                eigen_vec_info_img.clear();
                eigen_vec_info_img.setDisciplineEnumID(structure.getDisciplineEnumID());
                eigen_vec_info_img.setName(oss_eigvec_name_img.str());
                eigen_vec_info_img.setLoadCase(this->getCurrentLoadCase());
                eigen_vec_info_img.setModeInfo(i+1);
                
                // get the eigen pair
                eigen_solver->getEigenPair(i, &value_r, &value_i, *real_vec, *img_vec);

                // append the eigenvalue and the dynamic pressure to the eig_val_vec 
                eig_val_vec.push_back(value_r);
                eig_val_vec.push_back(value_i);
                eig_val_vec.push_back(dyn_press);
                
                std::string value;
                std::ostringstream oss_value;
                value = "EigenValue Number  ";
                value += oss.str();
                value += "   =   ";
                oss_value.setf(std::ios::showpoint);
                oss_value << std::setprecision(14) << value_r;
                oss_value << "  + i  ";
                oss_value << std::setprecision(14) << value_i;
                value += oss_value.str();
                
                FESystemUtility::Parallel::PrintAll(FESystem::COMM_WORLD, value);
                
                // copy the eigenvalue for storage
                eig_value_real(i) = value_r;
                eig_value_img(i) = value_i;
                                                
                // store the multiblock vector
                this->fesystem_controller.global_data_storage->storeVector
                (eigen_vec_info_real,*real_vec);
                this->fesystem_controller.global_data_storage->storeVector
                (eigen_vec_info_img, *img_vec);

                // now extract the structural eigenvector and write it to the storage 
                std::ostringstream oss_str_eigvec_name_real;
                std::ostringstream oss_str_eigvec_name_img;
                
                oss_str_eigvec_name_real << "EigenVector_Real_FltCond_" << ii;
                oss_str_eigvec_name_img << "EigenVector_Img_FltCond_" << ii;
                
                eigen_vec_info_real.clear();
                eigen_vec_info_real.setDisciplineEnumID(structure.getDisciplineEnumID());
                eigen_vec_info_real.setName(oss_str_eigvec_name_real.str());
                eigen_vec_info_real.setLoadCase(this->getCurrentLoadCase());
                eigen_vec_info_real.setModeInfo(i+1);
                
                eigen_vec_info_img.clear();
                eigen_vec_info_img.setDisciplineEnumID(structure.getDisciplineEnumID());
                eigen_vec_info_img.setName(oss_str_eigvec_name_img.str());
                eigen_vec_info_img.setLoadCase(this->getCurrentLoadCase());
                eigen_vec_info_img.setModeInfo(i+1);
                
                // extract the vector from the multiblock vector and 
                // copy the subvector to the global vector
                val_index = 0;
                if (i % 2 == 0)
                  val_index = i;
                else
                  val_index = i-1;
                  
                tmp_vec->zero();
                roger_matrix->extractFromRightMultiblockVector(0, *(numeric_vectors[val_index]), *tmp_vec); 
                this->copyVectorToGlobalVector(*tmp_vec , *str_real_vec, local_to_global_map);
                
                
                // same for the imaginary vector
                tmp_vec->zero();
                roger_matrix->extractFromRightMultiblockVector(0, *(numeric_vectors[val_index+1]), *tmp_vec); 
                this->copyVectorToGlobalVector(*tmp_vec , *str_img_vec, local_to_global_map); 
                
                
                
                // now write the solutions to the log file
                if (FESystem::local_processor == 0)
                  {
                    this->fesystem_controller.log_file->write(load_case_string);
                    this->fesystem_controller.log_file->write("EigenSolution");
                    this->fesystem_controller.log_file->write(value);
                    this->fesystem_controller.log_file->write("Eigen Vector: ");
                  }
                
                this->writeComplexSolutionVectorToLog
                (eigen_vec_info_real, *str_real_vec, *str_img_vec, structure.getDisciplineEnumID());
                
                this->fesystem_controller.global_data_storage->storeVector(eigen_vec_info_real,
                                                                           *str_real_vec);
                this->fesystem_controller.global_data_storage->storeVector(eigen_vec_info_img,
                                                                           *str_img_vec);
              }
            
            std::ostringstream oss_eigval_name_real;
            std::ostringstream oss_eigval_name_img;
            oss_eigval_name_real << "RogerEigenValue_Real_FltCond_" << ii;
            oss_eigval_name_img << "RogerEigenValue_Img_FltCond_" << ii;

            // now register the eigen vector and store it in the database
            FESystemDatabase::TimeIndependentDataInfo data_info;
            data_info.clear();
            data_info.setDisciplineEnumID(structure.getDisciplineEnumID());
            data_info.setName(oss_eigval_name_real.str());
            data_info.setLoadCase(this->current_load_case);
            this->fesystem_controller.global_data_storage->storeVector(data_info, eig_value_real);
            
            data_info.clear();
            data_info.setDisciplineEnumID(structure.getDisciplineEnumID());
            data_info.setName(oss_eigval_name_img.str());
            data_info.setLoadCase(this->current_load_case);
            this->fesystem_controller.global_data_storage->storeVector(data_info, eig_value_img);
          }
        }
      
      // write the eigenvalues to an output file
      std::ostringstream file_name;
      file_name << "eig_vals_lc_";
      file_name << this->getCurrentLoadCase();
      file_name << ".txt";
      std::fstream output_file;
      output_file.open(file_name.str().c_str(), std::fstream::out);
      
      for (unsigned int ii=0; ii<(eig_val_vec.size()/3); ii++)
        {
          for (unsigned int i=0; i<3; i++)
            output_file << std::showpoint << std::setprecision(7) 
            << eig_val_vec[ii*3+i] << "    ";
          output_file << std::endl;
        }
      output_file.close();
      
      // destroy the index set
      this->destroySubMatrixIndices(row_matrix_index_set, col_matrix_index_set);
    }
      break;
      
      
    case SENSITIVITY_ANALYSIS_ENUM_ID:
    {
      Assert(false, ExcInternalError());
    }
      break;
			
    default:
      abort();
      break;
			
  }
}




void 
Driver::RogerApproximationAeroelasticityDriver::getSubmatrixAfterBoundaryCondition
(SparseMatrix<double>& matrix,
 SparseMatrix<double>& matrix_sub,
 IS row_matrix_index_set, IS col_matrix_index_set)
{  
  Assert(matrix.initialized(), ExcInternalError());
  
  if (matrix_sub.initialized()) 
    matrix_sub.clear(); 
  
  Mat global_mat = dynamic_cast<PetscMatrix<double>&> (matrix).mat(); 
  Mat* local_mat = (dynamic_cast<PetscMatrix<double>&> (matrix_sub).mat_ptr()); 
  
  PetscErrorCode ierr = 0; 
  
  // create the submatrix 
  ierr = MatGetSubMatrix(global_mat, row_matrix_index_set, 
                         col_matrix_index_set, PETSC_DECIDE, MAT_INITIAL_MATRIX, local_mat);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
  matrix_sub.initialized_matrix_externally();
  matrix_sub.close();
}





void 
Driver::RogerApproximationAeroelasticityDriver::copyVectorToGlobalVector
(NumericVector<double>& sub_vector,
 NumericVector<double>& global_vector,
 const std::vector<int>& local_to_global_map)
{  
  
  Assert(local_to_global_map.size() > 0, ExcInternalError());
  Assert(global_vector.size() > 0, ExcInternalError());
  Assert(local_to_global_map.size() == sub_vector.size(), ExcInternalError());
  
  sub_vector.close();
  global_vector.close();
  global_vector.zero();
  
  int ndofs = local_to_global_map.size();
  
  Vec global_vec = dynamic_cast<PetscVector<double>&> (global_vector).vec();
  Vec local_vec = dynamic_cast<PetscVector<double>&> (sub_vector).vec();
  
  std::vector<int> local_ids(ndofs);
  Utility::iota(local_ids.begin(), local_ids.end(), 0);

  std::vector<double> vals(ndofs);
  std::fill(vals.begin(), vals.end(), 0.0);
  
  PetscErrorCode ierr = 0;
  
  // create the submatrix 
  ierr = VecGetValues(local_vec, ndofs, &local_ids[0], &vals[0]);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
  ierr = VecSetValues(global_vec, ndofs, &local_to_global_map[0], &vals[0], INSERT_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);

  global_vector.close();
}






void 
Driver::RogerApproximationAeroelasticityDriver::createSubMatrixIndices
(Discipline::AnalysisDisciplineBase& discipline,
 IS* row_matrix_index_set, IS* col_matrix_index_set,
 std::vector<int>& local_to_global_map)
{  
  LoadDatabase& load_database = *(this->fesystem_controller.load_database.get());
  const MeshDS::FEMeshData& mesh_data = discipline.getAnalysisMeshData();
  
  
  // get the boundary conditions
  std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
  load_info(discipline.getBoundaryConditionLoadInfo().release());
  
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
  
  
  for (; bc_it != bc_end ; bc_it++)
    {
      constrained_node =  (**bc_it).getNodeID();
      constrained_dof = (**bc_it).getDofNumber();
      dof = 
      (mesh_data.getNodeFromForeignID(constrained_node))->dof_number(0,(constrained_dof-1),0);
      
      bc_rows.insert(dof);
    }
  
  std::vector<unsigned int> non_bc_rows(discipline.getNDofs() - bc_rows.size());
  
  // now create the vector of dofs that are not bcd
  unsigned int n_elem = 0;
  for (unsigned int i=0; i < discipline.getNDofs(); i++)
    if (bc_rows.count(i) == 0)
      {
        non_bc_rows[n_elem] = i;
        n_elem++;
      }
  
  
  int ierr=0;
  ISLocalToGlobalMapping mapping;
  
  // create the index set for the unconstrained rows
  ierr = ISCreateGeneral(FESystem::COMM_WORLD,
                         non_bc_rows.size(),
                         (int*) &non_bc_rows[0],
                         row_matrix_index_set); 
  CHKERRABORT(FESystem::COMM_WORLD,ierr);

  // create the index set for the unconstrained rows
  ierr = ISCreateGeneral(FESystem::COMM_WORLD,
                         non_bc_rows.size(),
                         (int*) &non_bc_rows[0],
                         col_matrix_index_set); 
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
  // create a mapping context for the unconstrained rows
  ierr = ISLocalToGlobalMappingCreateIS(*col_matrix_index_set, &mapping);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
  
  // the mapped indices will be stored in these vectors
  std::vector<int> local_ids(non_bc_rows.size());
  Utility::iota(local_ids.begin(), local_ids.end(), 0);

  local_to_global_map.resize(non_bc_rows.size());
  std::fill(local_to_global_map.begin(), local_to_global_map.end(), 0);
  
  // get the vector of mapped local to global ids 
  ierr = ISLocalToGlobalMappingApply(mapping, 
                                     (PetscInt) local_to_global_map.size(), 
                                     (PetscInt*) &local_ids[0], 
                                     (PetscInt*) &local_to_global_map[0]);  
  CHKERRABORT(FESystem::COMM_WORLD,ierr);

  // destroy this
  ierr = ISLocalToGlobalMappingDestroy(mapping);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);

}





void 
Driver::RogerApproximationAeroelasticityDriver::destroySubMatrixIndices
(IS row_matrix_index_set, IS col_matrix_index_set)
{
  PetscErrorCode ierr = 0;
  
  ierr = ISDestroy(row_matrix_index_set);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);

  ierr = ISDestroy(col_matrix_index_set);
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
}



