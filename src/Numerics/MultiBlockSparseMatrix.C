// $Id:$
/*
 *  MultiBlockSparseMatrix.C
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/28/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */


// C++ includes


// FESystem includes
#include "Numerics/MultiBlockSparseMatrix.h"
#include "FESystem/FESystemExceptions.h"
#include "FESystem/FESystemController.h"

// libmesh includes
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"
#include "numerics/petsc_matrix.h"
#include "numerics/petsc_vector.h"
#include "base/dof_map.h"
#include "utils/utility.h"

template<typename T>
FESystemNumerics::MultiBlockSparseMatrix<T>::MultiBlockSparseMatrix():
finalized_sparsity(false),
row_dimension(0),
column_dimension(0),
matrix(NULL)
{
  
}




template<typename T>
FESystemNumerics::MultiBlockSparseMatrix<T>::~MultiBlockSparseMatrix()
{
  this->clear();
}



template<typename T>
void
FESystemNumerics::MultiBlockSparseMatrix<T>::clear()
{
  this->finalized_sparsity = false;
  this->row_dimension = 0;
  this->column_dimension = 0;
  this->matrix.reset();
  this->sub_matrices.clear();
  this->rows_per_row_wise_block.clear();
  this->columns_per_column_wise_block.clear();
}




template<typename T>
void
FESystemNumerics::MultiBlockSparseMatrix<T>::setMatrixDimension
(const unsigned int row_size, const unsigned int column_size)
{
  // first make sure that the current dimensions of the system are zero
  Assert(this->finalized_sparsity == false, ExcInternalError());
  Assert(this->sub_matrices.size() == 0, ExcInternalError());
  // next, check for the sanity of values
  Assert(row_size > 0, ExcInternalError());
  Assert(column_size > 0, ExcInternalError());

  
  this->row_dimension = row_size;
  this->column_dimension = column_size;
  // now create space in the system to store the matrix pointers
  this->sub_matrices.resize(row_size * column_size);
}



template<typename T>
void
FESystemNumerics::MultiBlockSparseMatrix<T>::setSubMatrix
(const unsigned int row_num, const unsigned int column_num, SparseMatrix<T>& matrix)
{
  // make sure that this object dimensions have been set and that the given 
  // dimensions agree with them
  Assert(this->finalized_sparsity == false, ExcInternalError());
  Assert(this->sub_matrices.size() > 0, ExcInternalError());
  Assert(row_num < this->row_dimension, ExcInternalError());
  Assert(column_num < this->column_dimension, ExcInternalError());
  
  this->sub_matrices[row_num * this->column_dimension + column_num] = &matrix;
}




template<typename T>
void
FESystemNumerics::MultiBlockSparseMatrix<T>::setSubMatrixToNull
(const unsigned int row_num, const unsigned int column_num)
{
  // make sure that this object dimensions have been set and that the given 
  // dimensions agree with them
  Assert(this->finalized_sparsity == false, ExcInternalError());
  Assert(this->sub_matrices.size() > 0, ExcInternalError());
  Assert(row_num < this->row_dimension, ExcInternalError());
  Assert(column_num < this->column_dimension, ExcInternalError());
    
  this->sub_matrices[row_num * this->column_dimension + column_num] = NULL;
}





template<typename T>
SparseMatrix<T>& 
FESystemNumerics::MultiBlockSparseMatrix<T>::getSparseMatrixObject()
{
  // return a reference to the sparse matrix object
  Assert(this->finalized_sparsity == true, ExcInternalError());
  Assert(this->matrix.get() != NULL, ExcInternalError());
  
  return *(this->matrix.get());
}




template<typename T>
FESystemNumerics::PetscMultiBlockSparseMatrix<T>::PetscMultiBlockSparseMatrix():
MultiBlockSparseMatrix<T>()
{
  
}




template<typename T>
FESystemNumerics::PetscMultiBlockSparseMatrix<T>::~PetscMultiBlockSparseMatrix()
{
  
}



template<typename T>
void
FESystemNumerics::PetscMultiBlockSparseMatrix<T>::finalizeSparsity()
{
  // make sure that the sparsity has not been finalized. 
  Assert(this->finalized_sparsity == false, ExcInternalError());
  
  // get the dimensions of the finale matrix
  this->rows_per_row_wise_block.resize(this->row_dimension);
  this->columns_per_column_wise_block.resize(this->column_dimension);
  unsigned int total_n_rows = 0, total_n_columns = 0, elem_num = 0;
  
  // fill the vector with zeros
  for (unsigned int i=0; i < this->row_dimension; i++)
    this->rows_per_row_wise_block[i] = 0;
  for (unsigned int i=0; i < this->column_dimension; i++)
    this->columns_per_column_wise_block[i] = 0;
  
  // iterate over all the submatrices and close them
  for (unsigned int i=0; i<this->sub_matrices.size(); i++)
    if (this->sub_matrices[i] != NULL)
      this->sub_matrices[i]->close();
  
  // in the process, also make sure that the number of rows and columns of the matrices 
  // consistent
  for (unsigned int i=0; i < this->row_dimension; i++)
    for (unsigned int j=0; j < this->column_dimension; j++)
      {
        elem_num = i * this->column_dimension + j;
        SparseMatrix<T>* mat = this->sub_matrices[elem_num];
        
        // if the size is zero, check if you can set it. If it is non-zero, then 
        // make sure the size is consistent
        
        // check for row size
        if (this->rows_per_row_wise_block[i] == 0)
          {if (mat != NULL) this->rows_per_row_wise_block[i] = mat->m();}
        else
          {if (mat != NULL) 
            Assert(this->rows_per_row_wise_block[i] == mat->m(), ExcInternalError());}

        // check for column size
        if (this->columns_per_column_wise_block[j] == 0)
          {if (mat != NULL) this->columns_per_column_wise_block[j] = mat->n();}
        else
          {if (mat != NULL)
            Assert(this->columns_per_column_wise_block[i] == mat->n(), ExcInternalError());}
      }
  
  // calculate the dimensions
  for (unsigned int i=0; i < this->row_dimension; i++)
    total_n_rows += this->rows_per_row_wise_block[i];
  for (unsigned int i=0; i < this->column_dimension; i++)
    total_n_columns += this->columns_per_column_wise_block[i];
  
  // next calculate the number of non-zero columns per row
  std::vector<unsigned int> nonzeros_per_row(total_n_rows);
  std::fill (nonzeros_per_row.begin(), nonzeros_per_row.end(), 0);

  int starting_row_id = 0, starting_col_id=0, row_num = 0,
  n_rows =0, n_cols = 0, next_elem_in_id_vec = 0;
  PetscErrorCode ierr = 0;

  // next, set the number of nonzero entries in the matrix 
  for (unsigned int i=0; i<this->row_dimension; i++)
    {
      for (unsigned int k=0; k<this->rows_per_row_wise_block[i]; k++)
        {
          // iterate over each matrix for this row and get the values to fill 
          // into the multiblock matrix
          for (unsigned int j=0; j<this->column_dimension; j++)
            {
              elem_num = i * this->column_dimension + j;
              SparseMatrix<T>* mat = this->sub_matrices[elem_num];
              if (mat != NULL)
                {
                  Mat pmat = dynamic_cast<PetscMatrix<T>*> (mat)->mat();
                  // get the row ij indices
                  ierr = MatGetSize(pmat, &n_rows, &n_cols);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                  Assert(n_rows == this->rows_per_row_wise_block[i], ExcInternalError());
                  Assert(n_cols == this->columns_per_column_wise_block[j], ExcInternalError());
                  
                  const PetscInt *petsc_col_ids = PETSC_NULL;
                  const PetscScalar *petsc_vals = PETSC_NULL;
                  // iterate over the rows in this matrix
                  int n_cols_in_row = 0;
                  ierr = MatGetRow(pmat, k, &n_cols_in_row, &petsc_col_ids, &petsc_vals);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                  
                  nonzeros_per_row[starting_row_id + k] += n_cols_in_row;
                  
                  // now restore the matrix row
                  ierr = MatRestoreRow(pmat, k, &n_cols_in_row, &petsc_col_ids, &petsc_vals);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                }
              else
                {
                  // add the diagonal entry, since that will be necessary for the 
                  // linear/eigen solvers 
                  nonzeros_per_row[starting_row_id + k] += 1; 
                }
            }
        }
      starting_row_id += this->rows_per_row_wise_block[i];
    }
  
  // now that the sparsity pattern is available, set it for the matrix 
  ierr = MatCreateSeqAIJ(FESystem::COMM_WORLD, total_n_rows, total_n_columns,
                         PETSC_DEFAULT, (PetscInt*)&(nonzeros_per_row[0]), &(this->petsc_mat));
  CHKERRABORT(FESystem::COMM_WORLD,ierr);
  
  // now create the sparse matrix object
  this->matrix.reset(new PetscMatrix<double>(this->petsc_mat));
  //this->matrix->close();

  //
  // now that the matrix is created, we will update the sparsity pattern of the matrix
  //
  
  // this column is used to fill the values for each row. The size is kept constant
  // and all values are kept as 1.0. Once the sparsity is set, the matrix anyway will 
  // be zeroes out. 
  std::vector<double> vals(total_n_columns);
  std::fill(vals.begin(), vals.end(), 1.0);
  // this is the vector of column ids that will be used to store the matrix 
  std::vector<int> col_ids(total_n_columns);
  std::fill(col_ids.begin(), col_ids.end(), -1.0);
  
  starting_row_id = 0, starting_col_id=0, row_num = 0,
  n_rows =0, n_cols = 0, next_elem_in_id_vec = 0;
  
  // next, set the values of the matrix entries
  for (unsigned int i=0; i<this->row_dimension; i++)
    {
      for (unsigned int k=0; k<this->rows_per_row_wise_block[i]; k++)
        {
          std::fill(col_ids.begin(), col_ids.end(), -1.0);
          // iterate over each matrix for this row and get the values to fill 
          // into the multiblock matrix
          starting_col_id = 0;
          next_elem_in_id_vec = 0;
          for (unsigned int j=0; j<this->column_dimension; j++)
            {
              elem_num = i * this->column_dimension + j;
              SparseMatrix<T>* mat = this->sub_matrices[elem_num];
              if (mat != NULL)
                {
                  Mat pmat = dynamic_cast<PetscMatrix<T>*> (mat)->mat();
                  // get the row ij indices
                  ierr = MatGetSize(pmat, &n_rows, &n_cols);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                  Assert(n_rows == this->rows_per_row_wise_block[i], ExcInternalError());
                  Assert(n_cols == this->columns_per_column_wise_block[j], ExcInternalError());
                  
                  const PetscInt *petsc_col_ids = PETSC_NULL;
                  const PetscScalar *petsc_vals = PETSC_NULL;
                  // iterate over the rows in this matrix
                  int n_cols_in_row = 0;
                  ierr = MatGetRow(pmat, k, &n_cols_in_row, &petsc_col_ids, &petsc_vals);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                  for (unsigned int l=0; l<n_cols_in_row; l++)
                    col_ids[next_elem_in_id_vec+l] = starting_col_id + petsc_col_ids[l];
                  next_elem_in_id_vec += n_cols_in_row;
                  
                  // now restore the matrix row
                  ierr = MatRestoreRow(pmat, k, &n_cols_in_row, &petsc_col_ids, &petsc_vals);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                }
              else
                {
                  // add the diagonal entry, since that will be necessary for the 
                  // linear/eigen solvers 
                  col_ids[next_elem_in_id_vec] = starting_col_id + k; 
                  next_elem_in_id_vec += 1;
                }
              
              starting_col_id += this->columns_per_column_wise_block[j];
            }

          row_num = starting_row_id + k;
          // now tell the matrix to set the values for this row
          ierr = MatSetValues (this->petsc_mat,
                               1, (PetscInt*) &(row_num),
                               (PetscInt) (nonzeros_per_row[starting_row_id + k]),
                               (PetscInt*) &(col_ids[0]),
                               (PetscScalar*) &(vals[0]), INSERT_VALUES);
          CHKERRABORT(libMesh::COMM_WORLD,ierr);
          
        }
      starting_row_id += this->rows_per_row_wise_block[i];
    }
  
  ierr = MatSetOption (this->petsc_mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
  ierr = MatSetOption (this->petsc_mat, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);

  this->matrix->close();
  this->matrix->zero();
  
  // set the flag to true
  this->finalized_sparsity = true;
}



template<typename T>
void
FESystemNumerics::PetscMultiBlockSparseMatrix<T>::update()
{
  // update the sparsity if it has not been done so far
  if (!this->finalized_sparsity) 
    this->finalizeSparsity();
 
  // now set the values
  int elem_num = 0, starting_row_id = 0, starting_col_id = 0;
  int n_rows =0, n_cols = 0, my_rows=0;
  std::vector<int> col_ids_vec(50); 
  PetscErrorCode ierr = 0;
  
  this->matrix->zero();
  this->matrix->close();
  
  // next, set the values of the matrix entries
  for (unsigned int i=0; i<this->row_dimension; i++)
    {
      // the size is calculated based on the columns of each matrix in a block-row.
      // get the sparsity pattern for each row
      for (unsigned int j=0; j<this->column_dimension; j++)
        {
          elem_num = i * this->column_dimension + j;
          SparseMatrix<T>* mat = this->sub_matrices[elem_num];
          if (mat != NULL)
            {
              mat->close();
              Mat pmat = dynamic_cast<PetscMatrix<T>*> (mat)->mat();
              // get the row ij indices
              ierr = MatGetSize(pmat, &n_rows, &n_cols);
              CHKERRABORT(FESystem::COMM_WORLD,ierr);
              Assert(n_rows == this->rows_per_row_wise_block[i], ExcInternalError());
              Assert(n_cols == this->columns_per_column_wise_block[j], ExcInternalError());

              const PetscInt *col_ids = PETSC_NULL;
              const PetscScalar *vals = PETSC_NULL;
              // iterate over the rows in this matrix
              for (int k=0; k < n_rows; k++)
                {
                  int n_cols_in_row = 0;
                  ierr = MatGetRow(pmat, k, &n_cols_in_row, &col_ids, &vals);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                  
                  // next, copy the values of the array into a vector and shift them by the 
                  // appropriate value
                  // apply the actual shift to the row and column ids
                  my_rows = k + starting_row_id;
                  if (col_ids_vec.size() < n_cols_in_row)
                    col_ids_vec.resize(n_cols_in_row);
                  for (unsigned int l=0; l < n_cols_in_row; l++)
                    col_ids_vec[l] = col_ids[l] + starting_col_id;
                  // now set the values in the matrix
                  ierr = MatSetValues(this->petsc_mat, 1, &my_rows, 
                                      n_cols_in_row, &(col_ids_vec[0]), vals, INSERT_VALUES);
                  
                  // now restore the matrix row
                  ierr = MatRestoreRow(pmat, k, &n_cols_in_row, &col_ids, &vals);
                  CHKERRABORT(FESystem::COMM_WORLD,ierr);
                }
            }
          
          starting_col_id += this->columns_per_column_wise_block[j];
        }
      starting_row_id += this->rows_per_row_wise_block[i];
      starting_col_id = 0;
    }
  
  
  // not finalize the matrix
  this->matrix->close();

}


template<typename T>
void
FESystemNumerics::PetscMultiBlockSparseMatrix<T>::initRightVector(NumericVector<T>& v)
{
  Assert(this->finalized_sparsity, ExcInternalError());
  
  unsigned int total_columns = 0; 
  for (unsigned int i=0; i < this->column_dimension; i++) 
    total_columns += this->columns_per_column_wise_block[i]; 
      
  v.init(total_columns, false);
}



template<typename T>
void
FESystemNumerics::PetscMultiBlockSparseMatrix<T>::extractFromRightMultiblockVector
(unsigned int col_block_num, NumericVector<T>& in_v, NumericVector<T>& out_v)
{
  // it is assumed that the vectors have been correctly initialized with respect
  // to the global and the block matrices
  
  Assert(this->finalized_sparsity, ExcInternalError());
  Assert(col_block_num < this->row_dimension, ExcInternalError());
  
  // total number of rows and columns 
  unsigned int total_columns = 0, starting_id = 0; 
  for (unsigned int i=0; i < this->column_dimension; i++) 
    {
      total_columns += this->columns_per_column_wise_block[i]; 
      if (i < col_block_num) 
        starting_id += this->columns_per_column_wise_block[i]; 
    }
  
  Assert(in_v.size() == total_columns, ExcInternalError());
  Assert(out_v.size() == this->columns_per_column_wise_block[col_block_num],
         ExcInternalError());

  std::vector<unsigned int> ids(out_v.size());
  Utility::iota (ids.begin(), ids.end(), starting_id);
  
  in_v.create_subvector(out_v, ids);
}



template class FESystemNumerics::MultiBlockSparseMatrix<double>;
template class FESystemNumerics::PetscMultiBlockSparseMatrix<double>;

