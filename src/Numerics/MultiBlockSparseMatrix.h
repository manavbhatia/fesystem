// $Id:$
/*
 *  MultiBlockSparseMatrix.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/28/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */


#ifndef __fesystem_multi_block_sparse_matrix_h__
#define __fesystem_multi_block_sparse_matrix_h__

// C++ includes
#include <vector>
#include <memory>

// FESystem incudes


// libMesh includes

// petsc includes
#include "petscmat.h"



template <typename T> class SparseMatrix;
template <typename T> class NumericVector;

namespace FESystemNumerics{
  

  /// This class allows construction of a sparse matrix using multiple sparse matrix 
  /// blocks. The sparse matrix blocks should be initialized and given to this matrix.
  /// The sparsity pattern of each one of these matrices will be used to create the
  /// sparsity pattern of this matrix  
  
  template <typename T>
  class MultiBlockSparseMatrix {

  public:
    
    /// constructor
    MultiBlockSparseMatrix();
    
    /// destructor
    virtual ~MultiBlockSparseMatrix();
    
    /// this method clears the data structures
    void clear();
    
    /// sets the number of sub-matrices in the row and column directions respectively
    void setMatrixDimension(const unsigned int row_size, const unsigned int column_size);
    
    
    /// sets the matrix for the given row and colum index of the sub-matrix. For sub-matrices
    /// that will be zero, they should be explicitly set to NULL using the method \psetSubMatrixToNull().
    /// The row and column numbers start at zero.
    void setSubMatrix(const unsigned int row_num, const unsigned int column_num,
                      SparseMatrix<T>& matrix);
    
    /// sets the sub matrix to a null matrix. It is important to call this method for all 
    /// sub-matrices which will be left as zero. The row and column numbers start at zero. 
    void setSubMatrixToNull(const unsigned int row_num, const unsigned int column_num);
    
    /// this returns a reference to the sparse-matrix object that contains the 
    /// multi-block-sparse-matrix object
    SparseMatrix<T>& getSparseMatrixObject();
    
    /// tells the object to update all its values using the sub-matrices. It is important 
    /// to have initialized and updated all values in the sub-matrices before this method is
    /// called.
    virtual void update() = 0;
    
    /// initializes the vector to the same parallel distribution as the matrix so that a
    /// matrix vector product could be carried out
    virtual void initRightVector(NumericVector<T>& vec) = 0;
    
    /// extract the vector from the multiblock vector. The right vector is assumed here
    virtual void extractFromRightMultiblockVector
    (unsigned int col_block_num, NumericVector<T>& in_v, NumericVector<T>& out_v) = 0;

  protected:
    
    /// finalize sparsity so that only value changes are allowed after this
    virtual void finalizeSparsity() = 0;

    /// boolean indicating if sparsity pattern has been finalized
    bool finalized_sparsity;
    
    /// number of sub-matrices in the row direction
    unsigned int row_dimension;
    
    /// number of sub-matrices in the column direction
    unsigned int column_dimension;
    
    /// pointer to the sparse-matrix object that stores the actual matrix
    std::auto_ptr<SparseMatrix<T> > matrix;
    
    /// Vector of sparse matrices that form the sub-matrices of this object.
    /// The matrices are stored here in a logical 2-D array in a row domimated format.
    std::vector<SparseMatrix<T>*> sub_matrices;
    
    /// stores the number of rows in each row-wise block
    std::vector<unsigned int> rows_per_row_wise_block;
    
    /// stores the number of columns in each column-wise block
    std::vector<unsigned int> columns_per_column_wise_block;
        
  };

  
  /// this class works with Petsc matrices to create a compound matrix
  template <typename T>
  class PetscMultiBlockSparseMatrix : public MultiBlockSparseMatrix<T> {

  public:
    /// constructor
    PetscMultiBlockSparseMatrix();
    
    /// destructor
    ~PetscMultiBlockSparseMatrix();
    
    /// updates the matrix from its original submatrices
    virtual void update();
    
    /// initializes the vector to the same parallel distribution as the matrix so that a
    /// matrix vector product could be carried out
    virtual void initRightVector(NumericVector<T>& vec);

    /// extract the vector from the multiblock vector. The right vector is assumed here
    virtual void extractFromRightMultiblockVector
    (unsigned int col_block_num, NumericVector<T>& in_v, NumericVector<T>& out_v);

  protected:
    
    /// finalize sparsity so that only value changes are allowed after this
    virtual void finalizeSparsity();

    /// petsc matrix object for this matrix
    Mat petsc_mat;
    
  };
  
}



#endif // __fesystem_multi_block_sparse_matrix_h__

