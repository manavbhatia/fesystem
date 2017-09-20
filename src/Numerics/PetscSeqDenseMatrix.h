// $Id: PetscSeqDenseMatrix.h,v 1.1.4.2 2007-03-04 02:46:42 manav Exp $
#ifndef __fesystem_petsc_seq_dense_matrix_h__
#define __fesystem_petsc_seq_dense_matrix_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Numerics/PetscMatrixBase.h"


namespace FESystemNumerics
{
  
  ///
  /// this class provides a base class for different kinds of matrices
  ///
  template <typename T>
  class PetscSeqDenseMatrix: public FESystemNumerics::PetscMatrixBase<T>
  {
public:
    /// constructor
    PetscSeqDenseMatrix(const unsigned int m=0, 
                        const unsigned int n=0);

    // copy constructor
    PetscSeqDenseMatrix(const FESystemNumerics::PetscSeqDenseMatrix<T>& mat);
    
    /**
    * Destructor. Empty.
     */     
    virtual ~PetscSeqDenseMatrix();
    
    /**
      * Copies \pmatrix to this matrix
     */
    virtual void copy(const FESystemNumerics::MatrixBase<T>& matrix);
    
    
    /**
      * Performs the operation: (*this) <- M2 * (*this). This is reimplemented
     * for the dense matrix class since it needs a different handling of one 
     * parameter
     */
    virtual void leftMultiply (const FESystemNumerics::MatrixBase<T>& matrix);
    
    /**
      * Performs the operation: (*this) <- M2^T * (*this). This is reimplemented
     * for the dense matrix class since it needs a different handling of one 
     * parameter 
     */
    virtual void leftMultiplyTranspose (const FESystemNumerics::MatrixBase<T>& matrix);
    
    
    /**
      * Performs the operation: (*this) <- (*this) * M3. This is reimplemented
     * for the dense matrix class since it needs a different handling of one 
     * parameter 
     */
    virtual void rightMultiply (const FESystemNumerics::MatrixBase<T>& matrix);

    
    /// resize the matrix to the given dimensions
    virtual void resize(const unsigned int m,
                        const unsigned int n);

protected:
          
  };
  
}



template<typename T>
FESystemNumerics::PetscSeqDenseMatrix<T>::PetscSeqDenseMatrix(const unsigned int m,
                                                              const unsigned int n):
FESystemNumerics::PetscMatrixBase<T>()
{
  PetscErrorCode ierr=0;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,m,n,PETSC_NULL,&(this->petsc_mat));
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  
  this->zero();
}



template<typename T>
FESystemNumerics::PetscSeqDenseMatrix<T>::PetscSeqDenseMatrix
(const FESystemNumerics::PetscSeqDenseMatrix<T>& mat):
FESystemNumerics::PetscMatrixBase<T>()
{
  PetscErrorCode ierr=0;
  ierr = MatDuplicate(mat.getMat(),MAT_COPY_VALUES, &(this->petsc_mat));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}



template<typename T>
FESystemNumerics::PetscSeqDenseMatrix<T>::~PetscSeqDenseMatrix()
{
  PetscErrorCode ierr=0;
  ierr = MatDestroy(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




template<typename T>
void
FESystemNumerics::PetscSeqDenseMatrix<T>::resize(const unsigned int m, 
                                                 const unsigned int n)
{
  if (this->m() == m && this->n() == n)
    {
    this->zero();
    return;
    }
  
  PetscErrorCode ierr=0;
  ierr = MatDestroy(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  ierr = MatCreateSeqDense(PETSC_COMM_SELF,m,n,PETSC_NULL,&(this->petsc_mat));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->zero();
}





template<typename T>
inline
void
FESystemNumerics::PetscSeqDenseMatrix<T>::copy
(const FESystemNumerics::MatrixBase<T>& matrix)
{
  Assert(matrix.getMatrixPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  
  
  const FESystemNumerics::PetscMatrixBase<T>& mat = 
    dynamic_cast<const FESystemNumerics::PetscMatrixBase<T>&>(matrix);
  const_cast<FESystemNumerics::PetscMatrixBase<T>&>(mat).finalAssemble();
  
  // create a temporary matrix to store the result in
  PetscErrorCode ierr = 0;

  ierr = MatDestroy(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  ierr = MatDuplicate(mat.getMat(),MAT_COPY_VALUES, &(this->petsc_mat));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}





template<typename T>
inline
void
FESystemNumerics::PetscSeqDenseMatrix<T>::leftMultiply
(const FESystemNumerics::MatrixBase<T>& matrix)
{
  Assert(matrix.getMatrixPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  
  
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      this->finalAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscMatrixBase<T>& mat = 
    dynamic_cast<const FESystemNumerics::PetscMatrixBase<T>&>(matrix);
  const_cast<FESystemNumerics::PetscMatrixBase<T>&>(mat).finalAssemble();
  
  // create a temporary matrix to store the result in
  Mat result_mat;
  PetscErrorCode ierr = 0;
  ierr = MatMatMult (mat.getMat(), this->petsc_mat, MAT_INITIAL_MATRIX, 1.0, 
                     &result_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // since Petsc will itself create the result matrix, we can destroy the matrix 
  // of this object, and store the pointer of result_mat
  ierr = MatDestroy(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->petsc_mat = result_mat;
}




template<typename T>
inline
void
FESystemNumerics::PetscSeqDenseMatrix<T>::leftMultiplyTranspose
(const FESystemNumerics::MatrixBase<T>& matrix)
{
  Assert(matrix.getMatrixPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  
  
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      this->finalAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscMatrixBase<T>& mat = 
    dynamic_cast<const FESystemNumerics::PetscMatrixBase<T>&>(matrix);
  const_cast<FESystemNumerics::PetscMatrixBase<T>&>(mat).finalAssemble();
  
  // create a temporary matrix to store the result in
  Mat result_mat;
  PetscErrorCode ierr = 0;
  ierr = MatMatMultTranspose(mat.getMat(), this->petsc_mat, 
                             MAT_INITIAL_MATRIX, 1.0, &result_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // since Petsc will itself create the result matrix, we can destroy the matrix 
  // of this object, and store the pointer of result_mat
  ierr = MatDestroy(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->petsc_mat = result_mat;
}





template<typename T>
inline
void
FESystemNumerics::PetscSeqDenseMatrix<T>::rightMultiply
(const FESystemNumerics::MatrixBase<T>& matrix)
{
  
  Assert(matrix.getMatrixPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      this->finalAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscMatrixBase<T>& mat = 
    dynamic_cast<const FESystemNumerics::PetscMatrixBase<T>&>(matrix);
  const_cast<FESystemNumerics::PetscMatrixBase<T>&>(mat).finalAssemble();
  
  // create a temporary matrix to store the result in
  Mat result_mat;
  PetscErrorCode ierr = 0;
  ierr = MatMatMult (this->petsc_mat, mat.getMat(), MAT_INITIAL_MATRIX, 1.0, &result_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // since Petsc will itself create the result matrix, we can destroy the matrix 
  // of this object, and store the pointer of result_mat
  ierr = MatDestroy(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->petsc_mat = result_mat;
}


//template <typename T>
//void PetscMatrix<T>::add_matrix(const FESystemNumerics::PetscSeqDenseMatrix<T>& dm,
//                                const std::vector<unsigned int>& dof_indices)
//{
//  this->add_matrix (dm, dof_indices, dof_indices);
//}




//
//template <typename T>
//void 
//PetscMatrix<T>::add_matrix(const FESystemNumerics::PetscSeqDenseMatrix<T>& dm,
//                           const std::vector<unsigned int>& rows,
//                           const std::vector<unsigned int>& cols)
//{
//  assert (this->initialized());
//  
//  const unsigned int m = dm.m();
//  const unsigned int n = dm.n();
//  
//  assert (rows.size() == m);
//  assert (cols.size() == n);
//  
//  int ierr=0;
//  
//  // get the array of the matrix
//  T *array;
//  ierr = MatGetArray(dm.getMat(), &array);
//  CHKERRABORT(FESystem::COMM_WORLD,ierr);
//  
//  
//  // These casts are required for PETSc <= 2.1.5
//  ierr = MatSetValues(_mat,
//                      m, (int*) &rows[0],
//                      n, (int*) &cols[0],
//                      (PetscScalar*) &dm.get_values()[0],
//                      ADD_VALUES);
//  CHKERRABORT(libMesh::COMM_WORLD,ierr);
//  
//  // restore the array
//  ierr = MatRestoreArray(dm.getMat(), &array);
//  CHKERRABORT(FESystem::COMM_WORLD,ierr);
//}
//

#endif // __fesystem_petsc_seq_dense_matrix_h__
