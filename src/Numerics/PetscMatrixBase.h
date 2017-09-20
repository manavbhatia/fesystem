// $Id: PetscMatrixBase.h,v 1.1.4.2 2007-03-04 02:46:42 manav Exp $
#ifndef __fesystem_petsc_matrix_base_h__
#define __fesystem_petsc_matrix_base_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Numerics/MatrixBase.h"
#include "FESystem/FESystemController.h"
#include "Numerics/PetscVectorBase.h"

// Petsc includes
#include "petscmat.h"


namespace FESystemNumerics
{
  
  ///
  /// this class provides a base class for different kinds of matrices
  ///
  template <typename T>
  class PetscMatrixBase: public MatrixBase<T>
  {
public:
    /// constructor
    PetscMatrixBase();
    
    
    /**
    * Destructor. Empty.
     */     
    virtual ~PetscMatrixBase();
    

    /// @returns the matrix package
    virtual unsigned int getMatrixPackage() const;
    
    /**
      * resize the matrix in the given dimensions
     */
    virtual void resize(const unsigned int m,
                        const unsigned int n) = 0;
    
    /**
      * Set every element in the matrix to 0.  You must redefine
     * what you mean by zeroing the matrix since it depends on
     * how your values are stored.
     */
    virtual void zero();
    
    /**
      * @returns the \p (i,j) element of the matrix.
     * Since internal data representations may differ, you
     * must redefine this function.
     */
    virtual T el(const unsigned int i,
                 const unsigned int j) const;
    

    /// returns the \p i th row of the matrix in the given vector
    virtual void getRow(const unsigned int i, 
                        FESystemNumerics::VectorBase<T>& vec);
    

    /// returns the diagonal of the matrix in the given vector
    virtual void getDiagonal(FESystemNumerics::VectorBase<T>& vec);

    
    // adds \p val  to all diagonal elements of the matrix
    virtual void shift(const T val);

    
    /**
      * sets the given element to \pval 
     */
    virtual void set(const unsigned int i,
                     const unsigned int j,
                     const T val);
    
    /**
      * sets multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void set(const std::vector<unsigned int>& i,
                     const std::vector<unsigned int>& j,
                     const T* val);

    /**
      * Adds \pmatrix to this matrix
     */
    virtual void setRow(const unsigned row_num,
                        const FESystemNumerics::PetscVectorBase<T>& vector);
    
    /**
      * Adds \pmatrix to this matrix
     */
    virtual void setRow(const unsigned row_num,
                        const std::vector<unsigned int>& col_ids,
                        const FESystemNumerics::PetscVectorBase<T>& vector);

    
    /// multiplies the matrix with a diagonal matrix on the left hand side. 
    /// The diagonal matrix is represented as a vector.
    virtual void leftDiagonalScale (const FESystemNumerics::VectorBase<T>& vector);
    
    
    /// multiplies the matrix with a diagonal matrix on the right hand side. 
    /// The diagonal matrix is represented as a vector.
    virtual void rightDiagonalScale (const FESystemNumerics::VectorBase<T>& vector);

    
    /**
      * Performs the operation: (*this) <- M2 * (*this) 
     */
    virtual void leftMultiply (const FESystemNumerics::MatrixBase<T>& matrix);
    
    /**
      * Performs the operation: (*this) <- M2^T * (*this) 
     */
    virtual void leftMultiplyTranspose (const FESystemNumerics::MatrixBase<T>& matrix);
    
    
    /**
      * Performs the operation: (*this) <- (*this) * M3
     */
    virtual void rightMultiply (const FESystemNumerics::MatrixBase<T>& matrix);
    
    
    /**
      * Performs the operation: (*this) <- M2 * (*this) 
     */
    virtual void leftMultiplyVector (const FESystemNumerics::VectorBase<T>& vector_in,
                                     FESystemNumerics::VectorBase<T>& result);
    
    
    /**
      * Performs the operation: (*this) <- M2 * (*this) 
     */
    virtual void rightMultiplyVector (const FESystemNumerics::VectorBase<T>& vector_in,
                                      FESystemNumerics::VectorBase<T>& result);
    
    
    /**
      * @returns the row-dimension of the matrix.
     */
    virtual unsigned int m() const;
    
    /**
      * @returns the column-dimension of the matrix.
     */
    virtual unsigned int n() const;
    

    /**
      * Adds \pmatrix to this matrix
     */
    virtual void add (const T val, 
                      const FESystemNumerics::MatrixBase<T>& matrix,
                      const bool same_nonzer);

    
    /**
      * Adds \p factor to every element in the matrix.
     */
    virtual void add (const unsigned int i, 
                      const unsigned int j,
                      const T factor);
    
    /**
      * adds multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void add(const std::vector<unsigned int>& i,
                     const std::vector<unsigned int>& j,
                     const T* val);
    
    /**
      * adds a diagonal matrix times a factor to this matrix. The diagonal matrix
     * is given as a vector.
     */
    virtual void diagonalAdd(const T val,
                             const FESystemNumerics::VectorBase<T>& diag);

    
    //    /**
//      * Adds \pmatrix to this matrix
//     */
//    virtual void addRowVector (const unsigned row_num,
//                               const T val,
//                               const FESystemNumerics::PetscVectorBase<T>& vector);
//    
//    /**
//      * Adds \pmatrix to this matrix
//     */
//    virtual void addColumnVector (const unsigned col_num,
//                                  const T val,
//                                  const FESystemNumerics::PetscVectorBase<T>& vector);
    

    
    /**
      * multiplies \p factor to every element in the matrix.
     */
    virtual void scale (const T factor);

    /// returns the Petsc matrix data structure
    Mat getMat();
    
    /// @returns a constant Petsc matrix data structure
    const Mat getMat() const;

    /// method to flush assemble the matrix
    void flushAssemble();
    
    /// method to perform a final assembly of the matrix
    void finalAssemble();
    

protected:
      
    
    /// an enumeration to keep track of action performed
    enum LastAction {ADD, INSERT, FLUSH_ASSEMBLE, FINAL_ASSEMBLE};
    
    /// last action performed, to decide whether to flush the matrix operations 
    /// or not
    LastAction matrix_last_action;
      
    /// Petsc matrix data structure
    Mat petsc_mat;
    
  };
  
}



template<typename T>
FESystemNumerics::PetscMatrixBase<T>::PetscMatrixBase():
FESystemNumerics::MatrixBase<T>()
{}



template<typename T>
FESystemNumerics::PetscMatrixBase<T>::~PetscMatrixBase()
{}



template<typename T>
unsigned int
FESystemNumerics::PetscMatrixBase<T>::getMatrixPackage() const
{
  return FESystem::PETSC_PACKAGE::num();
}



template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::zero()
{
  PetscErrorCode ierr;
  ierr = MatZeroEntries(this->petsc_mat);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->finalAssemble();
}



template<typename T>
inline
T
FESystemNumerics::PetscMatrixBase<T>::el(const unsigned int i,
                                         const unsigned int j) const
{
  // if the data has not been flushed, do that before accessing the value
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      const_cast<FESystemNumerics::PetscMatrixBase<T>*>(this)->finalAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  PetscInt i_num = i;
  PetscInt j_num = j;
  T val;

  ierr = MatGetValues (this->petsc_mat, 1, &i_num, 1, &j_num, (PetscScalar*) &val);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  return val;
}



template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::getRow(const unsigned int i, 
                                             FESystemNumerics::VectorBase<T>& vec)
{
  Assert(vec.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());

  Assert (vec.size() == this->n(), ExcInternalError());
    
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

  PetscErrorCode ierr;
  unsigned int n_cols, *cols;
  T *row_vals, *vec_vals;
  
  // get the row
  ierr = MatGetRow (this->petsc_mat, i, (PetscInt*) &(n_cols), (const PetscInt**) &(cols), 
                    (const PetscScalar**) &(row_vals));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // copy values in the vector
  vec.zero();
  Vec pvec = dynamic_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).getVec();
  ierr = VecGetArray (pvec, &vec_vals);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  for (unsigned int j=0; j < n_cols; j++)
    vec_vals[cols[j]] = row_vals[j];

  ierr = VecRestoreArray (pvec, &vec_vals);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  ierr = MatRestoreRow (this->petsc_mat, i, (PetscInt*) &(n_cols), (const PetscInt**) &(cols), 
                        (const PetscScalar**) &(row_vals));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::getDiagonal(FESystemNumerics::VectorBase<T>& vec)
{
  Assert(vec.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  
  Assert (vec.size() == this->n(), ExcInternalError());
  
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
  
  PetscErrorCode ierr;

  Vec pvec = dynamic_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).getVec();

  // get the row
  ierr = MatGetDiagonal(this->petsc_mat, pvec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::shift(const T val)
{
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

  PetscErrorCode ierr;
  ierr = MatShift (this->petsc_mat, val);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}





template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::set(const unsigned int i,
                                       const unsigned int j,
                                       const T val)
{
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      break;
    }

  PetscErrorCode ierr;
  ierr = MatSetValue (this->petsc_mat, i, j,
                       val, INSERT_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::INSERT;
}




template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::set(const std::vector<unsigned int>& i,
                                       const std::vector<unsigned int>& j,
                                       const T* val)
{
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      break;
    }

  PetscErrorCode ierr;
  ierr = MatSetValues (this->petsc_mat, i.size(), (const PetscInt*) &(i[0]), 
                       j.size(),  (const PetscInt*) &(j[0]),
                      val, INSERT_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::INSERT;
}




template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::setRow
(const unsigned int row_num,
 const FESystemNumerics::PetscVectorBase<T>& vector)
{
  // make sure that the dimensions agree
  Assert (row_num < this->m(), ExcInternalError());
  Assert (vector.size() == this->n(), ExcInternalError());  
  
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscVectorBase<T>& vec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector);
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).assemble();
  Vec pvec = const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).getVec();
  
  // get the vector array
  PetscErrorCode ierr;
  T *array = NULL;
  
  ierr = VecGetArray(pvec, &array);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  // now set the rows in the matrix
  ierr = MatSetValuesRow(this->petsc_mat, row_num, array);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  // now replace the vector array
  ierr = VecRestoreArray(pvec, &array);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::INSERT;
}




template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::setRow
(const unsigned int row_num,
 const std::vector<unsigned int>& col_ids,
 const FESystemNumerics::PetscVectorBase<T>& vector)
{
  // make sure that the dimensions agree
  Assert (row_num < this->m(), ExcInternalError());
  Assert (vector.size() == col_ids.size(), ExcInternalError());  
  
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscVectorBase<T>& vec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector);
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).assemble();
  Vec pvec = const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).getVec();
  
  // get the vector array
  PetscErrorCode ierr;
  T *array = NULL;
  
  ierr = VecGetArray(pvec, &array);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  std::vector<unsigned int> row_id_vec;
  row_id_vec.push_back(row_num);
  this->set(row_id_vec, col_ids, array);
  
  // now replace the vector array
  ierr = VecRestoreArray(pvec, &array);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::INSERT;
}



template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::leftDiagonalScale
(const FESystemNumerics::VectorBase<T>& vector)
{
  Assert(vector.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
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
  
  const FESystemNumerics::PetscVectorBase<T>& vec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector);
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).assemble();
  
  // create a temporary matrix to store the result in
  PetscErrorCode ierr = 0;
  ierr = MatDiagonalScale (this->petsc_mat, vec.getVec(), PETSC_NULL);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}



template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::rightDiagonalScale
(const FESystemNumerics::VectorBase<T>& vector)
{
  Assert(vector.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
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
  
  const FESystemNumerics::PetscVectorBase<T>& vec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector);
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).assemble();
  
  // create a temporary matrix to store the result in
  PetscErrorCode ierr = 0;
  ierr = MatDiagonalScale (this->petsc_mat, PETSC_NULL, vec.getVec());
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}





template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::leftMultiply(const FESystemNumerics::MatrixBase<T>& matrix)
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
  ierr = MatMatMult (mat.getMat(), this->petsc_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, 
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
FESystemNumerics::PetscMatrixBase<T>::leftMultiplyTranspose
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
                             MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result_mat);
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
FESystemNumerics::PetscMatrixBase<T>::rightMultiply(const FESystemNumerics::MatrixBase<T>& matrix)
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
  ierr = MatMatMult (this->petsc_mat, mat.getMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result_mat);
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
FESystemNumerics::PetscMatrixBase<T>::rightMultiplyVector
(const FESystemNumerics::VectorBase<T>& vector_in,
 FESystemNumerics::VectorBase<T>& result)
{
  
  Assert(vector_in.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  Assert(result.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
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
  
  const FESystemNumerics::PetscVectorBase<T>& vec_in = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector_in);
  FESystemNumerics::PetscVectorBase<T>& vec_res = 
    dynamic_cast<FESystemNumerics::PetscVectorBase<T>&>(result);

  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec_in).assemble();
  vec_res.zero();
  
  // create a temporary matrix to store the result in
  PetscErrorCode ierr = 0;
  ierr = MatMult(this->petsc_mat, vec_in.getVec(), vec_res.getVec());
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::leftMultiplyVector
(const FESystemNumerics::VectorBase<T>& vector_in,
 FESystemNumerics::VectorBase<T>& result)
{
  
  Assert(vector_in.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  Assert(result.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
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
  
  const FESystemNumerics::PetscVectorBase<T>& vec_in = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector_in);
  FESystemNumerics::PetscVectorBase<T>& vec_res = 
    dynamic_cast<FESystemNumerics::PetscVectorBase<T>&>(result);
  
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec_in).assemble();
  vec_res.zero();

  // create a temporary matrix to store the result in
  PetscErrorCode ierr = 0;
  ierr = MatMultTranspose(this->petsc_mat, vec_in.getVec(), vec_res.getVec());
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




template<typename T>
inline
unsigned int
FESystemNumerics::PetscMatrixBase<T>::m() const
{
  PetscErrorCode ierr;
  PetscInt num1=0, num2=0;
  ierr = MatGetSize (this->petsc_mat, &num1, &num2);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  return num1;
}



template<typename T>
inline
unsigned int
FESystemNumerics::PetscMatrixBase<T>::n() const
{
  PetscErrorCode ierr;
  PetscInt num1=0, num2=0;
  ierr = MatGetSize (this->petsc_mat, &num1, &num2);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  return num2;
}



template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::add (const T val, 
                                           const FESystemNumerics::MatrixBase<T>& matrix,
                                           const bool same_nonzero)
{
  Assert(matrix.getMatrixPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());

  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      this->finalAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscMatrixBase<T>& mat = 
    dynamic_cast<const FESystemNumerics::PetscMatrixBase<T>&>(matrix);
  const_cast<FESystemNumerics::PetscMatrixBase<T>&>(mat).finalAssemble();

  PetscErrorCode ierr;
  if (same_nonzero)
    ierr = MatAXPY(this->petsc_mat, val, mat.getMat(), SAME_NONZERO_PATTERN); 
  else
    ierr = MatAXPY(this->petsc_mat, val, mat.getMat(), DIFFERENT_NONZERO_PATTERN); 
  
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}





template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::add (const unsigned int i, 
                                           const unsigned int j,
                                           const T val)
{
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }

  PetscErrorCode ierr;
  ierr = MatSetValue (this->petsc_mat, i, j,
                      val, ADD_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::ADD;
}





template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::add(const std::vector<unsigned int>& i,
                                       const std::vector<unsigned int>& j,
                                       const T* val)
{
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }

  PetscErrorCode ierr;
  ierr = MatSetValues (this->petsc_mat, i.size(), (const PetscInt*) &(i[0]), 
                       j.size(), (const PetscInt*) &(j[0]),
                       val, ADD_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::ADD;
}




template<typename T>
inline
void 
FESystemNumerics::PetscMatrixBase<T>::diagonalAdd
(const T val,
 const FESystemNumerics::VectorBase<T> & vec)
{
  Assert (vec.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  Assert (this->m() == this->n(), ExcInternalError());
  Assert (vec.size() == this->m(), ExcInternalError());
  
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
      this->flushAssemble();
      break;
      
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscVectorBase<T> & pvec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T> &>(vec);
  const_cast<FESystemNumerics::PetscVectorBase<T> &>(pvec).assemble();
  Vec petsc_vector =   const_cast<FESystemNumerics::PetscVectorBase<T> &>(pvec).getVec();
  
  PetscErrorCode ierr;
  T *vec_vals;

  ierr = VecGetArray (petsc_vector, &vec_vals);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  for (unsigned int i=0; i<this->m(); i++)
    this->add(i,i, val * vec_vals[i]);
  
  ierr = VecRestoreArray (petsc_vector, &vec_vals);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}



template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::scale (const T factor)
{
  switch (this->matrix_last_action)
    {
    case FESystemNumerics::PetscMatrixBase<T>::INSERT:
    case FESystemNumerics::PetscMatrixBase<T>::ADD:
    case FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE:
      this->finalAssemble();
      break;

    case FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE:
      break;
    }
  

  PetscErrorCode ierr;
  ierr = MatScale (this->petsc_mat, factor);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}




template<typename T>
inline
Mat
FESystemNumerics::PetscMatrixBase<T>::getMat()
{
  return this->petsc_mat;
}



template<typename T>
inline
const Mat
FESystemNumerics::PetscMatrixBase<T>::getMat() const
{
  return this->petsc_mat;
}


template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::flushAssemble ()
{
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(this->petsc_mat,MAT_FLUSH_ASSEMBLY);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  ierr = MatAssemblyEnd(this->petsc_mat,MAT_FLUSH_ASSEMBLY);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::FLUSH_ASSEMBLE;
}


template<typename T>
inline
void
FESystemNumerics::PetscMatrixBase<T>::finalAssemble ()
{
  PetscErrorCode ierr;
  ierr = MatAssemblyBegin(this->petsc_mat,MAT_FINAL_ASSEMBLY);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  ierr = MatAssemblyEnd(this->petsc_mat,MAT_FINAL_ASSEMBLY);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  this->matrix_last_action = FESystemNumerics::PetscMatrixBase<T>::FINAL_ASSEMBLE;
}

#endif // __fesystem_petsc_matrix_base_h__
