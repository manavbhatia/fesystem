// $Id: MatrixBase.h,v 1.1.4.5 2008-04-06 04:01:55 manav Exp $
#ifndef __fesystem_matrix_base_h__
#define __fesystem_matrix_base_h__

#include "FESystem/FESystemExceptions.h"

namespace FESystemNumerics
{
 
  // forward declerations
  template<typename T> class VectorBase;
  
  ///
  /// this class provides a base class for different kinds of matrices
  ///
  template <typename T>
  class MatrixBase
  {
public:
    /// constructor
    MatrixBase();

    
    /**
    * Destructor. Empty.
     */     
    virtual ~MatrixBase();
    
    /// @returns the matrix package
    virtual unsigned int getMatrixPackage() const = 0;
    

    /**
      * Set every element in the matrix to 0.  You must redefine
     * what you mean by zeroing the matrix since it depends on
     * how your values are stored.
     */
    virtual void resize(const unsigned int m,
                        const unsigned int n) = 0;

    /**
    * Set every element in the matrix to 0.  You must redefine
     * what you mean by zeroing the matrix since it depends on
     * how your values are stored.
     */
    virtual void zero() = 0;
    
    /**
      * @returns the \p (i,j) element of the matrix.
     * Since internal data representations may differ, you
     * must redefine this function.
     */
    virtual T el(const unsigned int i,
                 const unsigned int j) const = 0;
    

    /// returns the \p i th row of the matrix in the given vector
    virtual void getRow(const unsigned int i, 
                        FESystemNumerics::VectorBase<T>& vec) = 0;


    /// returns the diagonal of the matrix in the given vector
    virtual void getDiagonal(FESystemNumerics::VectorBase<T>& vec) = 0;
    

    /// adds \p val  to all diagonal elements of the matrix
    virtual void shift(const T val) = 0;

    
    /**
      * sets the given element to \pval 
     */
    virtual void set(const unsigned int i,
                     const unsigned int j,
                     const T val) = 0;
    
    /**
      * sets multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void set(const std::vector<unsigned int>& i,
                     const std::vector<unsigned int>& j,
                     const T* val) = 0;

    /// multiplies the matrix with a diagonal matrix on the left hand side. 
    /// The diagonal matrix is represented as a vector.
    virtual void leftDiagonalScale (const FESystemNumerics::VectorBase<T>& vector) = 0;

    
    /// multiplies the matrix with a diagonal matrix on the right hand side. 
    /// The diagonal matrix is represented as a vector.
    virtual void rightDiagonalScale (const FESystemNumerics::VectorBase<T>& vector) = 0;

    
    /**
      * Performs the operation: (*this) <- M2 * (*this) 
     */
    virtual void leftMultiply (const FESystemNumerics::MatrixBase<T>& matrix) = 0;
    
    /**
      * Performs the operation: (*this) <- M2^T * (*this) 
     */
    virtual void leftMultiplyTranspose (const FESystemNumerics::MatrixBase<T>& matrix) = 0;

    
    /**
      * Performs the operation: (*this) <- (*this) * M3
     */
    virtual void rightMultiply (const FESystemNumerics::MatrixBase<T>& matrix) = 0;
    

    /**
      * Performs the operation: (*this) <- M2 * (*this) 
     */
    virtual void leftMultiplyVector (const FESystemNumerics::VectorBase<T>& vector_in,
                                     FESystemNumerics::VectorBase<T>& result) = 0;

    
    /**
      * Performs the operation: (*this) <- M2 * (*this) 
     */
    virtual void rightMultiplyVector (const FESystemNumerics::VectorBase<T>& vector_in,
                                      FESystemNumerics::VectorBase<T>& result) = 0;
    
    
    /**
      * @returns the row-dimension of the matrix.
     */
    virtual unsigned int m() const=0;
    
    /**
      * @returns the column-dimension of the matrix.
     */
    virtual unsigned int n() const =0;
    

    /**
      * Copies \pmatrix to this matrix
     */
    void operator=(const FESystemNumerics::MatrixBase<T>& matrix)
      {Assert(false,ExcInternalError());}

    
    /**
      * Copies \pmatrix to this matrix
     */
    virtual void copy(const FESystemNumerics::MatrixBase<T>& matrix)=0;

    
    /**
      * Adds \pmatrix to this matrix
     */
    virtual void add (const T val,
                      const FESystemNumerics::MatrixBase<T>& matrix,
                      const bool same_nonzero)=0;

    /**
      * Adds \p factor to given element in the matrix.
     */
    virtual void add (const unsigned int i,
                      const unsigned int j, 
                      const T val)=0;
    
    /**
      * adds multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void add(const std::vector<unsigned int>& i,
                     const std::vector<unsigned int>& j,
                     const T* val) = 0;
    
    /**
      * adds a diagonal matrix times a factor to this matrix. The diagonal matrix
     * is given as a vector.
     */
    virtual void diagonalAdd(const T val,
                             const FESystemNumerics::VectorBase<T>& diag) = 0;

    
    /**
      * multiplies \p factor to every element in the matrix.
     */
    virtual void scale (const T factor) =0;
    
    
    /**
     * prints the values of this matrix to the output
     */
    virtual void print (std::ostream& output);
    
  };
  
}



template<typename T>
FESystemNumerics::MatrixBase<T>::MatrixBase()
{}



template<typename T>
FESystemNumerics::MatrixBase<T>::~MatrixBase()
{}



template<typename T>
void
FESystemNumerics::MatrixBase<T>::print(std::ostream& output)
{
  unsigned int m= this->m(), n = this->n();
  for (unsigned int i=0; i< m; i++)
    {
      for (unsigned int j=0; j< n; j++)
        output << this->el(i,j) << "  ";
      output << std::endl;
    }
}


#endif // __fesystem_matrix_base_h__
