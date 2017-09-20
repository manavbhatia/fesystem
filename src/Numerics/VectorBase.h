// $Id: VectorBase.h,v 1.1.4.4 2008-04-06 04:01:55 manav Exp $
#ifndef __fesystem_vector_base_h__
#define __fesystem_vector_base_h__



namespace FESystemNumerics
{
  
  ///
  /// this class provides a base class for different kinds of matrices
  ///
  template <typename T>
  class VectorBase
  {
public:
    /// constructor
    VectorBase();
    
    
    /**
    * Destructor. Empty.
     */     
    virtual ~VectorBase();
    
    
    /// @returns the matrix package
    virtual unsigned int getVectorPackage() const = 0;
    
    
    /**
      * Copies \p vector to this vector
     */
    void operator= (const FESystemNumerics::VectorBase<T>& vector)
      {Assert(false, ExcInternalError());}

    
    /**
      * Copies \p vector to this vector
     */
    virtual void copy(const FESystemNumerics::VectorBase<T>& vector) = 0;

    
    /**
      * Set every element in the matrix to 0.  You must redefine
     * what you mean by zeroing the matrix since it depends on
     * how your values are stored.
     */
    virtual void resize(const unsigned int n) = 0;
    
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
    virtual T el(const unsigned int i) const = 0;

    
    /// raises each element of the vector by the given power
    virtual void power(const double val) = 0;
    

    /// sets each element of the vector to its reciprocal
    virtual void invertValues() = 0;


    /// returns the dot product of this vector with the given vector
    virtual T dot(const FESystemNumerics::VectorBase<T>& vec) const = 0;
    
    
    /// adds the value \p val to all elements of the vector
    virtual void shift(const T val)=0;
    
    
    /// sets (*this)[i] = (*this)[i] * vec[i]
    virtual void pointwiseMultiply(const FESystemNumerics::VectorBase<T>& vec) = 0;
    
    /**
      * sets the given element to \pval 
     */
    virtual void set(const unsigned int i,
                     const T val) = 0;
    
    /**
      * sets multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void set(const std::vector<unsigned int>& i,
                     const T* val) = 0;
    
    
    /**
      * @returns the row-dimension of the matrix.
     */
    virtual unsigned int size() const=0;
    
    
    /**
      * Adds \pmatrix to this matrix
     */
    virtual void add (const T val,
                      const FESystemNumerics::VectorBase<T>& vector)=0;
    
    /**
      * Adds \p factor to given element in the matrix.
     */
    virtual void add (const unsigned int i, 
                      const T val)=0;
    
    /**
      * adds multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void add(const std::vector<unsigned int>& i,
                     const T* val) = 0;
    
    /**
      * multiplies \p factor to every element in the matrix.
     */
    virtual void scale (const T factor) =0;
    
    
    /**
     * prints the values of the vector to the output
     */
    virtual void print (std::ostream& output);
    
  };
  
}



template<typename T>
FESystemNumerics::VectorBase<T>::VectorBase()
{}



template<typename T>
FESystemNumerics::VectorBase<T>::~VectorBase()
{}


template<typename T>
void
FESystemNumerics::VectorBase<T>::print(std::ostream& output)
{
  unsigned int m= this->size();
  for (unsigned int i=0; i< m; i++)
    output << this->el(i) << std::endl;
}





#endif // __fesystem_vector_base_h__
