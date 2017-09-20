// $Id: PetscVectorBase.h,v 1.1.4.2 2007-03-04 02:46:42 manav Exp $
#ifndef __fesystem_petsc_vector_base_h__
#define __fesystem_petsc_vector_base_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Numerics/VectorBase.h"
#include "FESystem/FESystemController.h"

// petsc includes
#include "petscvec.h"


namespace FESystemNumerics
{
  
  ///
  /// this class provides a base class for different kinds of matrices
  ///
  template <typename T>
  class PetscVectorBase: public FESystemNumerics::VectorBase<T>
  {
public:
    /// constructor
    PetscVectorBase();
    
    
    /**
    * Destructor. Empty.
     */     
    virtual ~PetscVectorBase();
    
    
    /// @returns the matrix package
    virtual unsigned int getVectorPackage() const;
    
        
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
    virtual T el(const unsigned int i) const;
    
    /**
      * sets the given element to \pval 
     */
    virtual void set(const unsigned int i,
                     const T val);
    

    /// raises each element of the vector by the given power
    virtual void power(const double val);
    
    
    /// sets each element of the vector to its reciprocal
    virtual void invertValues();
    
    
    /// returns the dot product of this vector with the given vector
    virtual T dot(const FESystemNumerics::VectorBase<T>& vec) const;

    
    /// adds the value \p val to all elements of the vector
    virtual void shift(const T val);
    
    
    /// sets (*this)[i] = (*this)[i] * vec[i]
    virtual void pointwiseMultiply(const FESystemNumerics::VectorBase<T>& vec);

    
    /**
      * sets multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void set(const std::vector<unsigned int>& i,
                     const T* val);
    
    
    /**
      * @returns the row-dimension of the matrix.
     */
    virtual unsigned int size() const;
    
    
    /**
      * Adds \pvector to this matrix
     */
    virtual void add (const T val,
                      const FESystemNumerics::VectorBase<T>& vector);
    
    
    /**
      * Adds \p factor to given element in the matrix.
     */
    virtual void add (const unsigned int i, 
                      const T val);
    
    /**
      * adds multiple elements with the indeces given in i and j. The 
     * val should be a 2-D matrix of dimensions equal to size of vectors 
     * i and j.
     */
    virtual void add(const std::vector<unsigned int>& i,
                     const T* val);
    
    /**
      * multiplies \p factor to every element in the matrix.
     */
    virtual void scale (const T factor);
    
    /// returns the Petsc matrix data structure
    Vec getVec();
    
    /// @returns a constant Petsc matrix data structure
    const Vec getVec() const;
    
    /// method to flush assemble the matrix
    void assemble();
        

protected:
      
    
    /// an enumeration to keep track of action performed
    enum LastAction {ADD, INSERT, ASSEMBLE};
    
    /// last action performed, to decide whether to flush the matrix operations 
    /// or not
    LastAction vector_last_action;
    
    /// Petsc matrix data structure
    Vec petsc_vec;

  };
  
}



template<typename T>
FESystemNumerics::PetscVectorBase<T>::PetscVectorBase():
FESystemNumerics::VectorBase<T>()
{}



template<typename T>
FESystemNumerics::PetscVectorBase<T>::~PetscVectorBase()
{}


template<typename T>
unsigned int
FESystemNumerics::PetscVectorBase<T>::getVectorPackage() const
{
  return FESystem::PETSC_PACKAGE::num();
}




template<typename T>
void 
FESystemNumerics::PetscVectorBase<T>::zero()
{
  PetscErrorCode ierr;
  ierr = VecZeroEntries(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->assemble();
}


template<typename T>
T 
FESystemNumerics::PetscVectorBase<T>::el(const unsigned int i) const
{
  // if the data has not been flushed, do that before accessing the value
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
    case FESystemNumerics::PetscVectorBase<T>::ADD:
      const_cast<FESystemNumerics::PetscVectorBase<T>*>(this)->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  PetscInt i_num = i;
  T val;
  
  ierr = VecGetValues (this->petsc_vec, 1, &i_num, &val);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  return val;  
}



template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::set(const unsigned int i,
                                          const T val)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  ierr = VecSetValue (this->petsc_vec, i,
                      val, INSERT_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->vector_last_action = FESystemNumerics::PetscVectorBase<T>::INSERT;
}



template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::set(const std::vector<unsigned int>& i,
                                          const T* val)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  ierr = VecSetValues (this->petsc_vec, i.size(), (const PetscInt*) &(i[0]), 
                       val, INSERT_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->vector_last_action = FESystemNumerics::PetscVectorBase<T>::INSERT;  
}



template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::power(const double val)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  T *vec_vals;
  PetscErrorCode ierr;
  ierr = VecGetArray (this->petsc_vec, &vec_vals);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
    
  for (unsigned int i=0; i<this->size(); i++)
    vec_vals[i] = pow(vec_vals[i], val);
  
  ierr = VecGetArray (this->petsc_vec, &vec_vals);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}






template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::invertValues()
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }

  PetscErrorCode ierr;
  ierr = VecReciprocal(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}





template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::shift(const T val)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  ierr = VecShift(this->petsc_vec, val);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}




template<typename T>
T
FESystemNumerics::PetscVectorBase<T>::dot(const FESystemNumerics::VectorBase<T>& vec) const
{
  Assert(vec.getVectorPackage() == FESystem::PETSC_PACKAGE::num(), 
         ExcInternalError());
  
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      const_cast<FESystemNumerics::PetscVectorBase<T>*>(this)->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscVectorBase<T>& pvec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vec);
  
  T val;
  PetscErrorCode ierr;
  ierr = VecDot(this->petsc_vec, pvec.getVec(), &val);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  return val;
}




template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::pointwiseMultiply
(const FESystemNumerics::VectorBase<T>& vec)
{
  Assert(vec.getVectorPackage() == FESystem::PETSC_PACKAGE::num(), 
         ExcInternalError());
  
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscVectorBase<T>& pvec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vec);

  PetscErrorCode ierr;
  ierr = VecPointwiseMult(this->petsc_vec, this->petsc_vec, pvec.getVec());
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
}




template<typename T>
unsigned int
FESystemNumerics::PetscVectorBase<T>::size() const
{
  PetscErrorCode ierr;
  PetscInt num1=0;
  ierr = VecGetSize (this->petsc_vec, &num1);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  return num1;  
}




template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::add (const T val,
                                           const FESystemNumerics::VectorBase<T>& vector)
{
  Assert(vector.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());

  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
    case FESystemNumerics::PetscVectorBase<T>::ADD:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  const FESystemNumerics::PetscVectorBase<T>& vec = 
  dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector);
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).assemble();

  PetscErrorCode ierr;
  ierr = VecAXPY (this->petsc_vec, val, vec.getVec());
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}





template<typename T>
void 
FESystemNumerics::PetscVectorBase<T>::add (const unsigned int i, 
                                           const T val)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  ierr = VecSetValue (this->petsc_vec, i,
                      val, ADD_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->vector_last_action = FESystemNumerics::PetscVectorBase<T>::ADD;    
}





template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::add(const std::vector<unsigned int>& i,
                                          const T* val)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ADD:
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  PetscErrorCode ierr;
  ierr = VecSetValues (this->petsc_vec, i.size(), (const PetscInt*) &(i[0]), 
                      val, ADD_VALUES);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->vector_last_action = FESystemNumerics::PetscVectorBase<T>::ADD;      
}






template<typename T>
void
FESystemNumerics::PetscVectorBase<T>::scale (const T factor)
{
  switch (this->vector_last_action)
    {
    case FESystemNumerics::PetscVectorBase<T>::INSERT:
    case FESystemNumerics::PetscVectorBase<T>::ADD:
      this->assemble();
      break;
      
    case FESystemNumerics::PetscVectorBase<T>::ASSEMBLE:
      break;
    }
  
  
  PetscErrorCode ierr;
  ierr = VecScale (this->petsc_vec, factor);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);    
}



template<typename T>
Vec
FESystemNumerics::PetscVectorBase<T>::getVec()
{
  return this->petsc_vec;
}



template<typename T>
const Vec
FESystemNumerics::PetscVectorBase<T>::getVec() const
{
  return this->petsc_vec;
}



template<typename T>
inline
void
FESystemNumerics::PetscVectorBase<T>::assemble ()
{
  PetscErrorCode ierr;
  ierr = VecAssemblyBegin(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  ierr = VecAssemblyEnd(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);  
  
  this->vector_last_action = FESystemNumerics::PetscVectorBase<T>::ASSEMBLE;
}



#endif // __fesystem_petsc_vector_base_h__
