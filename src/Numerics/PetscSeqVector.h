// $Id: PetscSeqVector.h,v 1.1.4.2 2007-03-04 02:46:42 manav Exp $
#ifndef __fesystem_petsc_seq_vector_h__
#define __fesystem_petsc_seq_vector_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Numerics/PetscVectorBase.h"


namespace FESystemNumerics
{
  
  ///
  /// this class provides a base class for different kinds of matrices
  ///
  template <typename T>
  class PetscSeqVector: public FESystemNumerics::PetscVectorBase<T>
  {
public:
    /// constructor
    PetscSeqVector(const unsigned int n=0);
    
    // copy constructor
    PetscSeqVector(const FESystemNumerics::PetscSeqVector<T>& mat);
    
    /**
      * Destructor. Empty.
     */     
    virtual ~PetscSeqVector();
    
    
    /**
      * Copies \p vector to this vector
     */
    virtual void copy (const FESystemNumerics::VectorBase<T>& vector);

    
    
    /// resize the matrix to the given dimensions
    virtual void resize(const unsigned int n);
    
protected:
      
  };
  
}



template<typename T>
FESystemNumerics::PetscSeqVector<T>::PetscSeqVector(const unsigned int n):
FESystemNumerics::PetscVectorBase<T>()
{
  PetscErrorCode ierr=0;
  ierr = VecCreateSeq(PETSC_COMM_SELF,n,&(this->petsc_vec));
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  
  this->zero();
}



template<typename T>
FESystemNumerics::PetscSeqVector<T>::PetscSeqVector
(const FESystemNumerics::PetscSeqVector<T>& vec):
FESystemNumerics::PetscVectorBase<T>()
{
  PetscErrorCode ierr=0;
  ierr = VecDuplicate(vec.getVec(), &(this->petsc_vec));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);

  ierr = VecCopy(vec.getVec(), this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}



template<typename T>
FESystemNumerics::PetscSeqVector<T>::~PetscSeqVector()
{
  PetscErrorCode ierr=0;
  ierr = VecDestroy(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}




template<typename T>
void
FESystemNumerics::PetscSeqVector<T>::resize(const unsigned int n)
{
  if (this->size() == n)
    {
    this->zero();
    return;
    }
  
  PetscErrorCode ierr=0;
  ierr = VecDestroy(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  ierr = VecCreateSeq(PETSC_COMM_SELF,n,&(this->petsc_vec));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  this->zero();
}




template<typename T>
inline
void
FESystemNumerics::PetscSeqVector<T>::copy
(const FESystemNumerics::VectorBase<T>& vector)
{
  Assert(vector.getVectorPackage() == FESystem::PETSC_PACKAGE::num(),
         ExcInternalError());
  
  
  const FESystemNumerics::PetscVectorBase<T>& vec = 
    dynamic_cast<const FESystemNumerics::PetscVectorBase<T>&>(vector);
  const_cast<FESystemNumerics::PetscVectorBase<T>&>(vec).assemble();
  
  // create a temporary matrix to store the result in
  PetscErrorCode ierr = 0;
  
  ierr = VecDestroy(this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  ierr = VecDuplicate(vec.getVec(), &(this->petsc_vec));
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  ierr = VecCopy(vec.getVec(), this->petsc_vec);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
}



//
//template <typename T>
//void 
//PetscVector<T>::add_vector (const FESystemNumerics::PetscSeqVector<T>& V,
//                            const std::vector<unsigned int>& dof_indices)
//{
//  assert (V.size() == dof_indices.size());
//  
//  for (unsigned int i=0; i<V.size(); i++)
//    this->add (dof_indices[i], V.el(i));
//}
//
//


#endif // __fesystem_petsc_seq_vector_h__
