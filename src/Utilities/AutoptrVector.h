// $Id: AutoptrVector.h,v 1.1.2.2 2007-05-08 05:19:13 manav Exp $

#ifndef __fesystem_auto_ptr_vector_h__
#define __fesystem_auto_ptr_vector_h__

// C++ includes
#include <vector>

// FESystem includes
#include "FESystem/FESystemExceptions.h"

namespace FESystemUtility
{
  template <typename T> class AutoPtrVector
    {
    public:
      AutoPtrVector():
	vec(NULL)
	{
	  vec = new std::vector<T*>(); 
	}


      AutoPtrVector(const unsigned int n):
	vec(NULL)
	{
	  vec = new std::vector<T*>();
	  vec->resize(n);
	}
      
      
      
      
      AutoPtrVector(FESystemUtility::AutoPtrVector<T>& auto_vec):
	vec(NULL)
	{
	  Assert(auto_vec.get() != NULL, ExcInternalError());
	  this->vec = auto_vec.release();
	}
      


      AutoPtrVector(std::vector<T*>* vec_ptr):
	vec(NULL)
	{
	  Assert(vec_ptr != NULL, ExcInternalError());
	  this->vec = vec_ptr;
	}



      ~AutoPtrVector()
	{
	  if (this->vec != NULL)
	    {
	      // iterate over all elements and delete them
	      for (unsigned int i=0; i < this->vec->size(); i++)
		delete (*this->vec)[i];
	    
	      // now delete the vector
	      delete this->vec;
	      this->vec = NULL;
	    }
	}


      /// returns the number of elements in the vector. This cannot be called if the
      /// vector has been released from this object.
      unsigned int size()
	{
	  Assert(this->vec != NULL, ExcInternalError());
	  return this->vec->size();
	}

      
      /// clears all the elements in the vector, and reduces its size to zero
      void clear()
	{
	  // this cannot be called if the vector has been released
	  Assert (this->vec != NULL, ExcInternalError());

	  // iterate over all elements and delete them
	  for (unsigned int i=0; i < this->vec->size(); i++)
	    delete (*this->vec)[i];
	  this->vec->clear();
	}

      
      /// if the vector was released from this object, this method can be used 
      /// to attach another vector to this object
      void reset(std::vector<T*>* new_vec)
	{
	  Assert(new_vec != NULL, ExcInternalError());
	
	  if (this->vec != NULL)
	    {
	      this->clear();
	      delete this->vec;
	    }
	
	  this->vec = new_vec;
	}


      /// this releases the ownership of the vector from this object, and sets the 
      /// pointer of this object to NULL
      std::vector<T*>* release()
	{
	  Assert(this->vec != NULL, ExcInternalError());
	  std::vector<T*>* new_vec = this->vec;
	
	  this->vec = NULL;
	  return new_vec;
	}


      /// this clears the contents of the vector and resizes it to the given dimensions. This cannot
      /// be called if the vector has been released
      void resize(const unsigned int n)
	{
	  Assert(this->vec != NULL, ExcInternalError());
	  this->clear();

	  this->vec->resize(n);
	}

      /// this attaches another pointer to this end of the vector
      void push_back(T* elem)
	{
	  Assert(this->vec != NULL, ExcInternalError());
	  Assert(elem != NULL, ExcInternalError());

	  this->vec->push_back(elem);
	}


      /// this returns the nth element from the vector. 
      T* operator[](const unsigned int n)
	{
	  Assert(this->vec != NULL, ExcInternalError());
	  Assert(n < this->vec->size(), ExcInternalError());

	  return (*this->vec)[n];
	}
      

      /// this resets the data at the ith location of the vector
      void reset(unsigned int i, T* ptr)
	{
	  Assert(i < this->size(), ExcInternalError());
	  Assert(ptr != NULL, ExcInternalError());
	  Assert(this->vec != NULL, ExcInternalError());
	  
	  if ((*this->vec)[i] != NULL)
	    delete (*this->vec)[i];
	  
	  (*this->vec)[i] = ptr;
	}

      
      /// this returns a constant pointer to the vector of this object. 
      const std::vector<T*>* get()
	{
	  return this->vec;
	}

      
      /// this returns a reference to the vector of this object. 
      std::vector<T*>& getReference()
	{
	  return *(this->vec);
	}

    protected:
     
      /// the pointer of vector for this object
      std::vector<T*>* vec;
    };
}



#endif // __fesystem_auto_ptr_vector_h__
