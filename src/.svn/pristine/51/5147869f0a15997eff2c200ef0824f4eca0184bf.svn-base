// $Id: TensorBase.h,v 1.2 2006-09-05 20:41:45 manav Exp $

#ifndef __fesystem__tensor_base_h__
#define __fesystem__tensor_base_h__

// C / C++ includes
#include <vector>
#include <cmath>

// FESystem includes
#include "FESystem/FESystemExceptions.h"


// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


/**
*  This is a base class for tensors, and provides some pure abstract methods.
 * All tensors will inherit from this base classs.
 * Through the definition of this class, all tensors can be referred to by a 
 * pointer or reference to this base class.
 */
class TensorBase
{
public:

  /// destrutor
  virtual ~TensorBase() = 0;
  
  /// @returns the dimension of this tensor
  virtual unsigned int getDim() const = 0;
  
  /// @returns the rank of this tensor
  virtual unsigned int getRank() const = 0;
  
  /**
    * Return the Frobenius-norm of a
   * tensor, i.e. the square root
   * of the sum of squares of all
   * entries. For the present case
   * of rank-1 tensors, this equals
   * the usual
   * <tt>l<sub>2</sub></tt> norm of
   * the vector.
   */
  virtual double norm () const =0;
  
  /**
    * Return the square of the
   * Frobenius-norm of a tensor,
   * i.e. the square root of the
   * sum of squares of all entries.
   *
   * This function mainly exists
   * because it makes computing the
   * norm simpler recursively, but
   * may also be useful in other
   * contexts.
   */
  virtual double norm_square () const=0;
  
  /**
    * Reset all values to zero.
   *
   * Note that this is partly inconsistent
   * with the semantics of the @p clear()
   * member functions of the STL and of
   * several other classes within deal.II
   * which not only reset the values of
   * stored elements to zero, but release
   * all memory and return the object into
   * a virginial state. However, since the
   * size of objects of the present type is
   * determined by its template parameters,
   * resizing is not an option, and indeed
   * the state where all elements have a
   * zero value is the state right after
   * construction of such an object.
   */
  virtual void clear () =0;
  
};



#endif // __fesystem_tensor_base_h__


