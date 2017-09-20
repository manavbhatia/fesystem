// $Id: TensorRankZero.h,v 1.2 2006-09-05 20:41:45 manav Exp $

#ifndef __fesystem__tensor_rank_zero_h__
#define __fesystem__tensor_rank_zero_h__

// C / C++ includes
#include <vector>
#include <cmath>

// deal.II inlcudes
#include "base/exceptions.h"

// FESystem includes
#include "Numerics/TensorBase.h"

// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


/**
* This class is a specialized version of the <tt>Tensor<rank,dim></tt> class.
 * It handles tensors with no index, i.e. scalars, of fixed dimension and
 * provides the basis so that a uniform interface is available for tensor
 * quantities, even in a situation that it has zero rank. This provides the 
 * same interface as the tensor class of deal.II.
 */
template <>
class Tensor<0,0>: public TensorBase
{
public:
  /**
  * Provide a way to get the
   * dimension of an object without
   * explicit knowledge of it's
   * data type. Implementation is
   * this way instead of providing
   * a function <tt>dimension()</tt>
   * because now it is possible to
   * get the dimension at compile
   * time without the expansion and
   * preevaluation of an inlined
   * function; the compiler may
   * therefore produce more
   * efficient code and you may use
   * this value to declare other
   * data types.
   */
  static const unsigned int dimension = 0;
  
  /**
  * Publish the rank of this tensor to
   * the outside world.
   */
  static const unsigned int rank      = 0;
  
  /**
    * Type of stored objects. This
   * is a double for a rank 1 tensor.
   */
  
  typedef double value_type;
  
  /**
    * Declare an array type which can
   * be used to initialize statically
   * an object of this type.
   *
   * Avoid warning about zero-sized
   * array for <tt>dim==0</tt> by
   * choosing lunatic value that is
   * likely to overflow memory
   * limits.
   */
  typedef double array_type;
  
  /**
    * Constructor. Initialize all entries
   * to zero if <tt>initialize==true</tt>; this
   * is the default behaviour.
   */
  explicit Tensor (const bool initialize = true);
  
  /**
    * Copy constructor, where the
   * data is copied from a C-style
   * array.
   */
  Tensor (const array_type &initializer);
  
  /**
    * Copy constructor.
   */
  Tensor (const Tensor<0,0> &);
  
  /// @returns the dimension of this tensor
  unsigned int getDim() const;
  
  /// @returns the rank of this tensor
  unsigned int getRank() const;
  
  /**
    * Read access to the value of this tensor. Note that the value of 
   * the index should be zero in all cases, since this is a scalar variable
   */
  double   operator [] (const unsigned int index) const;
  
  /**
    * Read and write access to the value of this tensor. Note that the value of
   * this index should be zero in all cases, since this is a scalar variable.
   */
  double & operator [] (const unsigned int index);
  
  /**
    * Assignment operator.
   */
  Tensor<0,0> & operator = (const Tensor<0,0> &);
  
  /**
    * Test for equality of two
   * tensors.
   */
  bool operator == (const Tensor<0,0> &) const;
  
  /**
    * Test for inequality of two
   * tensors.
   */
  bool operator != (const Tensor<0,0> &) const;
  
  /**
    * Add another vector, i.e. move
   * this point by the given
   * offset.
   */
  Tensor<0,0> & operator += (const Tensor<0,0> &);

  /**
    * Add another vector, i.e. move
   * this point by the given
   * offset.
   */
  Tensor<0,0> & operator += (const double &);
  
  
  /**
    * Subtract another vector.
   */
  Tensor<0,0> & operator -= (const Tensor<0,0> &);
  
  
  /**
    * Subtract another vector.
   */
  Tensor<0,0> & operator -= (const double &);
  
  
  /**
    * Scale the vector by
   * <tt>factor</tt>, i.e. multiply all
   * coordinates by <tt>factor</tt>.
   */
  Tensor<0,0> & operator *= (const double factor);
  
  /**
    * Scale the vector by <tt>1/factor</tt>.
   */
  Tensor<0,0> & operator /= (const double factor);
  
  /**
    * Returns the scalar product of
   * two zero rank tensors.
   */
  double          operator * (const Tensor<0,0> &) const;
  
  /**
    * Add two tensors. If possible,
   * use <tt>operator +=</tt> instead
   * since this does not need to
   * copy a point at least once.
   */
  Tensor<0,0>   operator + (const Tensor<0,0> &) const;
  
  /**
    * Subtract two tensors. If
   * possible, use <tt>operator +=</tt>
   * instead since this does not
   * need to copy a point at least
   * once.
   */
  Tensor<0,0>   operator - (const Tensor<0,0> &) const;
  
  /**
    * Tensor with inverted entries.
   */
  Tensor<0,0>   operator - () const;
  
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
  double norm () const;
  
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
  double norm_square () const;
  
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
  void clear ();
  
//  /**
//    * Fill a vector with all tensor elements.
//   *
//   * This function unrolls all
//   * tensor entries into a single,
//   * linearly numbered vector. As
//   * usual in C++, the rightmost
//   * index marches fastest.
//   */
//  void unroll (Vector<double> &result) const;
  
  /**
    * Determine an estimate for
   * the memory consumption (in
                             * bytes) of this
   * object.
   */
  static unsigned int memory_consumption ();
  
  /** @addtogroup Exceptions
    * @{ */
  
  /**
    * Exception
   */
  DeclException1 (ExcDimNotZero,
                  int,
                  << "Dimension must be == 0, but was " << arg1);
  
  //@}
private:
    /**
    * Store the values in a simple variable
     * array.
     */
    double value;
  
#ifdef DEAL_II_TEMPLATE_SPEC_ACCESS_WORKAROUND
public:
#endif
//    /**
//    * Help function for unroll. If
//     * we have detected an access
//     * control bug in the compiler,
//     * this function is declared
//     * public, otherwise private. Do
//     * not attempt to use this
//     * function from outside in any
//     * case, even if it should be
//     * public for your compiler.
//     */
//    void unroll_recursion (Vector<double> &result,
//                           unsigned int   &start_index) const;
  
private:
    /**
    * Make the following classes
     * friends to this class. In
     * principle, it would suffice if
     * otherrank==2, but that is not
     * possible in C++ at present.
     *
     * Also, it would be sufficient
     * to make the function
     * unroll_loops a friend, but
     * that seems to be impossible as
     * well.
     */
    template <int otherrank, int otherdim>  friend class Tensor;
  
//  /**
//    * Point is allowed access to
//   * the coordinates. This is
//   * supposed to improve speed.
//   */
//  friend class Point<dim>;
};


/**
*  Prints the values of this point in the
 *  form <tt>x1 x2 x3 etc</tt>.
 */
template <>
std::ostream & operator << (std::ostream &out, const Tensor<0,0> &p);

#ifndef DOXYGEN

/*------------------------------- Inline functions: Tensor ---------------------------*/


template <>
inline
Tensor<0,0>::Tensor (const bool initialize)
{
  if (initialize)
    this->value = 0;
}



template <>
inline
Tensor<0,0>::Tensor (const array_type &initializer)
{
  value = initializer;
}



template <>
inline
Tensor<0,0>::Tensor (const Tensor<0,0> &p)
{
  value = p.value;
}


template <>
inline
unsigned int Tensor<0,0>::getDim() const
{return 0;}


template <>
inline
unsigned int Tensor<0,0>::getRank() const
{return 0;}



template <>
inline
double Tensor<0,0>::operator [] (const unsigned int index) const
{
  Assert (index==dim, ExcDimensionMismatch (index, 0));
  return value;
}



template <>
inline
double & Tensor<0,0>::operator [] (const unsigned int index)
{
  Assert (index==dim, ExcDimensionMismatch (index, 0));
  return value;
}



template <>
inline
Tensor<0,0> & Tensor<0,0>::operator = (const Tensor<0,0> &p)
{
  value = p.value;
  return *this;
}



template <>
inline
bool Tensor<0,0>::operator == (const Tensor<0,0> &p) const
{
  if (value != p.value)
    return false;
  return true;
}



template <>
inline
bool Tensor<0,0>::operator != (const Tensor<0,0> &p) const
{
  return !((*this) == p);
}



template <>
inline
Tensor<0,0> & Tensor<0,0>::operator += (const Tensor<0,0> &p)
{
  value += p.value;
  return *this;
}


template <>
inline
Tensor<0,0> & Tensor<0,0>::operator += (const double &p)
{
  value += p;
  return *this;
}



template <>
inline
Tensor<0,0> & Tensor<0,0>::operator -= (const Tensor<0,0> &p)
{
  value -= p.value;
  return *this;
}


template <>
inline
Tensor<0,0> & Tensor<0,0>::operator -= (const double &p)
{
  value -= p;
  return *this;
}



template <>
inline
Tensor<0,0> & Tensor<0,0>::operator *= (const double s)
{
  value *= s;
  return *this;
}



template <>
inline
Tensor<0,0> & Tensor<0,0>::operator /= (const double s)
{
  value /= s;
  return *this;
}





template <>
inline
double Tensor<0,0>::operator * (const Tensor<0,0> &p) const
{
  double q = 0;
  q += value * p.value;
  return q;
}



template <>
inline
Tensor<0,0> Tensor<0,0>::operator + (const Tensor<0,0> &p) const
{
  return (Tensor<0,0>(*this) += p);
}



template <>
inline
Tensor<0,0> Tensor<0,0>::operator - (const Tensor<0,0> &p) const
{
  return (Tensor<0,0>(*this) -= p);
}



template <>
inline
Tensor<0,0> Tensor<0,0>::operator - () const
{
  Tensor<0,0> result;
  result.value = -value;
  return result;
}



template <>
inline
double Tensor<0,0>::norm () const
{
  return std::abs (value);
}



template <>
inline
double Tensor<0,0>::norm_square () const
{
  double s = value * value;
  
  return s;
}



template <>
inline
void Tensor<0,0>::clear ()
{
    value = 0;
}



template <>
inline
unsigned int
Tensor<0,0>::memory_consumption ()
{
  return sizeof(Tensor<0,0>);
}



/**
* Output operator for tensors of rank 0. Print the elements
 * consecutively, with a space in between.
 */
template <>
inline
std::ostream & operator << (std::ostream &out, const Tensor<0,0> &p)
{
  out << p[0] << std::endl;
  
  return out;
}






/**
* Multiplication of a tensor of rank 0 with a scalar double from the right.
 */
template <>
inline
Tensor<0,0>
operator * (const Tensor<0,0> &t,
            const double         factor)
{
  Tensor<0,0> tt;
  tt.value = t.value * factor;
  return tt;
}



/**
* Multiplication of a tensor of rank 0 with a scalar double from the left.
 */
template <>
inline
Tensor<0,0>
operator * (const double         factor,
            const Tensor<0,0> &t)
{
  Tensor<0,0> tt;
  tt.value = t.value * factor;
  return tt;
}



/**
* Division of a tensor of rank 0 by a scalar double.
 */
template <>
inline
Tensor<0,0>
operator / (const Tensor<0,0> &t,
            const double         factor)
{
  Tensor<0,0> tt;
  tt.value = t.value / factor;
  return tt;
}

#endif // DOXYGEN

#endif // __fesystem_tensor_zero_h__


