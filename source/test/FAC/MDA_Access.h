/*
 * File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/test/FAC/MDA_Access.h $
 * Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:	$LastChangedRevision: 1871 $
 * Modified:	$LastChangedDate: 2008-01-17 18:06:43 -0800 (Thu, 17 Jan 2008) $
 * Description:	Light-weight array class
 */

#ifndef include_MDA_Access_h
#define include_MDA_Access_h

#include <sys/types.h>
#include <assert.h>
#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif

#include "tbox/Utilities.h"

using namespace std;

/*!
 * @file
 * @brief Provides classes supporting Fortran-style
 * multidimensional array accessing in C++.
 *
 * The classes are written for performance (or at least
 * to not degrade performance), so they are almost all
 * inlined with no run-time toggle-able error checking.
 * It is possible that this approach leads to long
 * compile times and large binaries if you are using
 * a not-so-smart compiler.  In theory though, because
 * these classes are not doing any extraneous computations,
 * it generates codes that are as optimizable as any
 * other code doing similar functions, including Fortran
 * codes.
 *
 * Five classes are defined in this file:
 * -# MDA_IndexRange: a class to define and
 *    manipulate index range objects.
 * -# MDA_OrderRowMajor: class with functions
 *    to compute order-dependent info from the index
 *    range.  This version is for row-major order.
 * -# MDA_OrderColMajor: the column-major
 *    counterpart of MDA_OrderRowMajor.
 * -# MDA_Access: class to allow access to individually
 *    indexed elements of a multidimensional array.
 * -# MDA_AccessConst: the const counterpart of
 *    MDA_Access, to allow read-only access to data.
 *
 * To give the compiler the maximum amount of information
 * with which to perform optimization, always use locally
 * scoped objects.  These classes are very light-weight,
 * so copying them is cheap.
 */


/*!
 * @brief Defines index ranges for multidimensional arrays
 *
 * Defines the abstract index range and methods for setting
 * and accessing it.
 *
 * Nothing is known about the ordering of the array.
 */
template< unsigned short DIM >
class MDA_IndexRange {

public:

  /*!
   * @brief Type for the dimension counter.
   */
  typedef unsigned short dim_t;
  /*!
   * @brief Type for the index.
   */
  typedef int index_t;

  /*!
   * @brief Constructor for setting index data using size and
   * starting points.
   *
   * Any pointers that are NULL are not used.
   * The resulting default settings are:
   * - Array sizes are 0
   * - Array starting indices are 0
   * Since all arguments have default values,
   * this serves at the default constructor.
   *
   * There is another constructor which accepts the first and final
   * indices instead of the sizes and first indices.
   * @b NOTE: the place of the starting points is different than
   * it is for the constructor taking final indices instead of sizes.
   */
inline
MDA_IndexRange(
  /*! Array sizes		*/ const size_t *sz=((size_t*)0) ,
  /*! Array starting indices	*/ const index_t *st=((int*)0) )
{
  dim_t i;
  if ( st )	{ for ( i=0; i<DIM; ++i ) { d_start[i] = st[i]; } }
  else		{ for ( i=0; i<DIM; ++i ) { d_start[i] = 0; } }
  if ( sz )	{ for ( i=0; i<DIM; ++i ) { d_size[i] = sz[i]; } }
  else		{ for ( i=0; i<DIM; ++i ) { d_size[i] = 0; } }
  setDependentData();
  return;
}

  /*!
   * @brief Constructor for setting index data to a range.
   *
   * This version takes two @c index_t* arguments, for the initial
   * and final indices.  It does not support default arguments
   * until after the indices argument.  @b NOTE: the place of the
   * initial indices is different than it is for the constructor
   * taking sizes instead of final indices.
   *
   * If @c si is @c NULL, starting indices are set to 0.
   * If @c sf is @c NULL, sizes are set to zero.
   */
inline
MDA_IndexRange (
  /*! Array of initial indices	*/ const index_t *si ,
  /*! Array of final indices	*/ const index_t *sf )
{
  index_t i;
  if ( si )	{ for ( i=0; i<DIM; ++i ) d_start[i] = si[i]; }
  else		{ for ( i=0; i<DIM; ++i ) d_start[i] = 0; }
  if ( sf )	{ for ( i=0; i<DIM; ++i ) d_size[i] = 1 + sf[i] - d_start[i]; }
  else		{ for ( i=0; i<DIM; ++i ) d_size[i] = 0; }
  setDependentData();
  return;
}

  //@{ @name Functions to set indices
  
  /*!
   * Set size and starting indices.
   */
inline
void setSizeAndStart (
  /*! Array sizes (NULL for no change)		*/ const size_t *sz=((size_t*)0) ,
  /*! Starting indices (NULL for no change)	*/ const index_t *st=((index_t*)0) )
{
  if (sz) for ( dim_t i=0; i<DIM; ++i ) d_size[i] = sz[i];
  if (st) for ( dim_t i=0; i<DIM; ++i ) d_start[i] = st[i];
  setDependentData();
  return;
}

  /*!
   * Set first and final indices (inclusive).
   */
inline
void setInclusiveRange (
  /*! First valid indices (NULL for no change)	*/ const index_t first[DIM] ,
  /*! Final valid indices (NULL for no change)	*/ const index_t final[DIM] )
{
  if (first) for ( dim_t i=0; i<DIM; ++i ) d_start[i] = first[i];
  if (final) for ( dim_t i=0; i<DIM; ++i ) d_size[i] = final[i] - first[i] + 1;
  setDependentData();
  return;
}

  /*!
   * @brief Adjust the dimensions
   *
   * Adjust the first and final indices.
   * Set the dimension to adjust to >= DIM to adjust @em all dimensions.
   * @b Note: The third argument is the increment to the final index
   * and not the size.
   *
   * @b Note: No error checking is done, for example, to make sure that the
   * resulting size is non-negative.
   *
   * @return Adjusted MDA_IndexRange object
   */
inline
const MDA_IndexRange &adjustDim(
  /*! Dimension to adjust	*/ dim_t d ,
  /*! Increment to first index	*/ index_t first ,
  /*! Increment to final index	*/ index_t final )
{
  if ( d >= DIM ) {
    dim_t i;
    for ( i=0; i<DIM; ++i ) {
      d_start[i] += first;
      d_size[i] += final - first;
    }
  }
  else {
    d_start[d] += first;
    d_size[d] += final - first;
  }
  setDependentData();
  return *this;
}

  //@}



  //@{ @name Comparison functions

  /*!
   * @brief Equivalence comparison
   */
inline
bool operator==( const MDA_IndexRange &r ) const {
  dim_t d;
  for ( d=0; d<DIM; ++d ) {
    if ( d_start[d] != r.d_start[d] ) return false;
    if ( d_size[d] != r.d_size[d] ) return false;
  }
  return true;
}

  /*!
   * @brief Inequivalence comparison
   */
inline
bool operator!=( const MDA_IndexRange &r ) const {
  return !( (*this) == r );
}

  //@}



  //@{ @name IO functions
  /*!
   * @brief Output to ostream.
   */
std::ostream &streamPut( std::ostream &os ) const
{
  os << DIM;
  for ( dim_t i=0; i<DIM; ++i )
    os << ' ' << d_start[i] << ' ' << d_size[i];
  return os;
}
  /*!
   * @brief Input from istream.
   */
std::istream &streamGet( std::istream &is )
{
  dim_t dim;
  is >> dim;
  TBOX_ASSERT(dim == DIM);
  for ( dim_t i=0; i<DIM; ++i )
    is >> d_start[i] >> d_size[i];
  return is;
}

friend std::ostream &operator<<( std::ostream &os, const MDA_IndexRange<DIM> &r ) {
  return r.streamPut(os);
}
friend std::istream &operator<<( std::istream &is, MDA_IndexRange<DIM> &r ) {
  return r.streamGet(is);
}
  //@}



  //@{ @name Functions to facilitate looping

  /*!
   * @brief Give starting index of a given dimension.
   */
inline
const index_t &beg( /*! index of dimension */ size_t i ) const {
  return d_start[i];
}

  /*!
   * @brief Give ending index (one more than the last valid index)
   * of a given dimension.
   */
inline
const index_t &end( /*! index of dimension */ size_t i ) const {
  return d_stop[i];
}

  /*!
   * @brief Give size along a given dimension.
   */
inline
const size_t &size( /*! index of dimension */ size_t i ) const {
  return d_size[i];
}

  /*!
   * @brief Give size for all dimensions.
   */
inline
size_t totalSize() const {
  dim_t i=0;
  size_t total_size = d_size[i];
  for ( i=1; i<DIM; ++i ) total_size *= d_size[i];
  return total_size;
}

  //@}



  //@{ @name Error checking functions

  //! Check if indices are in range.
inline
bool inRange( index_t i0 ) const {
  return ( i0 <= d_start[0] )
    && ( i0 < (d_start[0]+d_size[0]) );
}

  //@}





private:


//! Set dependent data.
void setDependentData() {
  index_t i;
  for ( i=0; i<DIM; ++i ) d_stop[i] = d_start[i] + d_size[i] - 1;
}




protected:
  //! @brief Array of starting indices
  index_t d_start[DIM];

  //! @brief Array of stopping indices
  index_t d_stop[DIM];

  //! @brief Array of sizes
  size_t d_size[DIM];


};


/**********************************************************************/
/**********************************************************************/


/*!
 * @brief Performs computations based for row-major arrays.
 *
 * This class has no state.  But its type is important.
 * Its job is to compute things that are dependent on
 * element order in memory, in this case, for the
 * row-major order.
 */
template < unsigned short DIM >
class MDA_OrderRowMajor : private MDA_IndexRange<DIM> {
public:
  typedef int index_t;
  typedef MDA_IndexRange<DIM> range_t;
  typedef typename range_t::dim_t dim_t;
  typedef MDA_OrderRowMajor<DIM-1> reduced_order_t;
  //! @brief Similar to MDA_IndexRange constructor.
inline
MDA_OrderRowMajor (
  /*! Array sizes		*/ const size_t *sz=((size_t*)0) ,
  /*! Array starting indices	*/ const index_t *st=((index_t*)0) )
: MDA_IndexRange<DIM>(sz,st) {
  computeSizeDependentData();
  return;
}
  //! @brief Similar to MDA_IndexRange constructor.
inline
MDA_OrderRowMajor (
  /*! Array of initial indices	*/ const index_t *si ,
  /*! Array of final indices	*/ const index_t *sf )
: MDA_IndexRange<DIM>(si,sf) {
  computeSizeDependentData();
  return;
}

  /*!
   * @brief Constructor for specifying index range object
   */
inline
MDA_OrderRowMajor (
  /*! Array index object	*/ const range_t &r
)
: MDA_IndexRange<DIM>(r) {
  computeSizeDependentData();
  return;
}


  //@{
  //! @name Access to index range
  /*!
   * @brief Const access to the index range object.
   *
   * The index range cannot be modified through this reference.
   * To modify the index range, use other member functions.
   */
inline const range_t &range() const {
  return *this;
}
  //! @brief Similar to MDA_IndexRange::setSizeAndStart().
inline
const MDA_OrderRowMajor &setSizeAndStart (
  const size_t *sz=((size_t*)0) ,
  const index_t *st=((index_t*)0) )
{
  range_t::setSizeAndStart( sz, st );
  computeSizeDependentData();
  return *this;
}
  //! @brief Similar to MDA_IndexRange::setInclusiveRange().
inline
const MDA_OrderRowMajor &setInclusiveRange (
  const index_t first[DIM] ,
  const index_t final[DIM] )
{
  range_t::setInclusiveRange( first, final );
  computeSizeDependentData();
  return *this;
}
  //! @brief Similar to MDA_IndexRange::adjustDim().
inline
const MDA_OrderRowMajor &adjustDim (
  dim_t d ,
  index_t first ,
  index_t final )
{
  range_t::adjustDim(d,first,final);
  computeSizeDependentData();
  return *this;
}
//@}

//@{
//! @name Logical comparisons
  /*!
   * @brief Equivalence comparison.
   *
   * Only independent data is compared, not dependent (redundant) data.
   */
inline
bool operator==(const MDA_OrderRowMajor &r) const {
  return range() == r.range();
}
  /*!
   * @brief Inequivalence comparison.
   *
   * Only independent data is compared, not dependent (redundant) data.
   */
inline
bool operator!=(const MDA_OrderRowMajor &r) const {
  return range() != r.range();
}
//@}

//@{
//! @name Functions to compute offsets
inline index_t fixedOffset () const {
  return d_fixed_offset;
}
inline index_t offset (
                    index_t i0 ) /* SAMRAI_RESTRICT */ const {
  return i0- this -> beg(0);
  /*
  return
   * d_fixed_offset +
   * i0;
   */
}
inline index_t offset (
                    index_t i0 , index_t i1 ) /* SAMRAI_RESTRICT */ const {
  return
    (i0-this -> d_start[0])*d_total_size[1] +
    (i1-this -> d_start[1]);
  /*
   * return (i1-beg(1)) + size(1)*(i0-beg(0));
  return
   * d_fixed_offset +
   * i0*d_total_size[1] +
   * i1;
   */
}
inline index_t offset (
                    index_t i0 , index_t i1 , index_t i2 ) /* SAMRAI_RESTRICT */ const {
  return
    (i0-this -> d_start[0])*d_total_size[1] +
    (i1-this -> d_start[1])*d_total_size[2] +
    (i2-this -> d_start[2]);
  /*
   * return (i2-beg(2)) + size(2)*( (i1-beg(1)) + size(1)*(i0-beg(0)) );
  return
   * d_fixed_offset +
   * i0*d_total_size[1] +
   * i1*d_total_size[2] +
   * i2;
   */
}
inline index_t offset (
                    index_t i0 , index_t i1 , index_t i2 , index_t i3 ) /* SAMRAI_RESTRICT */ const {
  return
    (i0-this -> d_start[0])*d_total_size[1] +
    (i1-this -> d_start[1])*d_total_size[2] +
    (i2-this -> d_start[2])*d_total_size[3] +
    (i3-this -> d_start[3]);
  /*
   * return (i3-beg(3)) + size(3)*( (i2-beg(2)) + size(2)*( (i1-beg(1)) + size(1)*(i0-beg(0)) ) );
  return
   * d_fixed_offset +
   * i0*d_total_size[1] +
   * i1*d_total_size[2] +
   * i2*d_total_size[3] +
   * i3;
   */
}
/*!
 * @brief Compute offsets for arbitrary @c DIM
 *
 * This is flexible but not efficient!
 * You should use dimension-specific offset computations whenever possible.
 */
inline index_t offset ( const index_t (&i)[DIM] ) /* SAMRAI_RESTRICT */ const {
  int d;
  size_t o = i[DIM-1] - this -> beg(DIM-1);
  for ( d=DIM-2; d>=0; --d ) o += ( i[d] - this -> d_start[d] )*d_total_size[i+1];
  return o;
}
inline index_t variableOffset (
                    index_t i0 ) /* SAMRAI_RESTRICT */ const {
  return
    i0;
}
inline index_t variableOffset (
                    index_t i0 , index_t i1 ) /* SAMRAI_RESTRICT */ const {
  return
    i0*d_total_size[1] +
    i1;
}
inline index_t variableOffset (
                    index_t i0 , index_t i1 , index_t i2 ) /* SAMRAI_RESTRICT */ const {
  return
    i0*d_total_size[1] +
    i1*d_total_size[2] +
    i2;
}
inline index_t variableOffset (
                    index_t i0 , index_t i1 , index_t i2 , index_t i3 ) /* SAMRAI_RESTRICT */ const {
  return
    i0*d_total_size[1] +
    i1*d_total_size[2] +
    i2*d_total_size[3] +
    i3;
}
  //@}
/*!
 * @brief Return the total size of subarray starting with dimension d
 */
inline size_t totalSize ( unsigned short d ) const {
  return d_total_size[d];
}
/*!
 * @brief Computes the order object and offset for reducing the slowest
 * dimension.
 *
 * A reduced array is the subarray resulting from fixing the slowest
 * (first) index.  The reduced array has one fewer dimension, a different
 * ordering object and its data starts at a different point in memory.
 * The change in starting point is the returned offset value, and the
 * new order object is returned in the referenced argument.
 *
 * @return Pointer offset (always positive) to the reduced array pointer.
 */
inline size_t reduce( index_t i, reduced_order_t &new_order ) const {
  new_order.setSizeAndStart( &this -> size(1), &this -> beg(1) );
  return offset(i)*d_total_size[1];
}
private:
/*!
 * @brief Recompute the total sizes array, which is dependent on sizes.
 */
inline void computeSizeDependentData() {
  d_total_size[DIM-1] = this -> d_size[DIM-1];
  d_fixed_offset = -this -> d_start[DIM-1];
  int i=DIM-2;
  for ( ; i>0; --i ) {
    d_total_size[i] = this -> d_size[i]*d_total_size[i+1];
    d_fixed_offset -= this -> d_start[i]*d_total_size[i+1];
  }
  return;
}
/*!
 * @brief Total sizes of sub-dimensional arrays.
 *
 * @c d_total_size[i] is the size of the sub-matrix contained in the
 * last (fast) i dimensions of the array.  Incidentally, the stride
 * size of dimension @c i is @c d_total_size[i+1] (and the stride size
 * for dimension @c DIM-1 is 1.
 *
 * @c d_total_size[i] is really equal to
 * @f$ \prod_{j=i}^{DIM-1} size_j @f$
 *
 * This member simply caches size-dependent data.
 */
  size_t d_total_size[DIM];
/*!
 * @brief The fixed portions of offset calculations.
 *
 * Offsets can be separated into a fixed part (dependent only on range)
 * and a variable part (dependent on indices).  To prevent repeated
 * computation of the fixed part, it is saved in this variable.
 * Note that a good optimizing compiler should already do this,
 * so doing it in the code may not really be needed.
 *
 * This member simply caches size-dependent data.
 */
  index_t d_fixed_offset;
};	// end MDA_OrderRowMajor


/**********************************************************************/
/**********************************************************************/


/*!
 * @brief Performs computations based for column-major arrays.
 *
 * This class has no state.  But its type is important.
 * Its job is to compute things that are dependent on
 * element order in memory, in this case, for the
 * column-major order.
 */
template < unsigned short DIM >
class MDA_OrderColMajor : private MDA_IndexRange<DIM> {
public:
  typedef int index_t;
  typedef MDA_IndexRange<DIM> range_t;
  typedef typename range_t::dim_t dim_t;
  typedef MDA_OrderColMajor<DIM-1> reduced_order_t;
  //! @brief Similar to MDA_IndexRange constructor.
inline
MDA_OrderColMajor (
  /*! Array sizes		*/ const size_t *sz=((size_t*)0) ,
  /*! Array starting indices	*/ const index_t *st=((index_t*)0) )
: MDA_IndexRange<DIM>(sz,st) {
  computeSizeDependentData();
  return;
}
  //! @brief Similar to MDA_IndexRange constructor.
inline
MDA_OrderColMajor (
  /*! Array of initial indices	*/ const index_t *si ,
  /*! Array of final indices	*/ const index_t *sf )
: MDA_IndexRange<DIM>(si,sf) {
  computeSizeDependentData();
  return;
}

  /*!
   * @brief Constructor for specifying index range object
   */
inline
MDA_OrderColMajor (
  /*! Array index object	*/ const range_t &r
)
: MDA_IndexRange<DIM>(r) {
  computeSizeDependentData();
  return;
}

  //@{
  //! @name Access to index range (see MDA_IndexRange)
  /*!
   * @brief Const access to the index range object.
   *
   * The index range cannot be modified through this reference.
   * To modify the index range, use other member functions.
   */
inline const range_t &range() const {
  return *this;
}
  //! @brief Similar to MDA_IndexRange::setSizeAndStart().
inline
const MDA_OrderColMajor &setSizeAndStart (
  const size_t *sz=((size_t*)0) ,
  const index_t *st=((index_t*)0) )
{
  range_t::setSizeAndStart( sz, st );
  computeSizeDependentData();
  return *this;
}
  //! @brief Similar to MDA_IndexRange::setInclusiveRange().
inline
const MDA_OrderColMajor &setInclusiveRange (
  const index_t first[DIM] ,
  const index_t final[DIM] )
{
  range_t::setInclusiveRange( first, final );
  computeSizeDependentData();
  return *this;
}
  //! @brief Similar to MDA_IndexRange::adjustDim().
inline
const MDA_OrderColMajor &adjustDim (
  dim_t d ,
  index_t first ,
  index_t final )
{
  range_t::adjustDim(d,first,final);
  computeSizeDependentData();
  return *this;
}
//@}

//@{
//! @name Logical comparisons
  /*!
   * @brief Equivalence comparison.
   *
   * Only independent data is compared, not dependent (redundant) data.
   */
inline
bool operator==(const MDA_OrderColMajor &r) const {
  return range() == r.range();
}
  /*!
   * @brief Inequivalence comparison.
   *
   * Only independent data is compared, not dependent (redundant) data.
   */
inline
bool operator!=(const MDA_OrderColMajor &r) const {
  return range() != r.range();
}
//@}

//@{
//! @name Functions to compute offsets
inline index_t fixedOffset () const {
  return d_fixed_offset;
}
inline index_t offset (
                    index_t i0 ) /* SAMRAI_RESTRICT */ const {
  return i0-this->beg(0);
  /*
  return
   * d_fixed_offset +
   * i0;
   */
}
inline index_t offset (
                    index_t i0 , index_t i1 ) /* SAMRAI_RESTRICT */ const {
  return
    (i0-this -> d_start[0]) +
    (i1-this -> d_start[1])*d_total_size[0];
  /*
   * return (i0-this->beg(0)) + size(0)*(i1-this->beg(1));
  return
   * d_fixed_offset +
   * i0 +
   * i1*d_total_size[0];
   */
}
inline index_t offset (
                    index_t i0 , index_t i1 , index_t i2 ) /* SAMRAI_RESTRICT */ const {
  return
    (i0-this -> d_start[0]) +
    (i1-this -> d_start[1])*d_total_size[0] +
    (i2-this -> d_start[2])*d_total_size[1] ;
  /*
   * return (i0-this->beg(0)) + size(0)*( (i1-this->beg(1)) + size(1)*(i2-this->beg(2)) );
  return
   * d_fixed_offset +
   * i0 +
   * i1*d_total_size[0] +
   * i2*d_total_size[1] ;
   */
}
inline index_t offset (
                    index_t i0 , index_t i1 , index_t i2 , index_t i3 ) /* SAMRAI_RESTRICT */ const {
  return
    (i0-this -> d_start[0]) +
    (i1-this -> d_start[1])*d_total_size[0] +
    (i2-this -> d_start[2])*d_total_size[1] +
    (i3-this -> d_start[3])*d_total_size[2] ;
  /*
   * return (i0-this->beg(0)) + size(0)*( (i1-this->beg(1)) + size(1)*( (i2-this->beg(2)) + size(2)*(i3-this->beg(3)) ) );
  return
   * d_fixed_offset +
   * i0 +
   * i1*d_total_size[0] +
   * i2*d_total_size[1] +
   * i3*d_total_size[2] ;
   */
}
/*!
 * @brief Compute offsets for arbitrary @c DIM
 *
 * This is flexible but not efficient!
 * You should use dimension-specific offset computations whenever possible.
 */
inline index_t offset ( const index_t (&i)[DIM] ) /* SAMRAI_RESTRICT */ const {
  int d;
  size_t o = i[0] - this->beg(0);
  for ( d=1; d<DIM; ++d ) o += ( i[d] - this -> d_start[d] )*d_total_size[i-1];
  return o;
}
inline index_t variableOffset (
                    index_t i0 ) /* SAMRAI_RESTRICT */ const {
  return
    i0;
}
inline index_t variableOffset (
                    index_t i0 , index_t i1 ) /* SAMRAI_RESTRICT */ const {
  return
    i0 +
    i1*d_total_size[0];
}
inline index_t variableOffset (
                    index_t i0 , index_t i1 , index_t i2 ) /* SAMRAI_RESTRICT */ const {
  return
    i0 +
    i1*d_total_size[0] +
    i2*d_total_size[1] ;
}
inline index_t variableOffset (
                    index_t i0 , index_t i1 , index_t i2 , index_t i3 ) /* SAMRAI_RESTRICT */ const {
  return
    i0 +
    i1*d_total_size[0] +
    i2*d_total_size[1] +
    i3*d_total_size[2] ;
}
//@}
/*!
 * @brief Return the total size of subarray starting with dimension d
 */
inline size_t totalSize ( unsigned short d ) const {
  return d_total_size[d];
}
/*!
 * @brief Computes the order object and offset for reducing the slowest
 * dimension.
 *
 * A reduced array is the subarray resulting from fixing the slowest
 * (last) index.  The reduced array has one fewer dimension, a different
 * ordering object and its data starts at a different point in memory.
 * The change in starting point is the returned offset value, and the
 * new order object is returned in the referenced argument.
 *
 * @return Pointer offset (always positive) to the reduced array pointer.
 */
inline size_t reduce( index_t i, reduced_order_t &new_order ) const {
  new_order.setSizeAndStart( &this->size(0), &this->beg(0) );
  return offset(i)*d_total_size[DIM-2];
}
private:
/*!
 * @brief Recompute the total sizes array, which is dependent on sizes.
 */
inline void computeSizeDependentData() {
  d_total_size[0] = this -> d_size[0];
  d_fixed_offset = -this -> d_start[0];
  int i=0;
  for ( ; i<DIM-1; ++i ) {
    d_total_size[i+1] = this -> d_size[i+1]*d_total_size[i];
    d_fixed_offset -= this -> d_start[i+1]*d_total_size[i+1];
  }
  return;
}
/*!
 * @brief Total sizes of sub-dimensional arrays.
 *
 * @c d_total_size[i] is the size of the sub-matrix contained in the
 * first (fast) i+1 dimensions of the array.  Incidentally, the stride
 * size of dimension @c i is @c d_total_size[i-1] (and the stride size
 * for dimension @c 0 is 1.
 *
 * @c d_total_size[i] is really equal to
 * @f$ \prod_{j=0}^{i} size_j @f$
 *
 * This member simply caches size-dependent data.
 */
  size_t d_total_size[DIM];
/*!
 * @brief The fixed portions of offset calculations.
 *
 * Offsets can be separated into a fixed part (dependent only on range)
 * and a variable part (dependent on indices).  To prevent repeated
 * computation of the fixed part, it is saved in this variable.
 * Note that a good optimizing compiler should already do this,
 * so doing it in the code may not really be needed.
 *
 * This member simply caches size-dependent data.
 */
  index_t d_fixed_offset;
};	// end MDA_OrderColMajor


/**********************************************************************/
/**********************************************************************/


/*!
 * @brief Non-const multidimensional array access.
 *
 * This class @em never allocates or deallocates data.
 * It takes pointers to preallocated data and prodes an
 * interface to that data.  Member functions are used
 * to give that interface.
 *
 * This class provides functions for explicit index checking,
 * but it does @em NO implicit error checking on either the
 * dimensionality of the array or it size.
 * Such may be done through subclassing.
 *
 * The member functions should all be inlined for better
 * performance.
 *
 * This template class is set up to work with either
 * row-major or column-major data, depending on the
 * third template argument, which should be one of
 * -# @c MDA_OrderRowMajor<DIM>
 * -# @c MDA_OrderColMajor<DIM>
 *
 * The reduce() function return a new array of smaller
 * dimensional that require less integer arithmetic
 * to access individual array members.  This should help
 * in optimizing code.  (My preliminary performance tests
 * using gcc and gprof on i686 Linux showed that the
 * MDA_Access functions run at half to slightly
 * faster than the speed of Fortran, depending on use of
 * reduced arrays.  However, note that gcc is not great at
 * optimizing Fortran.)
 */
template < class TYPE, unsigned short DIM, class OrderType=MDA_OrderRowMajor<DIM> >
class MDA_Access {

public:

  /*!
   * @brief Type of data.
   */
  typedef TYPE value_t;
  typedef MDA_IndexRange<DIM> range_t;
  typedef typename range_t::dim_t dim_t;
  typedef typename range_t::index_t index_t;
  typedef OrderType order_t;
  typedef typename OrderType::reduced_order_t reduced_order_t;

  /*!
   * @brief Constructor for setting all data, with default values.
   *
   * Any pointers that are NULL are not used.  The resulting default
   * settings are:
   * - Data pointer is NULL
   * - Array sizes are 0
   * - Array starting indices are 0
   *
   * There is another constructor which accepts the first and final
   * indices instead of the sizes and first indices.
   * @b NOTE: the place of the initial indices is different than
   * it is for the constructor taking final indices instead of sizes.
   */
inline
MDA_Access (
  /*! Pointer to data		*/ value_t *p=((value_t*)0) ,
  /*! Array sizes		*/ const size_t *sz=((size_t*)0) ,
  /*! Array starting indices	*/ const index_t *st=((index_t*)0) )
: d_ptr(p),
  d_order(sz,st) {
  setPtr1();
  return;
}

  /*!
   * @brief Constructor for setting all data, with default values.
   *
   * Any pointers that are NULL are not used.
   * The resulting default settings are:
   * - Data pointer is NULL
   * - Array sizes are 0
   * - Array starting indices are 0
   *
   * This version takes two @c int* arguments, for the initial
   * and final indices.  It does not support default arguments
   * until after the indices argument.  @b NOTE: the place of the
   * initial indices is different than it is for the constructor
   * taking sizes instead of final indices.
   *
   * If @c si is @c NULL, starting indices are set to 0.
   * If @c sf is @c NULL, sizes are set to zero.
   */
inline
MDA_Access (
  /*! Pointer to data		*/ value_t *p ,
  /*! Array of initial indices	*/ const index_t *si ,
  /*! Array of final indices	*/ const index_t *sf )
: d_ptr(p),
  d_order(si,sf) {
  setPtr1();
  return;
}

  /*!
   * @brief Constructor for specifying pointer and ordering object
   */
inline
MDA_Access (
  /*! Pointer to data		*/ value_t *p ,
  /*! Array index object	*/ const order_t &r
)
: d_ptr(p),
  d_order(r) {
  setPtr1();
  return;
}

  /*!
   * @brief Copy constructor
   */
inline
MDA_Access (
  /*! Copyee object		*/ const MDA_Access &r
)
: d_ptr(r.d_ptr),
  d_order(r.d_order) {
  setPtr1();
  return;
}

  /*!
   * @brief Conversion into boolean.
   *
   * @return true iff data pointer is not NULL.
   */
inline
operator bool() const {
  return ( d_ptr != (value_t*)0 );
}

  /*!
   * @brief Conversion into pointer.
   *
   * @return the data pointer.
   */
inline
operator value_t*() const {
  return d_ptr;
}




  /*!
   * @brief Set the data pointer.
   */
inline
void setPointer (
  /*! Pointer value	*/ value_t *p
) {
  d_ptr = p;
  setPtr1();
}

  /*!
   * Set size and starting indices.
   *
   * @see MDA_IndexRange
   */
inline
void setSizeAndStart (
  /*! Array sizes (NULL for no change)		*/ const size_t *sz=((size_t*)0) ,
  /*! Starting indices (NULL for no change)	*/ const index_t *st=((index_t*)0) )
{
  d_order.setSizeAndStart( sz, st );
  setPtr1();
}

  /*!
   * Set first and final indices (inclusive).
   *
   * @see MDA_IndexRange
   */
inline
void setInclusiveRange (
  /*! First valid indices (NULL for no change)	*/ const index_t first[DIM] ,
  /*! Final valid indices (NULL for no change)	*/ const index_t final[DIM] )
{
  d_order.setInclusiveRange( first, final );
  setPtr1();
}

  /*!
   * @brief Adjust the dimensions
   *
   * @see MDA_IndexRange::adjustDim.
   */
inline
const range_t &adjustDim (
  /*! Dimension to adjust	*/ dim_t d ,
  /*! Increment to first index	*/ index_t first ,
  /*! Increment to final index	*/ index_t final )
{
  d_order.adjustDim(d,first,final);
  setPtr1();
  return d_order.range();
}



//@{ @name Comparison functions

  /*!
   * @name Equivalence comparison
   */
inline
bool operator==( const MDA_Access &r ) const {
  if ( d_order != r.d_order ) return false;
  if ( d_ptr != r.d_ptr ) return false;
  return true;
}

  /*!
   * @name Inequivalence comparison
   */
inline
bool operator!=( const MDA_Access &r ) const {
  return !( (*this) == r );
}

//@}



//@{ @name Functions for accessing items

inline
const range_t &range() const { return d_order.range(); };
inline
const index_t &beg(size_t i) const { return d_order.range().beg(i); }
inline
const index_t &end(size_t i) const { return d_order.range().end(i); }
inline
const size_t &size(size_t i) const { return d_order.range().size(i); }

  /*!
   * @brief Grant general access to item in a 1D array.
   */
inline
value_t &operator() ( index_t i0 ) /* SAMRAI_RESTRICT */ const {
  return d_ptr[d_order.offset(i0)];
  /*
  return d_ptr1[i0];
   */
}

  /*!
   * @brief Grant general access to item in a 2D array.
   */
inline
value_t &operator() ( index_t i0 , index_t i1 ) /* SAMRAI_RESTRICT */ const {
  return d_ptr[d_order.offset(i0,i1)];
  // return d_ptr1[d_order.variableOffset(i0,i1)];
}

  /*!
   * @brief Grant general access to item in a 3D array.
   */
inline
value_t &operator() ( index_t i0 , index_t i1 , index_t i2 ) /* SAMRAI_RESTRICT */ const {
  return d_ptr[d_order.offset(i0,i1,i2)];
  // return d_ptr1[d_order.variableOffset(i0,i1,i2)];
}

  /*!
   * @brief Grant general access to item in a 4D array.
   */
inline
value_t &operator() ( index_t i0 , index_t i1 , index_t i2 , index_t i3 ) /* SAMRAI_RESTRICT */ const {
  return d_ptr[d_order.offset(i0,i1,i2,i3)];
  // return d_ptr1[d_order.variableOffset(i0,i1,i2,i3)];
}

  /*!
   * @brief Special case for 1D arrays, identical to @c operator(index_t),
   * using pre-added fixed offsets.
   *
   * This @em may be more efficient than @c (i) but it only works in 1D.
   * It is not guaranteed to work if the fixed offset is negative and
   * has greater value than the pointer address, since the addition of
   * the two gives a negative address, which the C standard leaves as
   * undefined behavior.
   */
inline
value_t &operator[] ( index_t i0 ) /* SAMRAI_RESTRICT */ const {
  return d_ptr1[i0];
  /*
  return d_ptr[d_order.offset(i0)];
   */
}

//@}




//@{ @name Functions to extract reduced-dimensional arrays.

/*!
 * @brief Fix the index of the slowest dimension and return
 * the corresponding sub-array.
 *
 * This function is meant to facilitate optimization when using
 * this class.  In nested loops, the inner loops executes many
 * times with the indices corresponding to outer loops remaining
 * constants.  This leads to many many repeated integer arithmetics
 * operations that could be removed from the inner loop (but may
 * not be removed automatically by the compiler optimization).
 * To do this, reduce the array dimensionality one dimension at a time,
 * by fixing index corresponding to the slowest varying dimension.
 * (If you code is written to maximize cache data, this will NOT
 * be the index of the innermost loop.)
 *
 * To reduce multiple dimensions, string these calls together,
 * i.e. @c array.reduce(i).reduce(j).  However, since reduction
 * contains loops that cost O(DIM) and may be difficult for
 * compilers to optimize, you may want to save @c array.reduce(i)
 * and reuse it.
 *
 * @param i hier::Index in slowest dimension, which is the first
 * dimension in a row-major array and the last dimension in a
 * column-major array.
 *
 * @return The sub-array of dimension @c DIM-1, corresponding to
 * the index given.
 */
inline
MDA_Access<TYPE,DIM-1,typename OrderType::reduced_order_t> reduce( index_t i ) const {
  typename OrderType::reduced_order_t new_order;
  int ptr_offset;
  ptr_offset = d_order.reduce( i, new_order );
  return MDA_Access<TYPE,DIM-1,typename OrderType::reduced_order_t>(
    d_ptr + ptr_offset , new_order );
}

//@}




//! Pointer to data.
private: value_t * /* SAMRAI_RESTRICT */ d_ptr;
/*!
 * @brief Value of @c d_ptr-beg(0), used for optimizing 1D access.
 *
 * The use of precomputed @c d_ptr1=d_ptr-beg(0) speeds up 1D offset
 * computations by allowing us to compute  @c d_ptr1+i0 instead of
 * @c d_ptr+(i0-beg(0)), saving one integer subtraction for each
 * 1D access.  However, this could be a real problem if @c d_ptr<beg(0).
 * So far, that has not happened, and we are keeping our fingers
 * crossed.
 *
 * @see setPtr1()
 */
private: value_t * /* SAMRAI_RESTRICT */ d_ptr1;
private: void setPtr1() {
  /*
   * If the following assert fails, our d_ptr1 optimization may
   * give undefined result.
  TBOX_ASSERT( d_order.fixedOffset() > 0 ||
	  (unsigned long)d_ptr > (unsigned long)(-d_order.fixedOffset()) );
   */
  d_ptr1 = d_ptr + d_order.fixedOffset();
}
//! Offset computing object
private: order_t d_order;


};	// class MDA_Access


/**********************************************************************/
/**********************************************************************/


/*!
 * @brief Const data version of the multidimensional array access
 * template class MDA_Access.
 *
 * This class is almost exactly identical to its non-const
 * counterpart, MDA_Access.  It is used when the data
 * is const.
 *
 * This class differs only in that the value type is a const.
 * In fact, this class is trivial,
 * except for the public inheritance of MDA_Access
 * with the const type for the first template argument,
 * a constructor to build an object from a MDA_Access
 * object and an assignment operator to assign from a
 * MDA_Access object.
 * Other than that, see MDA_Access for documentations.
 *
 * The interfaces that are added by this class are trivial,
 * mirroring the interfaces defined in MDA_Access
 * with minor changes.
 *
 * @see MDA_Access
 */
template < class TYPE, unsigned short DIM, class OrderType=MDA_OrderRowMajor<DIM> >
class MDA_AccessConst : public MDA_Access< const TYPE, DIM, OrderType > {
public:

  /*!
   * @brief Type of data.
   *
   * This declaration is redundant because it should already be inherited,
   * but the xlC compiler on ASCI Blue does not get it.
   */
  typedef const TYPE value_t;
  typedef MDA_IndexRange<DIM> range_t;
  typedef typename range_t::dim_t dim_t;
  typedef typename range_t::index_t index_t;
  typedef OrderType order_t;

  /*!
   * @brief See the MDA_Access version of this function.
   * @see MDA_Access::MDA_Access(value_t*,const size_t*,const index_t*)
   */
inline
MDA_AccessConst (
  /*! Pointer to data		*/ value_t *p=((TYPE*)0) ,
  /*! Array sizes		*/ const size_t *sz=((size_t*)0) ,
  /*! Array starting indices	*/ const index_t *st=((index_t*)0) )
 : MDA_Access< const TYPE, DIM, OrderType > (p,sz,st)
{
  return;
}
  /*!
   * @brief See the MDA_Access version of this function.
   * @see MDA_Access::MDA_Access(value_t*,const index_t*,const index_t*)
   */
inline
MDA_AccessConst (
  /*! Pointer to data		*/ value_t *p ,
  /*! Array of initial indices	*/ const index_t *si ,
  /*! Array of final indices	*/ const index_t *sf )
: MDA_Access< const TYPE, DIM, OrderType > (p,si,sf)
{
  return;
}
  /*!
   * @brief See the MDA_Access version of this function.
   * @see MDA_Access::MDA_Access(value_t*,const MDA_IndexRange<DIM>&)
   */
inline
MDA_AccessConst (
  /*! Pointer to data		*/ value_t *p ,
  /*! Array index object	*/ const MDA_IndexRange<DIM> &r
)
: MDA_Access< const TYPE, DIM, OrderType > (p,r)
{
  return;
}
  /*!
   * @brief Construct from an object of the non-const version.
   * @see MDA_Access::MDA_Access(const MDA_Access<const TYPE,DIM>&)
   */
inline
MDA_AccessConst (
  const MDA_Access<TYPE,DIM,OrderType> &r
)
: MDA_Access< const TYPE, DIM,OrderType > ( (const TYPE*)(r), r.range() )
{
  return;
}
  /*!
   * @brief Assign value from an object of the non-const version.
   */
inline
const MDA_AccessConst &operator=( const MDA_Access<TYPE,DIM,OrderType> &r ) {
  (MDA_Access<TYPE,DIM,OrderType>&)(*this) = r;
  return *this;
}
};





#endif	// include_MDA_Access_h
