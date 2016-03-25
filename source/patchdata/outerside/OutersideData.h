//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outerside/OutersideData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated outerside centered patch data type
//

#ifndef included_pdat_OutersideData
#define included_pdat_OutersideData

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchData.h"
#include "ArrayData.h"
#include "SideData.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "tbox/Arena.h"
#include "tbox/Complex.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {


/*!
 * @brief Class OutersideData<DIM> provides an implementation for data defined
 * at cell sides (faces) on the boundaries of AMR patches.  It is derived from 
 * the hier::PatchData interface common to all SAMRAI patch data types.  Given
 * a CELL-centered AMR index space box, an outerside data object represents
 * data of some template TYPE and depth on the cell sides (faces) on the boundary
 * of the box.  Here, depth indicates the number of data values at each face
 * index location.  The OuteredgsideGeometry class provides the translation
 * between the standard SAMRAI cell-centered AMR index space and
 * outerside-centered data.
 *
 * Outerside data is stored in 2*DIM arrays, each of which contains data
 * associated with side (face) indices normal to a coordinate axis direction 
 * and an upper or lower box side (face) in the face normal direction.  
 * The data layout in the outerside data arrays matches the corresponding array 
 * sections provided by the side data implementation.  Also, in each of array, 
 * memory allocation is in column-major ordering (e.g., Fortran style) so that 
 * the leftmost index runs fastest in memory.  For example, a three-dimensional 
 * outerside data object created over a CELL-centered AMR index space 
 * [l0:u0,l1:u1,l2:u2] allocates six data arrays dimensioned as follows:
 * \verbatim

   face normal 0:
     lower face      [ l0:l0     , l1:u1     , l2:u2     , d ]
     upper face      [ u0+1:u0+1 , l1:u1     , l2:u2     , d ]

   face normal 1:
     lower face      [ l0:u0     , l1:l1     , l2:u2     , d ]
     upper face      [ l0:u0     , u1+1:u1+1 , l2:u2     , d ]

   face normal 2:
     lower face      [ l0:u0     , l1:u1     , l2:l2     , d ]
     upper face      [ l0:u0     , l1:u1     , u2+1:u2+1 , d ]

 * \endverbatim
 * Here the face normal directions 0, 1, 2 can be thought of as the x, y, and z
 * face normal directions, respectively, and d is the depth index (i.e., number
 * of values at each face index location).  Other spatial dimensions are
 * represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * IMPORTANT: The OuterfaceData<DIM> class provides the same storage
 * as this outerside data class, except that the coordinate directions of the
 * individual arrays are permuted; i.e., OuterfaceData is consistent
 * with the FaceData implementation.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::OutersideDataFactory
 * @see pdat::OutersideGeometry
 * @see pdat::SideIterator
 * @see pdat::SideIndex
 */

template<int DIM, class TYPE>
class OutersideData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outerside-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerside data object will be created.
    *            Note: the ghost cell width is assumed to be zero.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box,
                               int depth);

   /*!
    * @brief Constructor for an outerside data object.
    *
    * Note: Outerside data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerside data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   OutersideData(const hier::Box<DIM>& box,
                 int depth,
                 tbox::Pointer<tbox::Arena> pool =
                    tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief Virtual destructor for a outerside data object.
    */
   virtual ~OutersideData<DIM,TYPE>();

   /*!
    * @brief Return the depth (i.e., the number of data values for
    * each spatial location) of the array.
    */
   int getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular
    * side normal, side, and depth component of the outerside centered
    * array.
    *
    * @param side_normal  integer side normal direction for data,
    *              must satisfy 0 <= side_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   TYPE* getPointer(int side_normal,
                    int side,
                    int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular
    * side normal, side, and depth component of the outerside centered
    * array.
    *
    * @param side_normal  integer side normal direction for data,
    *              must satisfy 0 <= side_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   const TYPE* getPointer(int side_normal,
                          int side,
                          int depth = 0) const;

   /*!
    * @brief Return a reference to data entry corresponding
    * to a given side index, side location, and depth.
    *
    * @param i const reference to SideIndex, @em MUST be
    *          an index on the outerside of the box.
    * @param side  integer (lower/upper location of outerside data),
    *              must satisfy 0 <= side <= 1
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   TYPE& operator()(const SideIndex<DIM>& i,
                    int side,
                    int depth = 0);

   /*!
    * @brief Return a const reference to data entry corresponding
    * to a given side index, side location, and depth.
    *
    * @param i const reference to SideIndex, @em MUST be
    *          an index on the outerside of the box.
    * @param side  integer (lower/upper location of outerside data),
    *              must satisfy 0 <= side <= 1
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   const TYPE& operator()(const SideIndex<DIM>& i,
                          int side,
                          int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for
    * side normal and side location of the outerside centered array.
    *
    * @param side_normal  integer side normal direction for data,
    *              must satisfy 0 <= side_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    */
   ArrayData<DIM,TYPE>& getArrayData(int side_normal,
                                     int side);

   /*!
    * @brief Return a const reference to the array data object for
    * side normal and side location of the outerside centered array.
    *
    * @param side_normal  integer side normal direction for data,
    *              must satisfy 0 <= side_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    */
   const ArrayData<DIM,TYPE>& getArrayData(int side_normal,
                                           int side) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * SideData the same DIM and TYPE.  If not, then an unrecoverable error
    * results.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * SideData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * IMPORTANT: this routine is @b not @b yet @b implemented!
    */
   virtual void copy(const hier::PatchData<DIM>& src,
                     const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be SideData of the same DIM
    * and TYPE and the overlap must be a SideOverlap of the same
    * DIM.  If not, then an unrecoverable error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst,
                      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes side and outerside data objects are
    * defined over the same box) from the given side source data object to
    * this destination outerside data object at the specified depths.
    */
   void copyDepth(int dst_depth,
                  const SideData<DIM,TYPE>& src,
                  int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes side and outerside data objects are
    * defined over the same box) to the given side destination data object
    * from this source outerside data object at the specified depths.
    */
   void copyDepth2(int dst_depth,
                   SideData<DIM,TYPE>& dst,
                   int src_depth) const;

   /*!
    * @brief Return true if the patch data object can estimate the
    * stream size required to fit its data using only index
    * space information (i.e., a box).
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /*!
    * @brief Return the number of bytes needed to stream the data
    * in this patch data object lying in the specified box overlap
    * region.
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * SideOverlap of the same DIM.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * SideOverlap of the same DIM.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Fill all values at depth d with the value t.
    */
   void fill(const TYPE& t,
             int d = 0);

   /*!
    * @brief Fill all values at depth d within the box with the value t.
    */
   void fill(const TYPE& t,
             const hier::Box<DIM>& box,
             int d = 0);

   /*!
    * @brief Fill all depth components with value t.
    */
   void fillAll(const TYPE& t);

   /*!
    * @brief Fill all depth components within the box with value t.
    */
   void fillAll(const TYPE& t,
                const hier::Box<DIM>& box);

   /*!
    * @brief Print all outerside data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void print(const hier::Box<DIM>& box,
              std::ostream& os = tbox::plog,
              int prec = 12) const;

   /*!
    * @brief Print all outerside data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void print(const hier::Box<DIM>& box,
              int depth,
              std::ostream& os = tbox::plog,
              int prec = 12) const;

   /*!
    * @brief Print all outerside centered data values for specified
    * side_normal and side location residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.
    *
    * @param side_normal  integer side normal direction for data,
    *              must satisfy 0 <= side_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void printAxisSide(int side_normal,
                      int side,
                      const hier::Box<DIM>& box,
                      std::ostream& os = tbox::plog,
                      int prec = 12) const;

  /*!
    * @brief Print all outerside centered data values for specified
    * side_normal, side location, and depth residing in the specified box.
    *
    * @param side_normal  integer side normal direction for data,
    *              must satisfy 0 <= side_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void printAxisSide(int side_normal,
                      int side,
                      const hier::Box<DIM>& box,
                      int depth,
                      std::ostream& os = tbox::plog,
                      int prec = 12) const;

   /*!
    * @brief Check that class version and restart file version are equal.
    * If so, read data members from the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
        tbox::Pointer<tbox::Database> database);

   /*!
    * @brief Write out the class version number and other data members to
    * the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void putSpecializedToDatabase(
        tbox::Pointer<tbox::Database> database);

private:
   OutersideData(const OutersideData<DIM,TYPE>&); // not implemented
   void operator=(const OutersideData<DIM,TYPE>&);	  // not implemented

   int d_depth;
   ArrayData<DIM,TYPE> d_data[DIM][2];

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OutersideData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OutersideData.C"
#endif
