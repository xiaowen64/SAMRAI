//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/side/SideData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated side centered patch data type
//

#ifndef included_pdat_SideData
#define included_pdat_SideData

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_pdat_ArrayData
#include "ArrayData.h"
#endif
#ifndef included_pdat_SideIndex
#include "SideIndex.h"
#endif
#ifndef included_pdat_SideIterator
#include "SideIterator.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class SideData<DIM> provides an implementation for data defined
 * at cell sides (faces) on AMR patches.  It is derived from the hier::PatchData
 * interface common to all SAMRAI patch data types.  Given a CELL-centered
 * AMR index space box, a side data object represents data of some template
 * TYPE and depth on the sides (faces) of the cells in the box.  Here, depth
 * indicates the number of data values at each side index location.  The
 * SideGeometry class provides the translation between the standard SAMRAI
 * cell-centered AMR index space and side-centered data.
 *
 * IMPORTANT: The FaceData<DIM> class provides the same storage
 * as this side data class, except that the coordinate directions of the
 * individual arrays are permuted in the face data implementation.
 *
 * Side data is stored in DIM arrays, each of which contains the
 * data for the sides normal to a corresponding coordinate direction.
 * Memory allocation is in column-major ordering (e.g., Fortran
 * style) so that the leftmost index runs fastest in memory.
 * For example, a three-dimensional side data object created over a
 * CELL-centered AMR index space [l0:u0,l1:u1,l2:u2] allocates three data
 * arrays dimensioned as follows:
 * \verbatim

   side normal 0
     [ l0 : u0+1 ,
       l1 : u1 ,
       l2 : u2 , d ]   ,

   side normal 1
     [ l0 : u0 ,
       l1 : u1+1 ,
       l2 : u2 , d ]   ,

   side normal 2
     [ l0 : u0 ,
       l1 : u1 ,
       l2 : u2+1 , d ]   ,
      
 * \endverbatim
 * Here the side normal directions 0, 1, 2 can be thought of as the x, y, and z
 * side normal directions, respectively, and d is the depth index (i.e., number
 * of values at each side index location).  Other spatial dimensions are
 * represented similarly.
 *
 * Note also that it is possible to create a side data object for managing 
 * data at cell sides associated with a single coordinate direction only.
 * See the constructor for more information.  All operations are defined
 * only for the case where data storage is alike between two side data objects.  
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::SideDataFactory
 * @see pdat::SideIndex
 * @see pdat::SideIterator
 * @see pdat::SideGeometry
 */

template<int DIM, class TYPE>
class SideData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent side-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            side data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *               of the ghost cell region around the box over which
    *               the side data will be allocated.
    * @param directions const IntVector reference indicating which
    *                   coordinate directions are assumed to have data
    *                   for the purposes of the calculation.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box,
                               int depth,
                               const hier::IntVector<DIM>& ghosts,
                               const hier::IntVector<DIM>& directions);

   /*!
    * @brief The constructor for a side data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            side data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *               of the ghost cell region around the box over which
    *               the side data will be allocated.
    * @param directions const IntVector reference indicating which
    *                   coordinate directions will have data associated 
    *                   with them. The default state is to create data 
    *                   for all coordinate directions.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   SideData(const hier::Box<DIM>& box,
            int depth,
            const hier::IntVector<DIM>& ghosts,
            const hier::IntVector<DIM>& directions = hier::IntVector<DIM>(1), 
            tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief The virtual destructor for a side data object.
    */
   virtual ~SideData<DIM,TYPE>();

   /*!
    * @brief Return constant reference to vector describing which coordinate
    * directions have data associated with this side data object.
    * 
    * A vector entry of zero indicates that there is no data array
    * allocated for the corresponding coordinate direction.  A non-zero
    * value indicates that a valid data array is maintained for that
    * coordinate direction.
    */
   const hier::IntVector<DIM>& getDirectionVector() const;

   /*!
    * @brief Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular side normal and
    * depth component of the side centered array.
    */
   TYPE* getPointer(int side_normal, int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular side normal
    * and depth component of the side centered array.
    */
   const TYPE* getPointer(int side_normal, int depth = 0) const;

   /*!
    * @brief Return a reference to the data entry corresponding
    * to a given side index and depth.
    */
   TYPE& operator()(const SideIndex<DIM>& i, int depth = 0);

   /*!
    * @brief Return a const reference to the data entry corresponding
    * to a given side index and depth.
    */
   const TYPE& operator()(const SideIndex<DIM>& i, int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for the
    * given side normal of the side centered data object.
    */
   ArrayData<DIM,TYPE>& getArrayData(int side_normal);

   /*!
    * @brief Return a const reference to the array data object for the
    * given side normal of the side centered data object.
    */
   const ArrayData<DIM,TYPE>& getArrayData(int side_normal) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * an SideData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * an SideData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be SideData of the same DIM and TYPE
    * and the overlap must be a SideOverlap of the same DIM. If not,
    * then an unrecoverable error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src,
                     const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be SideData of the same DIM and TYPE
    * and the overlap must be a SideOverlap of the same DIM.  If not,
    * then an unrecoverable error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst,
                      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given CELL-centered AMR index box.
    */
   void copyOnBox(const SideData<DIM,TYPE>& src,
                  const hier::Box<DIM>& box);

   /*!
    * @brief Fast copy (i.e., source and this side data objects are
    * defined over the same box) to this destination side data object
    * from the given source side data object at the specified depths.
    */
   void copyDepth(int dst_depth,
                  const SideData<DIM,TYPE>& src,
                  int src_depth);

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
    * @brief Print all side data values residing in the specified box.
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
    * @brief Print all side data values at the given array depth in
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
    * @brief Print all side centered data values for specified side normal
    * direction residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param side_normal  integer side normal coordinate direction,
    *              must satisfy 0 <= side_normal < DIM
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
   void printAxis(int side_normal,
                  const hier::Box<DIM>& box,
                  std::ostream& os = tbox::plog,
                  int prec = 12) const;

   /*!
    * @brief Print all side centered data values for specified side normal
    * direction residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param side_normal  integer side normal coordinate direction,
    *              must satisfy 0 <= side_normal < DIM
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
   void printAxis(int side_normal,
                  const hier::Box<DIM>& box,
                  int depth,
                  std::ostream& os = tbox::plog,
                  int prec = 12) const;

   /*!
    * Check that class version and restart file version are equal.  If so,
    * read data members from the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
        tbox::Pointer<tbox::Database> database);

   /*!
    * Write out the class version number and other data members to
    * the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void putSpecializedToDatabase(
        tbox::Pointer<tbox::Database> database); 

   /*!
    * The side iterator iterates over the elements on one axis of a side
    * centered box geometry.  This typedef is a convenience for using the
    * SideIterator<DIM> class.
    */
   typedef SideIterator<DIM> Iterator;

private:
   SideData(const SideData<DIM,TYPE>&);	// not implemented
   void operator=(const SideData<DIM,TYPE>&);		// not implemented

   int d_depth;
   hier::IntVector<DIM> d_directions; 
   ArrayData<DIM,TYPE> d_data[DIM];

};

}
}
#ifndef DEBUG_NO_INLINE
#include "SideData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SideData.C"
#endif
