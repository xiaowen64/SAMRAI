//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outeredge/OuteredgeData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated outeredge centered patch data type
//

#ifndef included_pdat_OuteredgeData
#define included_pdat_OuteredgeData

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
#ifndef included_pdat_EdgeData
#include "EdgeData.h"
#endif
#ifndef included_pdat_EdgeIndex
#include "EdgeIndex.h"
#endif
#ifndef included_pdat_EdgeOverlap
#include "EdgeOverlap.h"
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
 * @brief Class OuteredgeData<DIM> provides an implementation for data defined
 * at cell edges on the boundaries of AMR patches.  It is derived from the 
 * hier::PatchData interface common to all SAMRAI patch data types.  Given 
 * a CELL-centered AMR index space box, an outeredge data object represents
 * data of some template TYPE and depth on the cell edges on the boundary
 * of the box.  Here, depth indicates the number of data values at each edge 
 * index location.  The OuteredgedgeGeometry class provides the translation 
 * between the standard SAMRAI cell-centered AMR index space and 
 * outeredge-centered data.
 *
 * Outeredge data is stored in 2*DIM*(DIM-2) arrays, each of which contains data
 * associated with edge indices in a coordinate axis direction, an outward 
 * pointing face normal direction, and an upper or lower box face in the 
 * face normal direction.  The data layout in the outernode data arrays matches 
 * the corresponding array sections provided by the node data implementation.
 * Note that outeredge data is NOT defined when the axis and face normal are equal.  
 * This is consistent with the edge data representation.  Where an edge index falls 
 * on more than one box face (patch boundary edges and corners), the outeredge 
 * data value belongs to only one data array so that there are no redundant data 
 * values. Specifically, when DIM > 2, outeredge data boxes are "trimmed" for each 
 * axis edge direction so that each edge index that lives on more than one 
 * face on the box boundary will be associated with the largest face normal 
 * direction and only that face.  Within each array, data is stored in  
 * (i,...,k,d) order, where i,...,k indicates a spatial index and the 
 * d indicates the component depth at that location.  Memory allocation is 
 * in column-major ordering (e.g., Fortran style) so that the leftmost 
 * index runs fastest in memory.
 *
 * To illustrate the outeredge data layout, in particular the "box trimming"
 * that prevents redundant data values, we describe the data for a
 * three-dimensional outeredge data object instantiated over a
 * box [l0:u0,l1:u1,l2:u2] in the standard SAMRAI cell-centered AMR index space.
 * Note: no boxes are trimmed when DIM < 3. 
 * 
 * \verbatim
  
      a = edge axis (corresponds to the standard EdgeData axis)
      f = face normal dir
      s = lower/upper face
      d = data depth (i.e., number of values at each edge index point).
  
      Here axis and face normal values 0, 1, and 2 can be thought of as X, Y, Z 
      respectively.
  
          (a,f,s)
      axis 0:
          (0,0,[0,1])  DATA IS NOT DEFINED when a == f
          (0,1,0)      [l0:u0, l1:l1,     l2+1:u2,   d]
          (0,1,1)      [l0:u0, u1+1:u1+1, l2+1:u2,   d]
          (0,2,0)      [l0:u0, l1:u1+1,   l2:l2,     d]
          (0,2,1)      [l0:u0, l1:u1+1,   u2+1:u2+1, d]
          Note: Edge indices are duplicated at the intersection of faces normal
                to directions 1 and 2.
                So boxes for face normal direction 1 are trimmed in direction 2
                so that edge indices shared with faces in direction 2 only
                appear in face normal direction 2 boxes.
  
      axis 1:
          (1,0,0)      [l0:l0,     l1:u1, l2+1:u2,   d]
          (1,0,1)      [u0+1:u0+1, l1:u1, l2+1:u2,   d]
          (1,1,[0,1])  DATA IS NOT DEFINED when a == f
          (1,2,0)      [l0:u0+1,   l1:u1, l2:l2,     d]
          (1,2,1)      [l0:u0+1,   l1:u1, u2+1:u2+1, d]
          Note: Edge indices are duplicated at the intersection of faces normal
                to directions 0 and 2.
                So boxes for face normal direction 0 are trimmed in direction 2
                so that edge indices shared with faces in direction 2 only
                appear in face normal direction 2 boxes.
  
      axis 2:
          (2,0,0)      [l0:l0,     l1+1:u1,   l2:u2, d]
          (2,0,1)      [u0+1:u0+1, l1+1:u1,   l2:u2, d]
          (2,1,0)      [l0:u0+1,   l1:l1,     l2:u2, d]
          (2,1,1)      [l0:u0+1,   u1+1:u1+1, l2:u2, d]
          (2,2,[0,1])  DATA IS NOT DEFINED when a == f
          Note: Edge indices are duplicated at the intersection of faces normal
                to directions 0 and 1.
                So boxes for face normal direction 0 are trimmed in direction 1
                so that edge indices shared with faces in direction 1 only
                appear in face normal direction 1 boxes.
  
 * \endverbatim
 * Other spatial dimensions are represented similarly. 
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::OuteredgeDataFactory
 * @see pdat::OuteredgeGeometry
 * @see pdat::EdgeIterator
 * @see pdat::EdgeIndex
 */

template <int DIM, class TYPE>
class OuteredgeData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outeredge-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    * 
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outeredge data object will be created.
    *            Note: the ghost cell width is assumed to be zero.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, 
                               int depth);

   /*!
    * @brief Constructor for an outeredge data object.
    * 
    * Note: Outeredge data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outeredge data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   OuteredgeData(const hier::Box<DIM>& box,
                 int depth,
                 tbox::Pointer<tbox::Arena> pool = 
                    tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief Virtual destructor for a outeredge data object.
    */
   virtual ~OuteredgeData<DIM,TYPE>();

   /*!
    * @brief Return the depth (i.e., the number of data values for 
    * each spatial location) of the array.
    */
   int getDepth() const;

   /*!
    * @brief Returns true if outeredge data exists for the given
    * axis and face normal direction; false otherwise.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    */
   bool dataExists(int axis, 
                   int face_normal) const;

   /*!
    * @brief Return the box of valid edge indices for
    *        outeredge data.  Note: the returned box
    *        will reside in the @em edge index space.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM 
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    */
   hier::Box<DIM> getDataBox(int axis,
                             int face_normal, 
                             int side);

   /*!
    * @brief Get a pointer to the beginning of a particular axis,
    * face normal, side, and depth component of the outeredge centered 
    * array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM 
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param depth integer depth component, must satisfy 
    *              0 <= depth < actual depth of data array
    */
   TYPE* getPointer(int axis,
                    int face_normal,
                    int side,
                    int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular axis,
    * face normal, side, and depth component of the outeredge centered
    * array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM 
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param depth integer depth component, must satisfy 
    *              0 <= depth < actual depth of data array
    */
   const TYPE* getPointer(int axis,
                          int face_normal,
                          int side,
                          int depth = 0) const;

   /*!
    * @brief Return a reference to data entry corresponding
    * to a given edge index and depth.
    *
    * @param i const reference to EdgeIndex, @em MUST be 
    *          an index on the outeredge of the box.
    * @param depth integer depth component, must satisfy 
    *              0 <= depth < actual depth of data array
    */
   TYPE& operator()(const EdgeIndex<DIM>& i, 
                    int depth = 0);

   /*!
    * @brief Return a const reference to data entry corresponding
    * to a given edge index and depth.
    *
    * @param i const reference to EdgeIndex, @em MUST be 
    *          an index on the outeredge of the box.
    * @param depth integer depth component, must satisfy 
    *              0 <= depth < actual depth of data array
    */
   const TYPE& operator()(const EdgeIndex<DIM>& i,
                          int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for 
    * given axis, face normal, and side index of the outeredge 
    * centered array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM 
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    */
   ArrayData<DIM,TYPE>& getArrayData(int axis,
                                     int face_normal,
                                     int side);

   /*!
    * @brief Return a const reference to the array data object for 
    * given axis, face normal, and side index of the outeredge 
    * centered array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM 
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    */
   const ArrayData<DIM,TYPE>& getArrayData(int axis,
                                           int face_normal,
                                           int side) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * either EdgeData or OuteredgeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * either EdgeData or OuteredgeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be either EdgeData or OuteredgeData
    * of the same DIM and TYPE and the overlap must be an EdgeOverlap 
    * of the same DIM.  If not, then an unrecoverable error
    * results.
    */
   virtual void copy(const hier::PatchData<DIM>& src,
                     const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be either EdgeData or OuteredgeData
    * of the same DIM and TYPE and the overlap must be an EdgeOverlap 
    * of the same DIM.  If not, then an unrecoverable error
    * results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst,
                      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes edge and outeredge data objects are
    * defined over the same box) from the given edge source data object to
    * this destination outeredge data object at the specified depths.
    */
   void copyDepth(int dst_depth,
                  const EdgeData<DIM,TYPE>& src,
                  int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes edge and outeredge data objects are
    * defined over the same box) to the given edge destination data object
    * from this source outeredge data object at the specified depths.
    */
   void copyDepth2(int dst_depth,
                   EdgeData<DIM,TYPE>& dst,
                   int src_depth) const;

   /*!
    * @brief Add data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be OuteredgeData of the same DIM and
    * TYPE and the overlap must be an EdgeOverlap of the same DIM.  
    * If not, then an unrecoverable error results.
    */
   virtual void sum(const hier::PatchData<DIM>& src,
                    const hier::BoxOverlap<DIM>& overlap);

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
    * EdgeOverlap of the same DIM.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * EdgeOverlap of the same DIM.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Unpack data from stream and add into this patch data object
    * over the specified box overlap region.
    */
   virtual void unpackStreamAndSum(tbox::AbstractStream& stream,
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
    * @brief Print all outeredge data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
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
    * @brief Print all outeredge data values at the given array depth in 
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
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
    * @brief Print all outeredge centered data values for specified axis,
    * face_normal, and side residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.  
    * 
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers 
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored. 
    */
   void printAxisSide(int axis,
                      int face_normal,
                      int side,
                      const hier::Box<DIM>& box,
                      std::ostream& os = tbox::plog,
                      int prec = 12) const;

   /*!
    * @brief Print all outeredge centered data values for specified axis,
    * face_normal, side, and depth residing in the specified box.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param depth integer depth component, must satisfy 
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers 
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored. 
    */
   void printAxisSide(int axis,
                      int face_normal,
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
   OuteredgeData<DIM,TYPE>(const OuteredgeData<DIM,TYPE>&); // not implemented
   void operator=(const OuteredgeData<DIM,TYPE>&);	  // not implemented

   //@{
   //! @name Internal implementations for data copy operations.
   void copyFromEdge( const EdgeData<DIM,TYPE> &src );
   void copyFromEdge( const EdgeData<DIM,TYPE> &src,
		      const EdgeOverlap<DIM> &overlap  );
   void copyToEdge( EdgeData<DIM,TYPE> &dst ) const;
   void copyToEdge( EdgeData<DIM,TYPE> &dst,
		    const EdgeOverlap<DIM> &overlap  ) const;
   void copyFromOuteredge( const OuteredgeData<DIM,TYPE> &src );
   void copyFromOuteredge( const OuteredgeData<DIM,TYPE> &src,
			   const EdgeOverlap<DIM> &overlap );
   void copyToOuteredge( OuteredgeData<DIM,TYPE> &dst ) const;
   void copyToOuteredge( OuteredgeData<DIM,TYPE> &dst,
			 const EdgeOverlap<DIM> &overlap ) const;
   //@}
		 

   int d_depth;
   ArrayData<DIM,TYPE> d_data[DIM][DIM][2];

};

}
}

#ifndef DEBUG_NO_INLINE
#include "OuteredgeData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuteredgeData.C"
#endif
