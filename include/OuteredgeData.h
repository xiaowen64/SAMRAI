//
// File:	OuteredgeData.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 692 $
// Modified:	$Date: 2005-10-28 13:37:46 -0700 (Fri, 28 Oct 2005) $
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
using namespace std;
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
 * @brief Class OuteredgeData<DIM> is a templated edge-centered patch data
 * structure derived from hier::PatchData<DIM>.  It differs from the
 * EdgeData<DIM> class in that, given a box, an outeredge data
 * object represents edge-centered data living only on the border of the
 * patch.
 *
 * It is templated on TYPE with a specified depth (that is, number of
 * components for each spatial location).  See the edge geometry and patch
 * data classes for more information about the translation between the
 * AMR index space and edge-centered data.
 *
 * Outeredge data is stored in DIM*DIM*2 arrays, containing the data for the 
 * patch boundary sides with each of the possible outward pointing normal 
 * directions. Where an outeredge falls on more than one side (patch edges 
 * and corners), the outeredge belongs to the array associated with the 
 * higher dimensional direction. In each of these arrays, memory allocation 
 * is in column-major ordering (e.g., Fortran style) so that the leftmost 
 * index runs fastest in memory.  
 *
 * The outeredge data is related to the edge data in the following way:
 *
 *    Outeredge box(axis,face_nrml,s) =  EdgeData<DIM>.getBox(axis) 
 *
 * where "axis" corresponds to the box of the standard edge datatype,
 * "face_nrml" is the normal face direction, and "s" indicates the upper
 * or lower face.  Note that when edge_dir = face_dir, there are no outside
 * edges so the data is NULL.
 *
 * A three-dimensional outeredge data object instantiated with 
 * a box [l0:u0,l1:u1,l2:u2] allocates 12 data (i.e., 3x2 pairs) arrays 
 * dimensioned as:
 *
 * \verbatim
 *
 *    a = edge axis
 *    f = face normal dir
 *    s = lower/upper face
 *
 *        (a,f,s) 
 *    0:  (0,0,[0,1])  NULL
 *        (0,1,[0,1])  [l0:u0,l1:l1,l2+1:u2,d],  [l0:u0,u1:u1,l2+1:u2,d]
 *        (0,2,[0,1])  [l0:u0,l1:u1+1,l2:l2,d],  [l0:u0,l1:u1+1,u2:u2,d]
 *        note: trimmed in 2, not trimmed in 1
 *
 *    1:  (1,0,[0,1])  [l0:l0,l1:u1,l2+1:u2,d],  [u0:u0,l1:u1,l2+1:u2,d]
 *        (1,1,[0,1])  NULL
 *        (1,2,[0,1])  [l0:u0+1,l1:u1,l2:l2,d],  [l0:u0+1,l1:u1,u2:u2,d]
 *        note: trimmed in 2, not trimmed in 0
 *
 *    2:  (2,0,[0,1])  [l0:l0,l1+1:u1,l2:u2,d],  [u0:u0,l1+1:u1,l2:u2,d]
 *        (2,1,[0,1])  [l0:u0+1,l1:l1,l2:u2,d],  [l0:u0+1,u1:u1,l2:u2,d]
 *        (2,2,[0,1])  NULL
 *        note: trimmed in 1, not trimmed in 0
 *
 * \endverbatim
 * where 0, 1, and 2 can be thought of as X, Y, Z respectively, and d is the
 * depth of the data.  One- and two-dimensional edge data arrays are managed 
 * similary.  The "a" dimension corresponds with the "axis" of standard 
 * EdgeData<DIM>.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see ArrayData<DIM>
 * @see hier::PatchData<DIM>
 * @see OuteredgeDataFactory<DIM>
 * @see OuteredgeIndex<DIM>
 * @see EdgeIterator<DIM>
 * @see OuteredgeGeometry<DIM>
 */

template <int DIM, class TYPE>
class OuteredgeData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Constructor for an outeredge data object.
    *
    * @param box describes the interior of the index space
		Note that the ghost cell width is currently fixed at zero.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   OuteredgeData(
      const hier::Box<DIM>& box,
      const int depth,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief Virtual destructor for a outeredge data object.
    */
   virtual ~OuteredgeData<DIM,TYPE>();

   /*!
    * @brief A fast copy between the source and destination (i.e., this)
    * patch data objects.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * EdgeData<DIM>.  If copy() does not understand the source object type,
    * then copy2() is called.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy between the source and destination.  
    *
    * Member function copy2() is similar to copy() except that the location 
    * of source and destination are reversed and copy2() throws an 
    * exception (aka dumps core) if it does not understand the type of 
    * the argument.
    *
    * Currently, the destination data must be EdgeData<DIM>.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief
    * Copy data from the source into the destination using the designated
    * overlap descriptor.
    *
    * The overlap description will have been computed
    * using the appropriate box geometry objects.
    * If copy() does not understand the source object type, then copy2()
    * is called.  Currently, this function does nothing since copying to
    * Outeredge from Outeredge is not defined on all box overlap configurations.
    */
   virtual void copy(
      const hier::PatchData<DIM>& src,
      const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief
    * Copy data from the source into the destination using the designated
    * overlap descriptor.
    *
    * Member function copy2() is similar to the copy()
    * member function except that the location of source and destination are
    * reversed and copy2() throws an exception (aka dumps core) if it does
    * not understand the type of the argument.
    * Currently, the destination data must be EdgeData<DIM>.
    */
   virtual void copy2(
      hier::PatchData<DIM>& dst,
      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Fast copy between the source and destination at the
    * specified depths.
    *
    * Data is copied from the source into the destination where there
    * is overlap in the underlying index space.  The copy is performed
    * on the interior plus the ghost cell width (for both the source
    * and destination).
    */
   void copyDepth(int dst_depth,
                  const EdgeData<DIM,TYPE>& src,
                  int src_depth);

   /*!
    * @brief
    * Determines whether the patch data subclass can estinate the necessary
    * stream size using only index space information.
    *
    * This routine is defined
    * for the standard built-in types (bool, char, double, float, and int).
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /*!
    * @brief
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.
    *
    * This routine is defined for the
    * standard built-in types (bool, char, double, float, and int).
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief
    * Pack data lying on the specified index set into the output stream.
    */
   virtual void packStream(
      tbox::AbstractStream& stream,
      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief
    * Unpack data from the message stream into the specified index set.
    */
   virtual void unpackStream(
      tbox::AbstractStream& stream,
      const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief
    * Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;

   /*!
    * @brief Returns whether outeredges exists for the face normal
    * to a given axis.
    *
    * The "axis" and "face_nrml" arguments correspond to the edge
    * axis and face normal, as discussed above.  
    */
   bool dataExists( int axis, 
                    int face_nrml ) const;

   /*!
    * @brief
    * Get a pointer to the beginning of a particular component of the
    * outeredge centered array.
    *
    * The "axis" and "face_nrml" arguments specify the X=0, Y=1, or Z=2
    * axis and "s" specifies the lower=0 or upper=1 side.
    * See class description for the size of the array returned.
    */
   TYPE *getPointer(
      const int axis,
      const int face_nrml,
      const int s,
      const int d = 0);

   /*!
    * @brief
    * Get a const pointer to the beginning of a particular component of 
    * the outeredge centered array.
    *
    * The "axis" and "face_nrml" arguments specify the X=0, Y=1, or Z=2
    * axis and "s" specifies the lower=0 or upper=1 side.
    * See class description for the size of the array returned.
    */
   const TYPE *getPointer(int axis,
                          int face_nrml,
                          int s,
                          int d = 0) const;

   /*!
    * @brief
    * Get a pointer to the array data object of a particular component of the
    * outeredge centered array.
    *
    * The "axis" and "face_nrml" arguments specify the X=0, Y=1, or Z=2
    * axis and "s" specifies the lower=0 or upper=1 side.
    */
   ArrayData<DIM,TYPE> &getArrayData(int axis,
                                     int face_nrml,
                                     int s);

   /*!
    * @brief
    * Get a const pointer to the array data object of a particular component of the
    * outeredge centered array.
    *
    * The "axis" and "face_nrml" arguments specify the X=0, Y=1, or Z=2
    * axis and "s" specifies the lower=0 or upper=1 side.
    */
   const ArrayData<DIM,TYPE> &getArrayData(int axis,
                                           int face_nrml,
                                           int s) const;

   /*!
    * @brief
    * Index<DIM> into the outeredge data array using a edge index.
    *
    * The index @em MUST be an index on the outeredge of the box.
    */
   TYPE& operator()(const EdgeIndex<DIM>& i, 
                    int axis,
                    int depth = 0);

   /*!
    * @brief
    * Index<DIM> into the outeredge data array (via a const reference) using
    * a edge index.
    *
    * The index @em MUST be an index on the outeredge of the box.
    */
   const TYPE& operator()(const EdgeIndex<DIM>& i,
                          int axis,
                          int depth = 0) const;

   /*
    * @brief Return the box of valid edge indices.
    *
    * @param axis should be in [0:DIM-1]
    * @param face_nrml should be in [0:DIM-1]
    * @param s side should be 0 for the lower indexed side
    *        or 1 for the higher indexed side.
    */
   hier::Box<DIM> getDataBox( int axis,
                              int face_nrml, 
                              int s );

   /*!
    * @brief
    * Fill all values of component depth d with the value t.
    */
   void fill(const TYPE& t, const int d = 0);

   /*!
    * @brief
    * Fill all values of component depth d within the box with the value t.
    */
   void fill(const TYPE& t, const hier::Box<DIM>& box, const int d = 0);

   /*!
    * @brief
    * Fill all components with value t.
    */
   void fillAll(const TYPE& t);

   /*!
    * @brief
    * Fill all components within the box with value t.
    */
   void fillAll(const TYPE& t, const hier::Box<DIM>& box);

   /*!
    * @brief
    * Calculate the amount of memory needed to represent the data in an
    * outeredge centered grid.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, const int depth);

   /*!
    * @brief
    * Print all outeredge data residing in the specified box.
    *
    * If the
    * depth of the array is greater than one, all components are printed.
    * Precision of floating point numbers (i.e., TYPE = float, double, or
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this
    * is ignored.
    */
   void print(const hier::Box<DIM>& box, ostream& os = tbox::plog, int prec = -1) const;

   /*!
    * @brief
    * Print the specified component of the outeredge data residing in
    * the specified box.
    *
    * Precision of floating point numbers (i.e.,
    * TYPE = float, double, or dcomplex) can be specified using the provided
    * argument.  The default is 12 decimal places for double and complex
    * floating point numbers, and the default is 6 decimal places floats.
    * For other types, this is ignored.
    */
   void print(const hier::Box<DIM>& box, int d, ostream& os = tbox::plog,
              int prec = -1) const;

   /*!
    * @brief
    * Print all outeredge centered data for specified (axis,face_nrml) index 
    * residing in the specified box and side.
    *
    * If the depth of the data is
    * greater than one, then all components are printed.  
    * Precision of floating point numbers (i.e., TYPE = float, double, or
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this
    * is ignored.
    */
   void printAxisSide(int axis,
                      int face_nrml,
                      int s,
                      const hier::Box<DIM>& box,
                      ostream& os = tbox::plog,
                      int prec = -1) const;

   /*!
    * @brief
    * Print all outeredge centered data for specified (axis,face_nrml) index 
    * residing in the specified box and side.
    *
    * If the depth of the data is
    * greater than one, then all components are printed.  
    * Precision of floating point numbers (i.e., TYPE = float, double, or
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this
    * is ignored.
    */
   void printAxisSide(int axis,
                      int face_nrml,
                      int s,
                      const hier::Box<DIM>& box,
                      int d,
                      ostream& os = tbox::plog,
                      int prec = -1) const;

   /*!
    * @brief
    * Check that class version and restart file version are equal.  If so,
    * read data members from the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
        tbox::Pointer<tbox::Database> database);

   /*!
    * @brief
    * Write out the class version number and other data members to
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
   //! @name Internal implementations for data copy interfaces.
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

   hier::IntVector<DIM> *d_no_ghosts;

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
