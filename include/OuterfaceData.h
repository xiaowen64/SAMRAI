//
// File:	OuterfaceData.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 601 $
// Modified:	$Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description:	Templated outerface centered patch data type
//

#ifndef included_pdat_OuterfaceData
#define included_pdat_OuterfaceData

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
#ifndef included_tbox_ArrayData
#include "ArrayData.h"
#endif
#ifndef included_pdat_FaceIndex
#include "FaceIndex.h"
#endif
#ifndef included_pdat_FaceData
#include "FaceData.h"
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

/**
 * Class OuterfaceData<DIM> is a templated face-centered patch data 
 * structure derived from hier::PatchData<DIM>.  It differs from the 
 * FaceData<DIM> class in that, given a box, an outerface data
 * object represents face-centered data living only on the boundary of the
 * box.  It is templated on TYPE with a specified depth (that is, number of
 * components for each spatial location).  See the face geometry and patch
 * data classes for more information about the translation between the 
 * AMR index space and face-centered data.
 *
 * Outerface data is stored in 2*DIM arrays, containing the data for the 
 * patch boundary faces with each of the possible outward pointing normal 
 * directions.  In each of these arrays, memory allocation is in column-major 
 * ordering (e.g., Fortran style) so that the leftmost index runs fastest 
 * in memory.  In addition, the indices are permuted to match the face data
 * storage.  For example, a three-dimensional outerface data object 
 * instantiated with a box [l0:u0,l1:u1,l2:u2] allocates six data (i.e.,
 * three pairs) arrays dimensioned as:
 * \verbatim

     [ l1 : u1 ,
       l2 : u2 , d ]   ,

     [ l2 : u2 ,
       l0 : u0 , d ]   ,

     [ l0 : u0 ,
       l1 : u1 , d ]   ,

 * \endverbatim
 * for the upper and lower x, y, and z (or 0, 1, 2) face directions, 
 * respectively, and where d is the depth index (i.e., number of values at 
 * each face index location).  One- and two-dimensional face data arrays 
 * are managed similary. 
 *
 * IMPORTANT: The OutersideData class provides the same storage
 * as this outerface data class, except that the individual arrays are not
 * permuted in the outerside data type.  Also, outerface and outerside data 
 * classes are intended to interact with their face-centered and side-centered 
 * data counterparts, respectively.  Mixing types, while allowed, is 
 * discouraged to prevent undesirable behavior.
 * 
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::OuterfaceDataFactory
 * @see pdat::OuterfaceIndex
 * @see pdat::OuterfaceIterator
 * @see pdat::OuterfaceGeometry
 */

template<int DIM, class TYPE>
class OuterfaceData : public hier::PatchData<DIM>
{
public:
   /**
    * The constructor for an outerface data object.  The box describes the
    * interior of the index space and the depth gives the number of components
    * for each spatial location in the array.  If the memory arena is not
    * given, then the standard arena is used.  Note that the ghost cell width
    * is currently fixed at zero.
    */
   OuterfaceData(
      const hier::Box<DIM>& box,
      const int depth,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /**
    * The virtual destructor for a outerface data object.
    */
   virtual ~OuterfaceData<DIM,TYPE>();

   /**
    * A fast copy between the source and destination (i.e., this) patch data objects.  
    * Data is copied where there is overlap in the underlying index space.  The copy 
    * is performed on the interior plus the ghost cell width (for both the source and 
    * destination).  Currently, source data must be FaceData.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /**
    * A fast copy between the source and destination.  Member function copy2()
    * is similar to copy() except that the location of source and destination
    * are reversed and copy2() throws an exception (aka dumps core) if it
    * does not understand the type of the argument.  Currently, the destination
    * data must be FaceData.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /**
    * @b Not @b yet @b implemented!
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects.  Currently, this function
    * does nothing since copying to Outerface from Outerface is not defined
    * on all box overlap configurations.
    */
   virtual void copy(
      const hier::PatchData<DIM>& src,
      const hier::BoxOverlap<DIM>& overlap);

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  Member function copy2() is similar to the copy()
    * member function except that the location of source and destination are
    * reversed and copy2() throws an exception (aka dumps core) if it does
    * not understand the type of the argument.  Currently, the destination
    * data must be FaceData.
    */
   virtual void copy2(
      hier::PatchData<DIM>& dst,
      const hier::BoxOverlap<DIM>& overlap) const;

   /**
    * @brief Fast copy between the source and destination at the
    *   specified depths.
    *
    * Data is copied from the source into the destination where there
    * is overlap in the underlying index space.  The copy is performed
    * on the interior plus the ghost cell width (for both the source
    * and destination).
    */
   void copyDepth(int dst_depth,
		  const FaceData<DIM,TYPE>& src,
		  int src_depth);

   /**
    * Determines whether the patch data subclass can estinate the necessary
    * stream size using only index space information.  This routine is defined
    * for the standard built-in types (bool, char, double, float, and int).
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /**
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.  This routine is defined for the
    * standard built-in types (bool, char, double, float, and int).
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /**
    * Pack data lying on the specified index set into the output stream.
    */
   virtual void packStream(
      tbox::AbstractStream& stream,
      const hier::BoxOverlap<DIM>& overlap) const;

   /**
    * Unpack data from the message stream into the specified index set.
    */
   virtual void unpackStream(
      tbox::AbstractStream& stream,
      const hier::BoxOverlap<DIM>& overlap);

   /**
    * Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;

   /**
    * Get a pointer to the beginning of a particular component of the
    * outerface centered array.  The axis gives the X=0, Y=1, or Z=2
    * axis and the face gives the lower=0 or upper=1 face.
    */
   TYPE *getPointer(
      const int axis,
      const int face,
      const int d = 0);

   /**
    * Get a const pointer to the beginning of a particular component of 
    * the outerface centered array.  The axis gives the X=0, Y=1, or Z=2
    * axis and the face gives the lower=0 or upper=1 face.
    */
   const TYPE *getPointer(
      const int axis,
      const int face,
      const int d = 0) const;

   /**
    * Get a pointer to the array data object of a particular component of the
    * outerface centered array.  The axis gives the X=0, Y=1, or Z=2
    * axis and the face gives the lower=0 or upper=1 face.
    */
   ArrayData<DIM,TYPE> &getArrayData(
      const int axis,
      const int face);

   /**
    * Get a const pointer to the array data object of a particular component of 
    * the outerface centered array.  The axis gives the X=0, Y=1, or Z=2
    * axis and the face gives the lower=0 or upper=1 face.
    */
   const ArrayData<DIM,TYPE> &getArrayData(
      const int axis,
      const int face) const;

   /**
    * hier::Index into the outerface data array using a face index.
    */
   TYPE& operator()(const FaceIndex<DIM>& i, const int face, const int d = 0);

   /**
    * hier::Index into the outerface data array (via a const reference) using
    * a face index.
    */
   const TYPE& operator()(
      const FaceIndex<DIM>& i, const int face, const int d = 0) const;

   /**
    * Fill all values of component d with the value t.
    */
   void fill(const TYPE& t, const int d = 0);

   /**
    * Fill all values of component d within the box with the value t.
    */
   void fill(const TYPE& t, const hier::Box<DIM>& box, const int d = 0);

   /**
    * Fill all components with value t.
    */
   void fillAll(const TYPE& t);

   /**
    * Fill all components within the box with value t.
    */
   void fillAll(const TYPE& t, const hier::Box<DIM>& box);

   /**
    * Calculate the amount of memory needed to represent the data in an
    * outerface centered grid.  This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).  If this is not the case,
    * then a specialized function must be defined.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, const int depth);

   /**
    * Print all face centered data residing in the specified box.  If the
    * depth of the array is greater than one, all components are printed.
    * Precision of floating point numbers (i.e., TYPE = float, double, or
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this
    * is ignored.
    */
   void print(const hier::Box<DIM>& box, ostream& os = tbox::plog, int prec = -1) const;

   /**
    * Print the specified component of the face centered data residing in
    * the specified box.  Precision of floating point numbers (i.e.,
    * TYPE = float, double, or dcomplex) can be specified using the provided
    * argument.  The default is 12 decimal places for double and complex
    * floating point numbers, and the default is 6 decimal places floats.
    * For other types, this is ignored.
    */
   void print(const hier::Box<DIM>& box, const int d, ostream& os = tbox::plog,
              int prec = -1) const;

   /**
    * Print all outerface centered data for specified axis index residing
    * in the specified box, axis, and face. If the depth of the data is
    * greater than one, then all components are printed. 
    * Precision of floating point numbers (i.e., TYPE = float, double, or
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this
    * is ignored. 
    */
   void printAxisFace(
      const int axis,
      const int face,
      const hier::Box<DIM>& box,
      ostream& os = tbox::plog,
      int prec = -1) const;

   /**
    * Print specified component for all outerface centered data for the
    * specified axis index, face index, and spatial component residing
    * in the specified box.  Precision of floating point numbers (i.e.,
    * TYPE = float, double, or dcomplex) can be specified using the provided
    * argument.  The default is 12 decimal places for double and complex
    * floating point numbers, and the default is 6 decimal places floats.
    * For other types, this is ignored.
    */
   void printAxisFace(
      const int axis,
      const int face,
      const hier::Box<DIM>& box,
      const int d,
      ostream& os = tbox::plog,
      int prec = -1) const;  

   /**
    * Check that class version and restart file version are equal.  If so,
    * read data members from the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
        tbox::Pointer<tbox::Database> database);

   /**
    * Write out the class version number and other data members to
    * the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void putSpecializedToDatabase(
        tbox::Pointer<tbox::Database> database);

private:
   OuterfaceData(const OuterfaceData<DIM,TYPE>&); // not implemented
   void operator=(const OuterfaceData<DIM,TYPE>&);	  // not implemented

   int d_depth;
   ArrayData<DIM,TYPE> d_data[DIM][2];

   hier::IntVector<DIM> *d_no_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OuterfaceData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuterfaceData.C"
#endif
