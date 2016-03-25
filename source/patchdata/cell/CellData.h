//
// File:	CellData.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated cell centered patch data type
//

#ifndef included_pdat_CellData
#define included_pdat_CellData

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
#ifndef included_pdat_CellIndex
#include "CellIndex.h"
#endif
#ifndef included_pdat_CellIterator
#include "CellIterator.h"
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
 * @brief Class CellData<DIM> manages data defined on cell centers in 
 * patches.  It is a templated cell-centered patch data structure
 * derived from hier::PatchData<DIM>.  Given a box, a cell data object represents
 * cell-centered data of some template TYPE with a specified depth (that is,
 * number of components for each spatial location).  See the cell geometry
 * class for more information about the translation between the AMR index
 * space and cell-centered data.
 *
 * A cell data array is stored in (i,...,k,d) order, where i,...,k indicates
 * spatial indices and d indicates the component depth at that location.
 * Memory allocation is in colum-major ordering (e.g., Fortran style)
 * so that the leftmost index runs fastest in memory.  For example, a 
 * three- dimensional cell data object instantiated with a box 
 * [l0:u0,l1:u1,l2:u2] allocates a data array dimensioned as 
 * \verbatim

     [ l0 : u0 ,
       l1 : u1 ,
       l2 : u2 , d ]

 * \endverbatim
 * One- and two-dimensional cell data arrays are managed similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::CellDataFactory
 * @see pdat::CellIndex
 * @see pdat::CellIterator
 * @see pdat::CellGeometry
 */

template<int DIM, class TYPE>
class CellData : public hier::PatchData<DIM>
{
public:
   /*!
    * The constructor for a cell data object.  The box describes the interior
    * of the index space and the ghosts vector describes the ghost cells in
    * each coordinate direction.  The depth gives the number of components
    * for each spatial location in the array.  If the memory arena is not
    * given, then the standard arena is used.
    */
   CellData(
      const hier::Box<DIM>& box,
      const int depth,
      const hier::IntVector<DIM>& ghosts,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * The virtual destructor for a cell data object.
    */
   virtual ~CellData<DIM,TYPE>();

   /*!
    * A fast copy between the source and destination.  Data is copied from
    * the source into the destination where there is overlap in the underlying
    * index space.  The copy is performed on the interior plus the ghost cell
    * width (for both the source and destination).  If copy() does not
    * understand the source object type, then copy2() is called.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * A fast copy between the source and destination.  Member function copy2()
    * is similar to copy() except that the location of source and destination
    * are reversed and copy2() throws an exception (aka dumps core) if it
    * does not understand the type of the argument.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects.  If copy() does not
    * understand the source object type, then copy2() is called.
    */
   virtual void copy(
      const hier::PatchData<DIM>& src, const hier::BoxOverlap<DIM>& overlap);

   /*!
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  Member function copy2() is similar to the copy()
    * member function except that the location of source and destination are
    * reversed and copy2() throws an exception (aka dumps core) if it does
    * not understand the type of the argument.
    */
   virtual void copy2(
      hier::PatchData<DIM>& dst, const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Fast copy between the source and destination at the
    *   specified depths.
    *
    * Data is copied from the source into the destination where there
    * is overlap in the underlying index space.  The copy is performed
    * on the interior plus the ghost cell width (for both the source
    * and destination).
    */
   virtual void copyDepth(int dst_depth,
			  const CellData<DIM,TYPE>& src,
			  int src_depth);

   /*!
    * Determines whether the patch data subclass can estimate the necessary
    * stream size using only index space information.  This routine is defined
    * for the standard built-in types (bool, char, double, float, and int).
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /*!
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.  This routine is defined for the
    * standard built-in types (bool, char, double, float, and int).
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * Pack data lying on the specified index set into the output stream.
    */
   virtual void packStream(
      tbox::AbstractStream& stream, const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * Unpack data from the message stream into the specified index set.
    */
   virtual void unpackStream(
      tbox::AbstractStream& stream, const hier::BoxOverlap<DIM>& overlap);

   /*!
    * Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;

   /*!
    * Return reference to the array data object.
    */
   ArrayData<DIM,TYPE>& getArrayData();

   /*!
    * Return a const reference to the array data object.
    */
   const ArrayData<DIM,TYPE>& getArrayData() const;

   /*!
    * Get a pointer to the beginning of a particular component of the
    * cell centered array.
    */
   TYPE *getPointer(const int d = 0);

   /*!
    * Get a const pointer to the beginning of a particular component
    * of the cell centered array.
    */
   const TYPE *getPointer(const int d = 0) const;

   /*!
    * hier::Index into the cell data array using a cell index.
    */
   TYPE& operator()(const CellIndex<DIM>& i, const int d = 0);

   /*!
    * hier::Index into the cell data array (via a const reference) using
    * a cell index.
    */
   const TYPE& operator()(const CellIndex<DIM>& i, const int d = 0) const;

   /*!
    * Fill all values of component d with the value t.
    */
   void fill(const TYPE& t, const int d = 0);

   /*!
    * Fill all values of component d within the box with the value t.
    */
   void fill(const TYPE& t, const hier::Box<DIM>& box, const int d = 0);

   /*!
    * Fill all components with value t.
    */
   void fillAll(const TYPE& t);

   /*!
    * Fill all components within the box with value t.
    */
   void fillAll(const TYPE& t, const hier::Box<DIM>& box);

   /*!
    * Copy data from supplied source over the supplied box.
    */
   void copyOnBox(const CellData<DIM,TYPE>& src, const hier::Box<DIM>& box);
   
   /*!
    * Calculate the amount of memory needed to represent the data in a
    * cell centered grid.  This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    */
   static size_t getSizeOfData(
      const hier::Box<DIM>& box, const int depth, const hier::IntVector<DIM>& ghosts);

   /*!
    * Print all cell centered data residing in the specified box.  If the
    * depth of the array is greater than one, all components are printed.
    * Precision of floating point numbers (i.e., TYPE = float, double, or 
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this 
    * is ignored.
    */
   void print(const hier::Box<DIM>& box, ostream& os = tbox::plog, int prec = -1) const;

   /*!
    * Print the specified component of the cell centered data residing in
    * the specified box.  Precision of floating point numbers (i.e., 
    * TYPE = float, double, or dcomplex) can be specified using the provided 
    * argument.  The default is 12 decimal places for double and complex 
    * floating point numbers, and the default is 6 decimal places floats.  
    * For other types, this is ignored.
    */
   void print(const hier::Box<DIM>& box, const int d, ostream& os = tbox::plog, 
              int prec = -1) const;

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
    * The cell iterator iterates over the elements of a cell
    * centered box geometry.  This typedef is a convenient link
    * to the CellIterator<DIM> class.
    */
   typedef CellIterator<DIM> Iterator;

private:
   CellData(const CellData<DIM,TYPE>&);	// not implemented
   void operator=(const CellData<DIM,TYPE>&);		// not implemented

   int d_depth;
   ArrayData<DIM,TYPE> d_data;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "CellData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellData.C"
#endif
