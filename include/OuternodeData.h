//
// File:	OuternodeData.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated outernode centered patch data type
//

#ifndef included_pdat_OuternodeData
#define included_pdat_OuternodeData

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
#ifndef included_pdat_NodeData
#include "NodeData.h"
#endif
#ifndef included_pdat_NodeIndex
#include "NodeIndex.h"
#endif
#ifndef included_pdat_NodeOverlap
#include "NodeOverlap.h"
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
 * @brief Class OuternodeData is a templated node-centered patch data
 * structure derived from hier::PatchData.  It differs from the
 * NodeData<DIM> class in that, given a box, an outernode data
 * object represents node-centered data living only on the border of the
 * patch.
 *
 * It is templated on TYPE with a specified depth (that is, number of
 * components for each spatial location).  See the node geometry and patch
 * data classes for more information about the translation between the
 * AMR index space and node-centered data.
 *
 * Outernode data is stored in 2*DIM arrays, containing the data for the 
 * patch boundary sides with each of the possible outward pointing normal 
 * directions. Where an outernode falls on more than one side (patch edges 
 * and corners), the outernode belongs to the array associated with the 
 * higher dimensional direction. In each of these arrays, memory allocation 
 * is in column-major ordering (e.g., Fortran style) so that the leftmost 
 * index runs fastest in memory.  For example, a three-dimensional outernode 
 * data object instantiated with a box [l0:u0,l1:u1,l2:u2] allocates six 
 * data (i.e., three pairs) arrays dimensioned as:
 * \verbatim
 *
 *    
 *    X:  [ l1+1 : u1   ,
 *          l2+1 : u2   , d ]   ,
 *
 *    Y:  [ l0   : u0+1 ,
 *          l2+1 : u2   , d ]   ,
 *
 *    Z:  [ l0   : u0+1 ,
 *          l1   : u1+1 , d ]   ,
 *
 * \endverbatim
 * for the upper and lower x, y, and z (or 0, 1, 2) face directions, 
 * respectively.  One- and two-dimensional node data arrays are managed 
 * similary.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::OuternodeDataFactory
 * @see pdat::OuternodeIndex
 * @see pdat::NodeIterator
 * @see pdat::OuternodeGeometry
 */

template <int DIM, class TYPE>
class OuternodeData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Constructor for an outernode data object.
    *
    * @param box describes the interior of the index space
		Note that the ghost cell width is currently fixed at zero.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   OuternodeData(
      const hier::Box<DIM>& box,
      const int depth,
      tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief Virtual destructor for a outernode data object.
    */
   virtual ~OuternodeData<DIM,TYPE>();

   /*!
    * @brief A fast copy between the source and destination (i.e., this)
    * patch data objects.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * NodeData<DIM>.  If copy() does not understand the source object type,
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
    * Currently, the destination data must be NodeData<DIM>.
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
    * Outernode from Outernode is not defined on all box overlap configurations.
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
    * Currently, the destination data must be NodeData<DIM>.
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
                  const NodeData<DIM,TYPE>& src,
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
    * @brief Returns whether outernodes exists for the side normal
    * to a given axis.
    *
    * The axis gives the I=0, J=1, or K=2 axis.
    *
    * Recall that outernodes that lie on corners or edges
    * are seen on multiple sides of a box.  By convention,
    * duplicated nodes belong to the side of the box normal
    * to the axis with the higher dimension.  If a side
    * consists of only outnerodes that are owned by a higher
    * dimension, no data exists on that side.  This method
    * is used to determine whether any data exists on the
    * given side.
    */
   bool dataExists( const int axis ) const;

   /*!
    * @brief
    * Get a pointer to the beginning of a particular component of the
    * outernode centered array.
    *
    * The axis gives the X=0, Y=1, or Z=2
    * axis and the side gives the lower=0 or upper=1 side.
    * See class description for the size of the array returned.
    */
   TYPE *getPointer(
      const int axis,
      const int side,
      const int d = 0);

   /*!
    * @brief
    * Get a const pointer to the beginning of a particular component of 
    * the outernode centered array.
    *
    * The axis gives the X=0, Y=1, or Z=2
    * axis and the side gives the lower=0 or upper=1 side.
    * See class description for the size of the array returned.
    */
   const TYPE *getPointer(
      const int axis,
      const int side,
      const int d = 0) const;

   /*!
    * @brief
    * Get a pointer to the array data object of a particular component of the
    * outernode centered array.
    *
    * The axis gives the X=0, Y=1, or Z=2
    * axis and the side gives the lower=0 or upper=1 side.
    */
   ArrayData<DIM,TYPE> &getArrayData(
      const int axis,
      const int side);

   /*!
    * @brief
    * Get a const pointer to the array data object of a particular component of 
    * the outernode centered array.
    *
    * The axis gives the X=0, Y=1, or Z=2
    * axis and the side gives the lower=0 or upper=1 side.
    */
   const ArrayData<DIM,TYPE> &getArrayData(
      const int axis,
      const int side) const;

   /*!
    * @brief
    * Index<DIM> into the outernode data array using a node index.
    *
    * The index @em MUST be an index on the outernode of the box.
    */
   TYPE& operator()(const NodeIndex<DIM>& i, const int depth = 0);

   /*!
    * @brief
    * Index<DIM> into the outernode data array (via a const reference) using
    * a node index.
    *
    * The index @em MUST be an index on the outernode of the box.
    */
   const TYPE& operator()(
      const NodeIndex<DIM>& i, const int depth = 0) const;

   /*
    * @brief Return the box of valid node indices.
    *
    * @param dim Dimension, should be in [0:DIM-1]
    * @param side Side, should be 0 for the lower indexed side
    *        or 1 for the higher indexed side.
    */
   hier::Box<DIM> getDataBox( int dim, int side );

   /*!
    * @brief
    * Fill all values of component d with the value t.
    */
   void fill(const TYPE& t, const int d = 0);

   /*!
    * @brief
    * Fill all values of component d within the box with the value t.
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
    * outernode centered grid.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, const int depth);

   /*!
    * @brief
    * Print all outernode data residing in the specified box.
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
    * Print the specified component of the outernode data residing in
    * the specified box.
    *
    * Precision of floating point numbers (i.e.,
    * TYPE = float, double, or dcomplex) can be specified using the provided
    * argument.  The default is 12 decimal places for double and complex
    * floating point numbers, and the default is 6 decimal places floats.
    * For other types, this is ignored.
    */
   void print(const hier::Box<DIM>& box, const int d, ostream& os = tbox::plog,
              int prec = -1) const;

   /*!
    * @brief
    * Print all outernode centered data for specified axis index residing
    * in the specified box, axis, and side.
    *
    * If the depth of the data is
    * greater than one, then all components are printed.  
    * Precision of floating point numbers (i.e., TYPE = float, double, or
    * dcomplex) can be specified using the provided argument.  The default
    * is 12 decimal places for double and complex floating point numbers,
    * and the default is 6 decimal places floats.  For other types, this
    * is ignored.
    */
   void printAxisSide(
      const int axis,
      const int side,
      const hier::Box<DIM>& box,
      ostream& os = tbox::plog,
      int prec = -1) const;

   /*!
    * @brief
    * Print specified component for all outernode centered data for the
    * specified axis index, side index, and spatial component residing
    * in the specified box.
    *
    * Precision of floating point numbers (i.e.,
    * TYPE = float, double, or dcomplex) can be specified using the provided
    * argument.  The default is 12 decimal places for double and complex
    * floating point numbers, and the default is 6 decimal places floats.
    * For other types, this is ignored.
    */
   void printAxisSide(
      const int axis,
      const int side,
      const hier::Box<DIM>& box,
      const int d,
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
   OuternodeData<DIM,TYPE>(const OuternodeData<DIM,TYPE>&); // not implemented
   void operator=(const OuternodeData<DIM,TYPE>&);	  // not implemented

   //@{
   //! @name Internal implementations for data copy interfaces.
   void copyFromNode( const NodeData<DIM,TYPE> &src );
   void copyFromNode( const NodeData<DIM,TYPE> &src,
		      const NodeOverlap<DIM> &overlap  );
   void copyToNode( NodeData<DIM,TYPE> &dst ) const;
   void copyToNode( NodeData<DIM,TYPE> &dst,
		    const NodeOverlap<DIM> &overlap  ) const;
   void copyFromOuternode( const OuternodeData<DIM,TYPE> &src );
   void copyFromOuternode( const OuternodeData<DIM,TYPE> &src,
			   const NodeOverlap<DIM> &overlap );
   void copyToOuternode( OuternodeData<DIM,TYPE> &dst ) const;
   void copyToOuternode( OuternodeData<DIM,TYPE> &dst,
			 const NodeOverlap<DIM> &overlap ) const;
   //@}
		 

   int d_depth;
   ArrayData<DIM,TYPE> d_data[DIM][2];

   hier::IntVector<DIM> *d_no_ghosts;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "OuternodeData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuternodeData.C"
#endif
