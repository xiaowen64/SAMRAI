//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/node/NodeIndex.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Index for node centered patch data types
//

#ifndef included_pdat_NodeIndex
#define included_pdat_NodeIndex

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Index
#include "Index.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class NodeIndex<DIM> implements a simple n-dimensional integer
 * vector for node centered variables.  Given a hier::Box in the AMR abstract
 * index space, the index space for a node-centered variable runs from the
 * lower corner of the box to the upper corner of the box plus one in each
 * dimension.  See the node box geometry class for more information about
 * the mapping between the AMR index space and the node indices.
 *
 * @see hier::Index
 * @see pdat::NodeData
 * @see pdat::NodeGeometry
 * @see pdat::NodeIterator
 */

template<int DIM> class NodeIndex : public hier::Index<DIM>
{
public:
   /**
    * The Corner enumerated type is used when converting from a cell centered
    * index to a node centered index.  In 1d, use Left and Right.  In 2d, use
    * LowerLeft, LowerRight, UpperLeft, and UpperRight.  In 3d, the naming is
    * less intuitive, and use names LLL through UUU, where L means lower and
    * U means upper.  Therefore, to get the box upper in X, lower in Y, and
    * lower in Z, use corner name ULL.
    */
   enum Corner {
      Left = 0, Right = 1,
      LowerLeft = 0, LowerRight = 1, UpperLeft = 2, UpperRight = 3,
      LLL = 0, ULL = 1, LUL = 2, UUL = 3, LLU = 4, ULU = 5, LUU = 6, UUU = 7
   };

   /**
    * The default constructor for a node index creates an uninitialized index.
    */
   NodeIndex();

   /**
    * Construct a node index from a regular index and a corner.
    *
    * The Corner enumerated type is only defined for 3D or lower, so use
    * the next constructor with an hier::IntVector argument when using higher
    * dimensions.
    */
   NodeIndex(const hier::Index<DIM>& rhs, const Corner corner);

   /**
    * Construct a node index from a regular index and an hier::IntVector.  The
    * hier::IntVector is binary--an assertion failure will result if it contains
    * any values other than 0 or 1.  For each dimension, if the hier::IntVector
    * contains a 0, the node index will represent a lower bound in that
    * dimensional direction, and if 1 will represent an upper bound in that
    * direction.
    */
   NodeIndex(const hier::Index<DIM>& rhs, const hier::IntVector<DIM>& corner);

   /**
    * The copy constructor creates a node index equal to the argument.
    */
   NodeIndex(const NodeIndex<DIM>& rhs);

   /**
    * The assignment operator sets the node index equal to the argument.
    */
   NodeIndex<DIM>& operator=(const NodeIndex<DIM>& rhs);

   /**
    * The node index destructor does nothing interesting.
    */
   ~NodeIndex<DIM>();

   /**
    * Plus-equals operator for a node index and an integer vector.
    */
   NodeIndex<DIM>& operator+=(const hier::IntVector<DIM>& rhs);

   /**
    * Plus operator for a node index and an integer vector.
    */
   NodeIndex<DIM> operator+(const hier::IntVector<DIM>& rhs) const;

   /**
    * Plus-equals operator for a node index and an integer.
    */
   NodeIndex<DIM>& operator+=(const int rhs);

   /**
    * Plus operator for a node index and an integer.
    */
   NodeIndex<DIM> operator+(const int rhs) const;

   /**
    * Minus-equals operator for a node index and an integer vector.
    */
   NodeIndex<DIM>& operator-=(const hier::IntVector<DIM>& rhs);

   /**
    * Minus operator for a node index and an integer vector.
    */
   NodeIndex<DIM> operator-(const hier::IntVector<DIM>& rhs) const;

   /**
    * Minus-equals operator for a node index and an integer.
    */
   NodeIndex<DIM>& operator-=(const int rhs);

   /**
    * Minus operator for a node index and an integer.
    */
   NodeIndex<DIM> operator-(const int rhs) const;

   /**
    * Times-equals operator for a node index and an integer vector.
    */
   NodeIndex<DIM>& operator*=(const hier::IntVector<DIM>& rhs);

   /**
    * Times operator for a node index and an integer vector.
    */
   NodeIndex<DIM> operator*(const hier::IntVector<DIM>& rhs) const;

   /**
    * Times-equals operator for a node index and an integer.
    */
   NodeIndex<DIM>& operator*=(const int rhs);

   /**
    * Times operator for a node index and an integer.
    */
   NodeIndex<DIM> operator*(const int rhs) const;

   /**
    * Returns true if two node index objects are equal.
    * All components must be the same for equality.
    */
   bool operator==(const NodeIndex<DIM>& rhs) const;

   /**
    * Returns true if two node index objects are not equal.
    * Any of the components may be different for inequality.
    */
   bool operator!=(const NodeIndex<DIM>& rhs) const;

private:
   /*
    * Initializes the offsets if it has not yet been done
    */
   void setOffsets();

   static hier::IntVector<DIM> s_offsets[2 << DIM];
   static bool s_offsets_are_set;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "NodeIndex.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NodeIndex.C"
#endif
