//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/cell/CellIndex.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Index for cell centered patch data types
//

#ifndef included_pdat_CellIndex
#define included_pdat_CellIndex

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
 * Class CellIndex<DIM> implements a simple n-dimensional integer
 * vector for cell centered variables.  Cell indices contain an integer
 * index location in AMR index space and are identical to the AMR indices.
 *
 * @see hier::Index
 * @see pdat::CellData
 * @see pdat::CellGeometry
 * @see pdat::CellIterator
 */

template<int DIM> class CellIndex : public hier::Index<DIM>
{
public:
   /**
    * The default constructor for a cell index creates an uninitialized index.
    */
   CellIndex();

   /**
    * Construct a cell index from a regular AMR index.
    */
   CellIndex(const hier::Index<DIM>& rhs);

   /**
    * The copy constructor creates a cell index equal to the argument.
    */
   CellIndex(const CellIndex<DIM>& rhs);

   /**
    * The assignment operator sets the cell index equal to the argument.
    */
   CellIndex<DIM>& operator=(const CellIndex<DIM>& rhs);

   /**
    * The cell index destructor does nothing interesting.
    */
   ~CellIndex<DIM>();

   /**
    * Plus-equals operator for a cell index and an integer vector.
    */
   CellIndex<DIM>& operator+=(const hier::IntVector<DIM>& rhs);

   /**
    * Plus operator for a cell index and an integer vector.
    */
   CellIndex<DIM> operator+(const hier::IntVector<DIM>& rhs) const;

   /**
    * Plus-equals operator for a cell index and an integer.
    */
   CellIndex<DIM>& operator+=(const int rhs);

   /**
    * Plus operator for a cell index and an integer.
    */
   CellIndex<DIM> operator+(const int rhs) const;

   /**
    * Minus-equals operator for a cell index and an integer vector.
    */
   CellIndex<DIM>& operator-=(const hier::IntVector<DIM>& rhs);

   /**
    * Minus operator for a cell index and an integer vector.
    */
   CellIndex<DIM> operator-(const hier::IntVector<DIM>& rhs) const;

   /**
    * Minus-equals operator for a cell index and an integer.
    */
   CellIndex<DIM>& operator-=(const int rhs);

   /**
    * Minus operator for a cell index and an integer.
    */
   CellIndex<DIM> operator-(const int rhs) const;

   /**
    * Times-equals operator for a cell index and an integer vector.
    */
   CellIndex<DIM>& operator*=(const hier::IntVector<DIM>& rhs);

   /**
    * Times operator for a cell index and an integer vector.
    */
   CellIndex<DIM> operator*(const hier::IntVector<DIM>& rhs) const;

   /**
    * Times-equals operator for a cell index and an integer.
    */
   CellIndex<DIM>& operator*=(const int rhs);

   /**
    * Times operator for a cell index and an integer.
    */
   CellIndex<DIM> operator*(const int rhs) const;

   /**
    * Returns true if two cell index objects are equal.  All components
    * must be the same for equality.
    */
   bool operator==(const CellIndex<DIM>& rhs) const;

   /**
    * Returns true if two cell index objects are not equal.  Any of
    * the components may be different for inequality.
    */
   bool operator!=(const CellIndex<DIM>& rhs) const;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "CellIndex.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellIndex.C"
#endif
