/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#ifndef included_pdat_CellIndex
#define included_pdat_CellIndex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class CellIndex implements a simple n-dimensional integer
 * vector for cell centered variables.  Cell indices contain an integer
 * index location in AMR index space and are identical to the AMR indices.
 *
 * @see hier::Index
 * @see pdat::CellData
 * @see pdat::CellGeometry
 * @see pdat::CellIterator
 */

class CellIndex:public hier::Index
{
public:
   /**
    * The default constructor for a cell index creates an uninitialized index.
    */
   explicit CellIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a cell index from a regular AMR index.
    *
    */
   explicit CellIndex(
      const hier::Index& rhs);

   /**
    * The copy constructor creates a cell index equal to the argument.
    */
   CellIndex(
      const CellIndex& rhs);

   /**
    * The assignment operator sets the cell index equal to the argument.
    */
   CellIndex&
   operator = (
      const CellIndex& rhs);

   /**
    * The cell index destructor does nothing interesting.
    */
   ~CellIndex();

   /**
    * Plus-equals operator for a cell index and an integer vector.
    */
   CellIndex&
   operator += (
      const hier::IntVector& rhs);

   /**
    * Plus operator for a cell index and an integer vector.
    */
   CellIndex
   operator + (
      const hier::IntVector& rhs) const;

   /**
    * Plus-equals operator for a cell index and an integer.
    */
   CellIndex&
   operator += (
      const int rhs);

   /**
    * Plus operator for a cell index and an integer.
    */
   CellIndex
   operator + (
      const int rhs) const;

   /**
    * Minus-equals operator for a cell index and an integer vector.
    */
   CellIndex&
   operator -= (
      const hier::IntVector& rhs);

   /**
    * Minus operator for a cell index and an integer vector.
    */
   CellIndex
   operator - (
      const hier::IntVector& rhs) const;

   /**
    * Minus-equals operator for a cell index and an integer.
    */
   CellIndex&
   operator -= (
      const int rhs);

   /**
    * Minus operator for a cell index and an integer.
    */
   CellIndex
   operator - (
      const int rhs) const;

   /**
    * Times-equals operator for a cell index and an integer vector.
    */
   CellIndex&
   operator *= (
      const hier::IntVector& rhs);

   /**
    * Times operator for a cell index and an integer vector.
    */
   CellIndex
   operator * (
      const hier::IntVector& rhs) const;

   /**
    * Times-equals operator for a cell index and an integer.
    */
   CellIndex&
   operator *= (
      const int rhs);

   /**
    * Times operator for a cell index and an integer.
    */
   CellIndex
   operator * (
      const int rhs) const;

   /**
    * Returns true if two cell index objects are equal.  All components
    * must be the same for equality.
    */
   bool
   operator == (
      const CellIndex& rhs) const;

   /**
    * Returns true if two cell index objects are not equal.  Any of
    * the components may be different for inequality.
    */
   bool
   operator != (
      const CellIndex& rhs) const;
};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/pdat/CellIndex.I"
#endif
#endif
