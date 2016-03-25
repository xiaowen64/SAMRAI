//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/cell/PatchCellDataOpsInteger.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Operations for integer cell-centered patch data.
//

#ifndef included_math_PatchCellDataOpsInteger
#define included_math_PatchCellDataOpsInteger

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_math_PatchCellDataBasicOps
#include "PatchCellDataBasicOps.h"
#endif
#ifndef included_math_ArrayDataNormOpsInteger
#include "ArrayDataNormOpsInteger.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_pdat_CellData
#include "CellData.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif


namespace SAMRAI {
    namespace math {

/**
 * Class PatchCellDataOpsInteger<DIM> provides a collection of operations
 * that may be used to manipulate integer cell-centered patch data.  The
 * operations include basic arithmetic, min, max, etc.  With the exception 
 * of a few basic routines, this class inherits its interface (and 
 * thus its functionality) from the base class PatchCellDataBasicOps<DIM>
 * from which it is derived.
 *
 * A more extensive set of operations is implemented for real (double and 
 * float) and complex patch data in the classes PatchCellDataOpsReal<DIM> 
 * and PatchCellDataOpsComplex<DIM>, respectively. 
 *
 * @see math::PatchCellDataBasicOps
 */

template<int DIM>
class PatchCellDataOpsInteger : 
   public tbox::DescribedClass,
   public PatchCellDataBasicOps<DIM,int>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchCellDataOpsInteger();

   virtual ~PatchCellDataOpsInteger<DIM>();

   /**
    * Return the number of data values for the cell-centered data object
    * in the given box.
    */
   int numberOfEntries(const tbox::Pointer< pdat::CellData<DIM,int> >& data,
                       const hier::Box<DIM>& box) const;

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::CellData<DIM,int> >& dst,
                 const tbox::Pointer< pdat::CellData<DIM,int> >& src,
                 const hier::Box<DIM>& box) const;

   /**
    * Swap pointers for patch data objects.  Objects are checked for 
    * consistency of depth, box, and ghost box.
    */
   void swapData(tbox::Pointer< hier::Patch<DIM> > patch,
                 const int data1_id,
                 const int data2_id) const;

   /**
    * Print data entries over given box to given output stream.
    */
   void printData(const tbox::Pointer< pdat::CellData<DIM,int> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::CellData<DIM,int> >& dst,
                    const int& alpha,
                    const hier::Box<DIM>& box) const;

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void abs(tbox::Pointer< pdat::CellData<DIM,int> >& dst,
            const tbox::Pointer< pdat::CellData<DIM,int> >& src,
            const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchCellDataOpsInteger(const PatchCellDataOpsInteger<DIM>&);
   void operator=(const PatchCellDataOpsInteger<DIM>&);

   ArrayDataNormOpsInteger<DIM> d_array_ops;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchCellDataOpsInteger.C"
#endif
