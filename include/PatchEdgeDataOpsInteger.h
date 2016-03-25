//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/edge/PatchEdgeDataOpsInteger.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Operations for integer edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsInteger
#define included_math_PatchEdgeDataOpsInteger

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchEdgeDataBasicOps.h"
#include "ArrayDataNormOpsInteger.h"
#include "Box.h"
#include "Patch.h"
#include "EdgeData.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchEdgeDataOpsInteger<DIM> provides a collection of operations
 * that may be used to manipulate integer edge-centered patch data.  The
 * operations include basic arithmetic, min, max, etc.  With the exception 
 * of a few basic routines, this class inherits its interface (and 
 * thus its functionality) from the base class PatchEdgeDataBasicOps<DIM>
 * from which it is derived.
 *
 * A more extensive set of operations is implemented for real (double and 
 * float) and complex patch data in the classes PatchEdgeDataOpsReal<DIM> 
 * and PatchEdgeDataOpsComplex<DIM>, repsectively. 
 *
 * @see math::PatchEdgeDataBasicOps
 */

template<int DIM>
class PatchEdgeDataOpsInteger : 
   public tbox::DescribedClass,
   public PatchEdgeDataBasicOps<DIM,int>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchEdgeDataOpsInteger();

   virtual ~PatchEdgeDataOpsInteger<DIM>();

   /**
    * Return the number of data values for the edge-centered data object
    * in the given box.  Note that it is assumed that the box refers to
    * the cell-centered index space corresponding to the patch hierarchy.
    */
   int numberOfEntries(const tbox::Pointer< pdat::EdgeData<DIM,int> >& data,
                       const hier::Box<DIM>& box) const;

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::EdgeData<DIM,int> >& dst,
                 const tbox::Pointer< pdat::EdgeData<DIM,int> >& src,
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
   void printData(const tbox::Pointer< pdat::EdgeData<DIM,int> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::EdgeData<DIM,int> >& dst,
                    const int& alpha,
                    const hier::Box<DIM>& box) const;

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void abs(tbox::Pointer< pdat::EdgeData<DIM,int> >& dst,
            const tbox::Pointer< pdat::EdgeData<DIM,int> >& src,
            const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchEdgeDataOpsInteger(const PatchEdgeDataOpsInteger<DIM>&);
   void operator=(const PatchEdgeDataOpsInteger<DIM>&);

   ArrayDataNormOpsInteger<DIM> d_array_ops;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchEdgeDataOpsInteger.C"
#endif
