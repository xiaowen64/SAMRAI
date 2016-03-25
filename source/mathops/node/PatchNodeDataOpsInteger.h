//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/node/PatchNodeDataOpsInteger.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Operations for integer node-centered patch data.
//

#ifndef included_math_PatchNodeDataOpsInteger
#define included_math_PatchNodeDataOpsInteger

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchNodeDataBasicOps.h"
#include "ArrayDataNormOpsInteger.h"
#include "Box.h"
#include "Patch.h"
#include "NodeData.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchNodeDataOpsInteger<DIM> provides a collection of operations
 * that may be used to manipulate integer node-centered patch data.  The
 * operations include basic arithmetic, min, max, etc.  With the exception 
 * of a few basic routines, this class inherits its interface (and 
 * thus its functionality) from the base class PatchNodeDataBasicOps<DIM>
 * from which it is derived.
 *
 * A more extensive set of operations is implemented for real (double and 
 * float) and complex patch data in the classes PatchNodeDataOpsReal<DIM> 
 * and PatchNodeDataOpsComplex<DIM>, repsectively. 
 *
 * @see math::PatchNodeDataBasicOps
 */

template<int DIM>
class PatchNodeDataOpsInteger  : 
   public tbox::DescribedClass,
   public PatchNodeDataBasicOps<DIM,int>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchNodeDataOpsInteger();

   virtual ~PatchNodeDataOpsInteger<DIM>();

   /**
    * Return the number of data values for the node-centered data object
    * in the given box.  Note that it is assumed that the box refers to
    * the cell-centered index space corresponding to the patch hierarchy.
    */
   int numberOfEntries(const tbox::Pointer< pdat::NodeData<DIM,int> >& data,
                       const hier::Box<DIM>& box) const;

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::NodeData<DIM,int> >& dst,
                 const tbox::Pointer< pdat::NodeData<DIM,int> >& src,
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
   void printData(const tbox::Pointer< pdat::NodeData<DIM,int> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::NodeData<DIM,int> >& dst,
                    const int& alpha,
                    const hier::Box<DIM>& box) const;

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void abs(tbox::Pointer< pdat::NodeData<DIM,int> >& dst,
            const tbox::Pointer< pdat::NodeData<DIM,int> >& src,
            const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchNodeDataOpsInteger(const PatchNodeDataOpsInteger<DIM>&);
   void operator=(const PatchNodeDataOpsInteger<DIM>&);

   ArrayDataNormOpsInteger<DIM> d_array_ops;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchNodeDataOpsInteger.C"
#endif
