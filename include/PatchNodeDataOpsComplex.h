//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/node/PatchNodeDataOpsComplex.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Operations for complex node-centered patch data.
//

#ifndef included_math_PatchNodeDataOpsComplex
#define included_math_PatchNodeDataOpsComplex

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "tbox/Complex.h"
#include "PatchNodeDataBasicOps.h"
#include "PatchNodeDataNormOpsComplex.h"
#include "Box.h"
#include "Patch.h"
#include "NodeData.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchNodeDataOpsComplex<DIM> provides a collection of operations
 * that may be used to manipulate complex node-centered patch data.  The
 * operations include basic arithmetic and norms.  With the 
 * exception of a few basic routines, this class inherits its interface (and 
 * thus its functionality) from the base classes PatchNodeDataBasicOps<DIM>, 
 * PatchNodeDataNormOpsComplex<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * node-centered patch data operations that it provides.  
 *
 * A similar set of operations is implemented for real (double and float) and 
 * integer patch data in the classes PatchNodeDataOpsReal<DIM> and 
 * PatchNodeDataOpsInteger<DIM>, repsectively. 
 *
 * @see math::PatchNodeDataBasicOps
 * @see math::PatchNodeDataNormOpsComplex
 */

template<int DIM>
class PatchNodeDataOpsComplex  : 
   public tbox::DescribedClass,
   public PatchNodeDataBasicOps<DIM,dcomplex>,
   public PatchNodeDataNormOpsComplex<DIM>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchNodeDataOpsComplex();

   virtual ~PatchNodeDataOpsComplex<DIM>();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& dst,
                 const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& src,
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
   void printData(const tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::NodeData<DIM,dcomplex> >& dst,
                    const dcomplex& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchNodeDataOpsComplex(const PatchNodeDataOpsComplex<DIM>&);
   void operator=(const PatchNodeDataOpsComplex<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchNodeDataOpsComplex.C"
#endif
