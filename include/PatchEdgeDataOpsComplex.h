//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/edge/PatchEdgeDataOpsComplex.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Operations for complex edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsComplex
#define included_math_PatchEdgeDataOpsComplex

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "tbox/Complex.h"
#include "PatchEdgeDataBasicOps.h"
#include "PatchEdgeDataNormOpsComplex.h"
#include "Box.h"
#include "Patch.h"
#include "EdgeData.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchEdgeDataOpsComplex<DIM> provides a collection of operations
 * that may be used to manipulate complex edge-centered patch data.  The
 * operations include basic arithmetic and norms.  With the 
 * exception of a few basic routines, this class inherits its interface (and 
 * thus its functionality) from the base classes PatchEdgeDataBasicOps<DIM>, 
 * PatchEdgeDataNormOpsComplex<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * edge-centered patch data operations that it provides.  
 *
 * A similar set of operations is implemented for real (double and float) and 
 * integer patch data in the classes PatchEdgeDataOpsReal<DIM> and 
 * PatchEdgeDataOpsInteger<DIM>, repsectively. 
 *
 * @see math::PatchEdgeDataBasicOps
 * @see math::PatchEdgeDataNormOpsComplex
 */

template<int DIM>
class PatchEdgeDataOpsComplex : 
   public tbox::DescribedClass,
   public PatchEdgeDataBasicOps<DIM,dcomplex>,
   public PatchEdgeDataNormOpsComplex<DIM>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchEdgeDataOpsComplex();

   virtual ~PatchEdgeDataOpsComplex<DIM>();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& dst,
                 const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& src,
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
   void printData(const tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::EdgeData<DIM,dcomplex> >& dst,
                    const dcomplex& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchEdgeDataOpsComplex(const PatchEdgeDataOpsComplex<DIM>&);
   void operator=(const PatchEdgeDataOpsComplex<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchEdgeDataOpsComplex.C"
#endif
