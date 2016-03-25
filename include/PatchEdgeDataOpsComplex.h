//
// File:	PatchEdgeDataOpsComplex.h
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Operations for complex edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsComplex
#define included_math_PatchEdgeDataOpsComplex

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_math_PatchEdgeDataBasicOps
#include "PatchEdgeDataBasicOps.h"
#endif
#ifndef included_math_PatchEdgeDataNormOpsComplex
#include "PatchEdgeDataNormOpsComplex.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_pdat_EdgeData
#include "EdgeData.h"
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
                  ostream& s = tbox::plog) const;

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
