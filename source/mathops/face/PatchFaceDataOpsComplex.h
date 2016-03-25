//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/face/PatchFaceDataOpsComplex.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Operations for complex face-centered patch data.
//

#ifndef included_math_PatchFaceDataOpsComplex
#define included_math_PatchFaceDataOpsComplex

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_math_PatchFaceDataBasicOps
#include "PatchFaceDataBasicOps.h"
#endif
#ifndef included_math_PatchFaceDataNormOpsComplex
#include "PatchFaceDataNormOpsComplex.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_pdat_FaceData
#include "FaceData.h"
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
 * Class PatchFaceDataOpsComplex<DIM> provides a collection of operations
 * that may be used to manipulate complex face-centered patch data.  The
 * operations include basic arithmetic and norms.  With the 
 * exception of a few basic routines, this class inherits its interface (and 
 * thus its functionality) from the base classes PatchFaceDataBasicOps<DIM>, 
 * PatchFaceDataNormOpsComplex<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * face-centered patch data operations that it provides.  
 *
 * A similar set of operations is implemented for real (double and float) and 
 * integer patch data in the classes PatchFaceDataOpsReal<DIM> and 
 * PatchFaceDataOpsInteger<DIM>, repsectively. 
 *
 * @see math::PatchFaceDataBasicOps
 * @see math::PatchFaceDataNormOpsComplex
 */

template<int DIM> class PatchFaceDataOpsComplex  : 
   public tbox::DescribedClass,
   public PatchFaceDataBasicOps<DIM,dcomplex>,
   public PatchFaceDataNormOpsComplex<DIM>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchFaceDataOpsComplex();

   virtual ~PatchFaceDataOpsComplex<DIM>();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& dst,
                 const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& src,
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
   void printData(const tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::FaceData<DIM,dcomplex> >& dst,
                    const dcomplex& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchFaceDataOpsComplex(const PatchFaceDataOpsComplex<DIM>&);
   void operator=(const PatchFaceDataOpsComplex<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchFaceDataOpsComplex.C"
#endif
