//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/face/PatchFaceDataOpsReal.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated operations for real face-centered patch data.
//

#ifndef included_math_PatchFaceDataOpsReal
#define included_math_PatchFaceDataOpsReal

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchFaceDataBasicOps.h"
#include "PatchFaceDataMiscellaneousOpsReal.h"
#include "PatchFaceDataNormOpsReal.h"
#include "Box.h"
#include "Patch.h"
#include "FaceData.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace math {

/**
 * Class PatchFaceDataOpsReal<DIM> provides a collection of operations
 * to manipulate float and double numerical face-centered patch data.  The
 * operations include basic arithmetic, norms and ordering, and assorted 
 * miscellaneous operations.  With the exception of a few basic routines, 
 * this class inherits its interface (and thus its functionality) from the 
 * base classes PatchFaceDataBasicOps<DIM>, PatchFaceDataNormOpsReal<DIM>,
 * and PatchFaceDataMiscellaneousOpsReal<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * face-centered patch data operations that it provides.  
 *
 * Note that this templated class should only be used to instantiate 
 * objects with double or float as the template parameter.  A similar set of 
 * operations is implemented for complex and integer patch data in the classes 
 * PatchFaceDataOpsComplex<DIM> and PatchFaceDataOpsInteger<DIM>, 
 * repsectively. 
 *
 * @see math::PatchFaceDataBasicOps
 * @see math::PatchFaceDataMiscellaneousOpsReal
 * @see math::PatchFaceDataNormOpsReal
 */

template<int DIM, class TYPE>
class PatchFaceDataOpsReal : 
   public tbox::DescribedClass,
   public PatchFaceDataBasicOps<DIM,TYPE>,
   public PatchFaceDataMiscellaneousOpsReal<DIM,TYPE>,
   public PatchFaceDataNormOpsReal<DIM,TYPE>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchFaceDataOpsReal();

   virtual ~PatchFaceDataOpsReal<DIM,TYPE>();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                 const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& src,
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
   void printData(const tbox::Pointer< pdat::FaceData<DIM,TYPE> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::FaceData<DIM,TYPE> >& dst,
                    const TYPE& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchFaceDataOpsReal(const PatchFaceDataOpsReal<DIM,TYPE>&);
   void operator=(const PatchFaceDataOpsReal<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchFaceDataOpsReal.C"
#endif
