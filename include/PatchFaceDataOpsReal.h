//
// File:	PatchFaceDataOpsReal.h
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated operations for real face-centered patch data.
//

#ifndef included_math_PatchFaceDataOpsReal
#define included_math_PatchFaceDataOpsReal

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_math_PatchFaceDataBasicOps
#include "PatchFaceDataBasicOps.h"
#endif
#ifndef included_math_PatchFaceDataMiscellaneousOpsReal
#include "PatchFaceDataMiscellaneousOpsReal.h"
#endif
#ifndef included_math_PatchFaceDataNormOpsReal
#include "PatchFaceDataNormOpsReal.h"
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
                  ostream& s = tbox::plog) const;

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
