//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mathops/edge/PatchEdgeDataOpsReal.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated operations for real edge-centered patch data.
//

#ifndef included_math_PatchEdgeDataOpsReal
#define included_math_PatchEdgeDataOpsReal

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_math_PatchEdgeDataBasicOps
#include "PatchEdgeDataBasicOps.h"
#endif
#ifndef included_math_PatchEdgeDataMiscellaneousOpsReal
#include "PatchEdgeDataMiscellaneousOpsReal.h"
#endif
#ifndef included_math_PatchEdgeDataNormOpsReal
#include "PatchEdgeDataNormOpsReal.h"
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
 * Class PatchEdgeDataOpsReal<DIM> provides a collection of operations
 * to manipulate float and double numerical edge-centered patch data.  The
 * operations include basic arithmetic, norms and ordering, and assorted 
 * miscellaneous operations.  With the exception of a few basic routines, 
 * this class inherits its interface (and thus its functionality) from the 
 * base classes PatchEdgeDataBasicOps<DIM>, PatchEdgeDataNormOpsReal<DIM>,
 * and PatchEdgeDataMiscellaneousOpsReal<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * edge-centered patch data operations that it provides.  
 *
 * Note that this templated class should only be used to instantiate 
 * objects with double or float as the template parameter.  A similar set of 
 * operations is implemented for complex and integer patch data in the classes 
 * PatchEdgeDataOpsComplex<DIM> and PatchEdgeDataOpsInteger<DIM>, 
 * repsectively. 
 *
 * @see math::PatchEdgeDataBasicOps
 * @see math::PatchEdgeDataMiscellaneousOpsReal
 * @see math::PatchEdgeDataNormOpsReal
 */

template<int DIM, class TYPE>
class PatchEdgeDataOpsReal  : 
   public tbox::DescribedClass,
   public PatchEdgeDataBasicOps<DIM,TYPE>,
   public PatchEdgeDataMiscellaneousOpsReal<DIM,TYPE>,
   public PatchEdgeDataNormOpsReal<DIM,TYPE>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchEdgeDataOpsReal();

   virtual ~PatchEdgeDataOpsReal<DIM,TYPE>();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
                 const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& src,
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
   void printData(const tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& data,
                  const hier::Box<DIM>& box,
                  std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::EdgeData<DIM,TYPE> >& dst,
                    const TYPE& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchEdgeDataOpsReal(const PatchEdgeDataOpsReal<DIM,TYPE>&);
   void operator=(const PatchEdgeDataOpsReal<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchEdgeDataOpsReal.C"
#endif
