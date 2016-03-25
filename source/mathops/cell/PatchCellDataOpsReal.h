//
// File:	PatchCellDataOpsReal.h
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Templated operations for real cell-centered patch data.
//

#ifndef included_math_PatchCellDataOpsReal
#define included_math_PatchCellDataOpsReal

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_math_PatchCellDataBasicOps
#include "PatchCellDataBasicOps.h"
#endif
#ifndef included_math_PatchCellDataMiscellaneousOpsReal
#include "PatchCellDataMiscellaneousOpsReal.h"
#endif
#ifndef included_math_PatchCellDataNormOpsReal
#include "PatchCellDataNormOpsReal.h"
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
 * Class PatchCellDataOpsReal<DIM> provides a collection of operations
 * to manipulate float and double numerical cell-centered patch data.  The
 * operations include basic arithmetic, norms and ordering, and assorted 
 * miscellaneous operations.  With the exception of a few basic routines, 
 * this class inherits its interface (and thus its functionality) from the 
 * base classes PatchCellDataBasicOps<DIM>, PatchCellDataNormOpsReal<DIM>,
 * and PatchCellDataMiscellaneousOpsReal<DIM> from which it is derived.  The 
 * name of each of these base classes is indicative of the set of 
 * cell-centered patch data operations that it provides.  
 *
 * Note that this templated class should only be used to instantiate 
 * objects with double or float as the template parameter.  A similar set of 
 * operations is implemented for complex and integer patch data in the classes 
 * PatchCellDataOpsComplex<DIM> and PatchCellDataOpsInteger<DIM>, 
 * respectively. 
 *
 * @see math::PatchCellDataBasicOps
 * @see math::PatchCellDataMiscellaneousOpsReal
 * @see math::PatchCellDataNormOpsReal
 */

template<int DIM, class TYPE>
class PatchCellDataOpsReal : 
   public tbox::DescribedClass,
   public PatchCellDataBasicOps<DIM,TYPE>,
   public PatchCellDataMiscellaneousOpsReal<DIM,TYPE>,
   public PatchCellDataNormOpsReal<DIM,TYPE>
{
public:
   /** 
    * Empty constructor and destructor.
    */
   PatchCellDataOpsReal();

   virtual ~PatchCellDataOpsReal<DIM,TYPE>();

   /**
    * Copy dst data to src data over given box.
    */
   void copyData(tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
                 const tbox::Pointer< pdat::CellData<DIM,TYPE> >& src,
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
   void printData(const tbox::Pointer< pdat::CellData<DIM,TYPE> >& data,
                  const hier::Box<DIM>& box,
                  ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void setToScalar(tbox::Pointer< pdat::CellData<DIM,TYPE> >& dst,
                    const TYPE& alpha,
                    const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   PatchCellDataOpsReal(const PatchCellDataOpsReal<DIM,TYPE>&);
   void operator=(const PatchCellDataOpsReal<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchCellDataOpsReal.C"
#endif
