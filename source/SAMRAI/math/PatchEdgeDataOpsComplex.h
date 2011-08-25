/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex edge-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataOpsComplex
#define included_math_PatchEdgeDataOpsComplex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/math/PatchEdgeDataBasicOps.h"
#include "SAMRAI/math/PatchEdgeDataNormOpsComplex.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"

#include <iostream>

namespace SAMRAI {
namespace math {

/**
 * Class PatchEdgeDataOpsComplex provides a collection of operations
 * that may be used to manipulate complex edge-centered patch data.  The
 * operations include basic arithmetic and norms.  With the
 * assertion of a few basic routines, this class inherits its interface (and
 * thus its functionality) from the base classes PatchEdgeDataBasicOps,
 * PatchEdgeDataNormOpsComplex from which it is derived.  The
 * name of each of these base classes is indicative of the set of
 * edge-centered patch data operations that it provides.
 *
 * A similar set of operations is implemented for real (double and float) and
 * integer patch data in the classes PatchEdgeDataOpsReal and
 * PatchEdgeDataOpsInteger, repsectively.
 *
 * @see math::PatchEdgeDataBasicOps
 * @see math::PatchEdgeDataNormOpsComplex
 */

class PatchEdgeDataOpsComplex:
   public tbox::DescribedClass,
   public PatchEdgeDataBasicOps<dcomplex>,
   public PatchEdgeDataNormOpsComplex
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchEdgeDataOpsComplex();

   virtual ~PatchEdgeDataOpsComplex();

   /**
    * Copy dst data to src data over given box.
    */
   void
   copyData(
      tbox::Pointer<pdat::EdgeData<dcomplex> >& dst,
      const tbox::Pointer<pdat::EdgeData<dcomplex> >& src,
      const hier::Box& box) const;

   /**
    * Swap pointers for patch data objects.  Objects are checked for
    * consistency of depth, box, and ghost box.
    */
   void
   swapData(
      tbox::Pointer<hier::Patch> patch,
      const int data1_id,
      const int data2_id) const;

   /**
    * Print data entries over given box to given output stream.
    */
   void
   printData(
      const tbox::Pointer<pdat::EdgeData<dcomplex> >& data,
      const hier::Box& box,
      std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void
   setToScalar(
      tbox::Pointer<pdat::EdgeData<dcomplex> >& dst,
      const dcomplex& alpha,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   PatchEdgeDataOpsComplex(
      const PatchEdgeDataOpsComplex&);
   void
   operator = (
      const PatchEdgeDataOpsComplex&);

};

}
}
#endif
