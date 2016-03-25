//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/clustering/BoxGeneratorStrategy.h $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Strategy interface for box generation routines.
//
 
#ifndef included_mesh_BoxGeneratorStrategy
#define included_mesh_BoxGeneratorStrategy
 
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace mesh {

/**
 * Class BoxGeneratorStrategy<DIM> is an abstract base class that defines 
 * a Strategy pattern interface for operations to build boxes that cover a 
 * collection of tagged cells on a single AMR patch hierarchy level.
 * 
 * @see hier::PatchLevel
 */

template<int DIM> class BoxGeneratorStrategy  : public tbox::DescribedClass
{
public:

   /**
    * Default constructor.
    */
   BoxGeneratorStrategy();

   /**
    * Virtual destructor.
    */
   virtual ~BoxGeneratorStrategy<DIM>();

   /**
    * Create list of boxes whose union covers all integer tags on the patch 
    * level that match the specified tag value. Each box must be at least
    * as large as the given minimum size and the tolerances must be met.
    */
   virtual void findBoxesContainingTags(
      hier::BoxList<DIM>& boxes,
      const tbox::Pointer< hier::PatchLevel<DIM> > level,
      const int index,
      const int tag_val,
      const hier::Box<DIM>& bound_box,
      const hier::IntVector<DIM>& min_box,
      const double efficiency_tol,
      const double combine_tol) const = 0;

private:
   // The following are not implemented:
   BoxGeneratorStrategy(const BoxGeneratorStrategy<DIM>&);
   void operator=(const BoxGeneratorStrategy<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxGeneratorStrategy.C"
#endif
