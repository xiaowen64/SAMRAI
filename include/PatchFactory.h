//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/PatchFactory.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Abstract factory class for creating patch classes
//

#ifndef included_hier_PatchFactory
#define included_hier_PatchFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "tbox/Pointer.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"


namespace SAMRAI {
    namespace hier {

/**
 * Class PatchFactory<DIM> is a factory object used to create new patches.
 * New types of patch objects can be introduced into the hierarchy through
 * derivation and re-defining the allocate member function.  There should
 * be no direct calls to the patch constructor (other than through the
 * patch factory).
 *
 * @see hier::Patch
 */

template<int DIM> class PatchFactory  : public tbox::DescribedClass
{
public:
   /**
    * Construct a patch factory object.
    */
   PatchFactory();

   /**
    * Virtual destructor for patch factory objects.
    */
   virtual ~PatchFactory<DIM>();

   /**
    * Allocate a patch with the specified domain and patch descriptor.
    */
   virtual tbox::Pointer< Patch<DIM> >
   allocate(const Box<DIM>& box,
            tbox::Pointer< PatchDescriptor<DIM> > descriptor) const;

private:
   PatchFactory(const PatchFactory<DIM>&);	// not implemented
   void operator=(const PatchFactory<DIM>&);		// not implemented

};

}
}
#ifndef DEBUG_NO_INLINE
#include "PatchFactory.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchFactory.C"
#endif
