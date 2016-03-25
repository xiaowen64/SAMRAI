//
// File:	PatchFactory.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Abstract factory class for creating patch classes
//

#ifndef included_hier_PatchFactory
#define included_hier_PatchFactory

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchDescriptor
#include "PatchDescriptor.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif


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
