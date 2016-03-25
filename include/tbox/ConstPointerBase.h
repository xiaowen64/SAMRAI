//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/ConstPointerBase.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: A smart const pointer base class with RTTI
//

#ifndef included_tbox_ConstPointerBase
#define included_tbox_ConstPointerBase

#include "SAMRAI_config.h"
#include "tbox/ReferenceCounter.h"
#include "tbox/DescribedClass.h"

namespace SAMRAI {
   namespace tbox {


/**
 * Class ConstPointerBase is an abstract base class used by template
 * class ConstPointer<TYPE> for type-safe conversion between various
 * pointer types.  It forms the base of the RTTI conversion hierarchy for
 * pointers.  Both the non-const pointer base class and the const pointer
 * class are subclasses of the const pointer base class.  This structure
 * ensures that RTTI conversion of pointers from const to const, non-const
 * to non-const, and non-const to const work as expected but that conversion
 * from const to non-const fail at compile-time.
 *
 * @see tbox::ConstPointer
 * @see tbox::PointerBase
 * @see tbox::Pointer
 */

class ConstPointerBase
{
public:
   ConstPointerBase();
   virtual ~ConstPointerBase();
   virtual ReferenceCounter *getSubclassReferenceCounter() const = 0;
   virtual const DescribedClass *getSubclassPointer() const = 0;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/ConstPointerBase.I"
#endif
#endif
