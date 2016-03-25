//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/PointerBase.h $
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: A smart pointer base class with RTTI
//

#ifndef included_tbox_PointerBase
#define included_tbox_PointerBase

#include "SAMRAI_config.h"
#include "tbox/ConstPointerBase.h"
#include "tbox/ReferenceCounter.h"

namespace SAMRAI {
   namespace tbox {

/**
 * Class PointerBase is a base class used by template class
 * PointerRef<TYPE> for type-safe conversion between non-const
 * pointer types.  It is a subclass of ConstPointerBase.  Since
 * the non-const pointer class only takes this as a base class (and not
 * the const pointer base class), const pointers cannot be converted
 * into non-const pointers.
 *
 * @see tbox::ConstPointerBase
 * @see tbox::Pointer
 */

class PointerBase : public ConstPointerBase
{
public:
   PointerBase();
   virtual ~PointerBase();
   virtual ReferenceCounter *getSubclassReferenceCounter() const = 0;
   virtual const DescribedClass *getSubclassPointer() const = 0;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/PointerBase.I"
#endif
#endif
