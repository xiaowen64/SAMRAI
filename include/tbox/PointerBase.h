//
// File:	PointerBase.h
// Package:	SAMRAI toolbox for memory management
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A smart pointer base class with RTTI
//

#ifndef included_tbox_PointerBase
#define included_tbox_PointerBase

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_ConstPointerBase
#include "tbox/ConstPointerBase.h"
#endif
#ifndef included_tbox_ReferenceCounter
#include "tbox/ReferenceCounter.h"
#endif

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
