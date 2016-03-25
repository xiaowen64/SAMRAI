/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart const pointer base class with RTTI
 *
 ************************************************************************/

#ifndef included_tbox_ConstPointerBase
#define included_tbox_ConstPointerBase

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/ReferenceCounter.h"
#include "SAMRAI/tbox/DescribedClass.h"

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
   virtual ReferenceCounter *
   getSubclassReferenceCounter() const = 0;
   virtual const DescribedClass *
   getSubclassPointer() const = 0;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/ConstPointerBase.I"
#endif
#endif
