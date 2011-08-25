/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer base class with RTTI
 *
 ************************************************************************/

#ifndef included_tbox_PointerBase
#define included_tbox_PointerBase

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/ConstPointerBase.h"
#include "SAMRAI/tbox/ReferenceCounter.h"

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

class PointerBase:public ConstPointerBase
{
public:
   PointerBase();
   virtual ~PointerBase();
   virtual ReferenceCounter *
   getSubclassReferenceCounter() const = 0;
   virtual const DescribedClass *
   getSubclassPointer() const = 0;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/PointerBase.I"
#endif
#endif
