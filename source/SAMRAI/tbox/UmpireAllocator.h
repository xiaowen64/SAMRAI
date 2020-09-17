/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2020 Lawrence Livermore National Security, LLC
 * Description:   A wrapper for Umpire allocators.
 *
 ************************************************************************/

#ifndef included_tbox_UmpireAllocator
#define included_tbox_UmpireAllocator

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_UMPIRE
#include "umpire/Allocator.hpp"
#include "umpire/TypedAllocator.hpp"
#endif

namespace SAMRAI {
namespace tbox {

/**
 * Class UmpireAllocator is a wrapper class for umpire::Allocator that
 * can also serve as a no-op class when the library is built without Umpire,
 * so that code bases do not need to use ifdef blocks to differentiate calls
 * that use Umpire and calls that don't.
 */

class UmpireAllocator
{
public:
   /**
    */
   UmpireAllocator();

#ifdef HAVE_UMPIRE
   explicit UmpireAllocator(const umpire::Allocator& allocator);
#endif



   /**
    */
   virtual ~UmpireAllocator();

#ifdef HAVE_UMPIRE
   umpire::Allocator getAllocator() const
   {
      return d_allocator;
   }
#endif

private:
//   UmpireAllocator(
//      const UmpireAllocator&);             // not implemented
//   UmpireAllocator&
//   operator = (
//      const UmpireAllocator&);             // not implemented

#ifdef HAVE_UMPIRE
   umpire::Allocator d_allocator;
#endif

};

}
}

#endif
