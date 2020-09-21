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
#endif

namespace SAMRAI {
namespace tbox {

/**
 * UmpireAllocator a type alias for umpire::Allocator that enables API
 * consistency when SAMRAI is built with or without the Umpire library.
 * When Umpire is availalbe UmpireAllocator is an alias for umpire::Allocator,
 * so calling codes can pass in an umpire::Allocator anywhere that
 * UmpireAllocator is required in the SAMRAI API.  If Umpire is not
 * available, UmpireAllocator is an empty struct.
 */

#ifdef HAVE_UMPIRE

using UmpireAllocator = umpire::Allocator;

#else

struct UmpireAllocator
{
};

#endif

}
}

#endif
