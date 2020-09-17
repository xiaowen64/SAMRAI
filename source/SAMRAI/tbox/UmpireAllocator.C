/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2020 Lawrence Livermore National Security, LLC
 * Description:   A wrapper for Umpire allocators
 *
 ************************************************************************/

#include "SAMRAI/tbox/UmpireAllocator.h"

namespace SAMRAI {
namespace tbox {

UmpireAllocator::UmpireAllocator()
{
}

#ifdef HAVE_UMPIRE
UmpireAllocator::UmpireAllocator(const umpire::Allocator& allocator)
: d_allocator(allocator)
{
}
#endif

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

UmpireAllocator::~UmpireAllocator()
{
}


}
}
