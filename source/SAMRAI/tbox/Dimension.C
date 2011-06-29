/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer template class with RTTI 
 *
 ************************************************************************/

#ifndef included_tbox_Dimension_C
#define included_tbox_Dimension_C

#include "SAMRAI/tbox/Dimension.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/tbox/Dimension.I"
#endif

namespace SAMRAI {
namespace tbox {

std::ostream& operator << (
   std::ostream& s,
   const Dimension& dim)
{
   s << dim.getValue() << 'D';
   return s;
}

}
}

#endif
