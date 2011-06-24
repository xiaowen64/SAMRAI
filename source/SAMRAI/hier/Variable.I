/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Base class for application-level variables 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
int Variable::getCurrentMaximumInstanceNumber()
{
   return s_instance_counter;
}

SAMRAI_INLINE_KEYWORD
int Variable::getInstanceIdentifier() const
{
   return d_instance;
}

SAMRAI_INLINE_KEYWORD
const std::string& Variable::getName() const
{
   return d_name;
}

SAMRAI_INLINE_KEYWORD
void Variable::setPatchDataFactory(
   tbox::Pointer<PatchDataFactory> factory)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *factory);
   d_factory = factory;
}

SAMRAI_INLINE_KEYWORD
tbox::Pointer<PatchDataFactory>
Variable::getPatchDataFactory() const
{
   return d_factory;
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& Variable::getDim() const
{
   return d_dim;
}

}
}
