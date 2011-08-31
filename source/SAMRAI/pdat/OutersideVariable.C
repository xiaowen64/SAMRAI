/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OutersideVariable_C
#define included_pdat_OutersideVariable_C

#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/OutersideDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *                                                                       *
 * Constructor and destructor for side variable objects                  *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
OutersideVariable<TYPE>::OutersideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  tbox::Pointer<SAMRAI::hier::PatchDataFactory>(new
                                                                OutersideDataFactory
                                                                <
                                                                   TYPE>(dim,
                                                                         depth)))
{
}

template<class TYPE>
OutersideVariable<TYPE>::~OutersideVariable()
{
}

template<class TYPE>
int OutersideVariable<TYPE>::getDepth() const
{
   tbox::Pointer<OutersideDataFactory<TYPE> > factory =
      this->getPatchDataFactory();
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

/*
 *************************************************************************
 *                                                                       *
 * These are private and should not be used.  They are defined here      *
 * because some template instantiation methods fail if some member       *
 * functions are left undefined.                                         *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
OutersideVariable<TYPE>::OutersideVariable(
   const OutersideVariable<TYPE>& foo):
   hier::Variable(NULL, tbox::Pointer<SAMRAI::hier::PatchDataFactory>(NULL))
{
   NULL_USE(foo);
}

template<class TYPE>
void OutersideVariable<TYPE>::operator = (
   const OutersideVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
