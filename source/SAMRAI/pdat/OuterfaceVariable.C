/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceVariable_C
#define included_pdat_OuterfaceVariable_C

#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OuterfaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *									*
 * Constructor and destructor for face variable objects			*
 *									*
 *************************************************************************
 */

template<class TYPE>
OuterfaceVariable<TYPE>::OuterfaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  tbox::Pointer<hier::PatchDataFactory>(new
                                                        OuterfaceDataFactory<
                                                           TYPE>(dim,
                                                                 depth)))
{
}

template<class TYPE>
OuterfaceVariable<TYPE>::~OuterfaceVariable()
{
}

template<class TYPE>
int OuterfaceVariable<TYPE>::getDepth() const
{
   tbox::Pointer<OuterfaceDataFactory<TYPE> > factory =
      this->getPatchDataFactory();
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

/*
 *************************************************************************
 *									*
 * These are private and should not be used.  They are defined here	*
 * because some template instantiation methods fail if some member	*
 * functions are left undefined.						*
 *									*
 *************************************************************************
 */

template<class TYPE>
OuterfaceVariable<TYPE>::OuterfaceVariable(
   const OuterfaceVariable<TYPE>& foo):
   hier::Variable(NULL, tbox::Pointer<SAMRAI::hier::PatchDataFactory>(NULL))
{
   NULL_USE(foo);
}

template<class TYPE>
void OuterfaceVariable<TYPE>::operator = (
   const OuterfaceVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
