/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_CellVariable_C
#define included_pdat_CellVariable_C

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for cell variable objects
 *
 *************************************************************************
 */

template<class TYPE>
CellVariable<TYPE>::CellVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  boost::shared_ptr<hier::PatchDataFactory>(
                     new CellDataFactory<TYPE>(depth,
                                               hier::IntVector::getZero(dim)))) // default zero ghost cells
{
}

template<class TYPE>
CellVariable<TYPE>::~CellVariable()
{
}

template<class TYPE>
int CellVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<CellDataFactory<TYPE> > cell_factory =
      this->getPatchDataFactory();
   TBOX_ASSERT(cell_factory);
   return cell_factory->getDepth();
}

/*
 *************************************************************************
 *
 * These are private and should not be used.  They are defined here
 * because some template instantiation methods fail if some member
 * functions are left undefined.
 *
 *************************************************************************
 */

template<class TYPE>
CellVariable<TYPE>::CellVariable(
   const CellVariable<TYPE>& foo):
   hier::Variable(NULL,
                  boost::shared_ptr<hier::PatchDataFactory>(NULL))
{
   NULL_USE(foo);
}

template<class TYPE>
void CellVariable<TYPE>::operator = (
   const CellVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
