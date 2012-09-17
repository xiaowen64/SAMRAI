/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_SideVariable_C
#define included_pdat_SideVariable_C

#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for side variable objects
 *
 *************************************************************************
 */

template<class TYPE>
SideVariable<TYPE>::SideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   const hier::IntVector& directions,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  boost::make_shared<SideDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     directions)),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(directions)
{
   TBOX_ASSERT(directions.getDim() == getDim());
}

template<class TYPE>
SideVariable<TYPE>::~SideVariable()
{
}

template<class TYPE>
const hier::IntVector& SideVariable<TYPE>::getDirectionVector() const
{
   return d_directions;
}

template<class TYPE>
int SideVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<SideDataFactory<TYPE> > factory(
      getPatchDataFactory(),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
