//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/SAMRAI/hier/RefineOperator.I $
// Package:     SAMRAI data transfer
// Copyright:   (c) 1997-2009 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2861 $
// Modified:    $LastChangedDate: 2009-02-02 15:22:36 -0800 (Mon, 02 Feb 2009) $
// Description: Abstract base class for spatial refinement operators.
//

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const std::string&
RefineOperator::getOperatorName() const
{
   return d_name;
}

/*
 *************************************************************************
 * Register this in the static look-up table.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineOperator::registerInLookupTable(
   const std::string& name)
{
   s_lookup_table.insert(std::pair<std::string, RefineOperator *>(name, this));
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
RefineOperator::getDim() const
{
   return d_dim;
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineOperator::finalizeCallback()
{
   s_lookup_table.clear();
}

}
}
