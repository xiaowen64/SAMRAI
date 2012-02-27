/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management in AMR hierarchy
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const BoxContainer&
GridGeometry::getPhysicalDomain() const
{
   return d_physical_domain;
}

SAMRAI_INLINE_KEYWORD
bool
GridGeometry::getDomainIsSingleBox(
   const BlockId& block_id) const
{
   return d_domain_is_single_box[block_id.getBlockValue()];
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
GridGeometry::getDim() const {
   return d_dim;
}

SAMRAI_INLINE_KEYWORD
const BoxContainer&
GridGeometry::getDomainSearchTree() const
{
   return d_physical_domain;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

SAMRAI_INLINE_KEYWORD
const BoxContainer&
GridGeometry::getPeriodicDomainSearchTree() const
{
   return d_domain_with_images;
}

SAMRAI_INLINE_KEYWORD
const std::string&
GridGeometry::getObjectName() const {
   return d_object_name;
}

SAMRAI_INLINE_KEYWORD
const BoxContainer&
GridGeometry::Neighbor::getTransformedDomain() const
{
   return d_transformed_domain;
}

SAMRAI_INLINE_KEYWORD
int
GridGeometry::getNumberOfBlockSingularities() const
{
   return d_number_of_block_singularities;
}

SAMRAI_INLINE_KEYWORD
const BoxContainer&
GridGeometry::getSingularityBoxContainer(
   const BlockId& block_id) const
{
   return d_singularity[block_id.getBlockValue()];
}

SAMRAI_INLINE_KEYWORD
const std::vector<int>&
GridGeometry::getSingularityIndices(
   const BlockId& block_id) const
{
   return d_singularity_indices[block_id.getBlockValue()];
}

SAMRAI_INLINE_KEYWORD
bool
GridGeometry::reducedConnectivityExists(
   const BlockId& block_id) const
{
   return d_reduced_connect[block_id.getBlockValue()];
}

/*
 * ************************************************************************
 *
 * Returns the number of neighbors a specific block of the multiblock
 * domain has.
 *
 * ************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
GridGeometry::getNumberOfNeighbors(
   const BlockId& block_id) const
{
   return d_block_neighbors[block_id.getBlockValue()].getNumberOfItems();
}

SAMRAI_INLINE_KEYWORD
const tbox::List<GridGeometry::Neighbor>&
GridGeometry::getNeighbors(
   const BlockId& block_id) const
{
   return d_block_neighbors[block_id.getBlockValue()];
}

SAMRAI_INLINE_KEYWORD
bool
GridGeometry::Neighbor::isSingularity() const
{
   return d_is_singularity;
}

SAMRAI_INLINE_KEYWORD
const BlockId&
GridGeometry::Neighbor::getBlockId() const
{
   return d_block_id;
}

SAMRAI_INLINE_KEYWORD
const Transformation&
GridGeometry::Neighbor::getTransformation() const
{
   return d_transformation;
}

SAMRAI_INLINE_KEYWORD
Transformation::RotationIdentifier
GridGeometry::Neighbor::getRotationIdentifier() const
{
   return d_transformation.getRotation();
}

SAMRAI_INLINE_KEYWORD
const IntVector&
GridGeometry::Neighbor::getShift() const
{
   return d_transformation.getOffset();
}

SAMRAI_INLINE_KEYWORD
int
GridGeometry::getNumberBlocks() const
{
   return d_number_blocks;
}

SAMRAI_INLINE_KEYWORD
bool
GridGeometry::hasEnhancedConnectivity() const
{
   return d_has_enhanced_connectivity;
}

/*
 *************************************************************************
 *
 * Add operator to appropriate lookup list.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
GridGeometry::addCoarsenOperator(
   const boost::shared_ptr<CoarsenOperator>& coarsen_op)
{
   d_transfer_operator_registry->addCoarsenOperator(coarsen_op);
}

SAMRAI_INLINE_KEYWORD
void
GridGeometry::addRefineOperator(
   const boost::shared_ptr<RefineOperator>& refine_op)
{
   d_transfer_operator_registry->addRefineOperator(refine_op);
}
SAMRAI_INLINE_KEYWORD
void
GridGeometry::addTimeInterpolateOperator(
   const boost::shared_ptr<TimeInterpolateOperator>& time_op)
{
   d_transfer_operator_registry->addTimeInterpolateOperator(time_op);
}

/*
 *************************************************************************
 *
 * Search operator lists for operator matching request.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<CoarsenOperator>
GridGeometry::lookupCoarsenOperator(
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   return d_transfer_operator_registry->lookupCoarsenOperator(
      var, op_name);
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<RefineOperator>
GridGeometry::lookupRefineOperator(
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   return d_transfer_operator_registry->lookupRefineOperator(
      var, op_name);
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<TimeInterpolateOperator>
GridGeometry::lookupTimeInterpolateOperator(
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   return d_transfer_operator_registry->lookupTimeInterpolateOperator(
      var, op_name);
}

/*
 *************************************************************************
 * Compute the max operator stencil width from all constructed
 * refine and coarsen operators and the user-specified minimum value.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
IntVector
GridGeometry::getMaxTransferOpStencilWidth()
{
   return d_transfer_operator_registry->getMaxTransferOpStencilWidth();
}

/*
 *************************************************************************
 * Set the mininum value to be returned by
 * getMaxTransferOpStencilWidth().
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
GridGeometry::setMinTransferOpStencilWidth(
   const IntVector& min_width)
{
   d_transfer_operator_registry->setMinTransferOpStencilWidth(min_width);
}

}
}
