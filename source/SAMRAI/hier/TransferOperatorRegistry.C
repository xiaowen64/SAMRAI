/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Singleton registry for all tranfer operators.
 *
 ************************************************************************/

#ifndef included_hier_TransferOperatorRegistry_C
#define included_hier_TransferOperatorRegistry_C

#include "SAMRAI/hier/TransferOperatorRegistry.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *                                                                       *
 * Constructor and destructor for TransferOperatorRegistry objects.      *
 *                                                                       *
 *************************************************************************
 */

TransferOperatorRegistry::TransferOperatorRegistry(
   const tbox::Dimension& dim):
   d_min_stencil_width(dim, 0),
   d_dim(dim),
   d_max_op_stencil_width_req(false)
{
}

TransferOperatorRegistry::~TransferOperatorRegistry()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Add operator to appropriate lookup list.                              *
 *                                                                       *
 *************************************************************************
 */

void TransferOperatorRegistry::addCoarsenOperator(
   tbox::Pointer<CoarsenOperator> coarsen_op)
{
   if (d_max_op_stencil_width_req &&
       (coarsen_op->getStencilWidth() > getMaxTransferOpStencilWidth())) {
      TBOX_WARNING(
         "Adding coarsen operator " << coarsen_op->getOperatorName()
                                    << "\nwith stencil width greater than current maximum\n"
                                    << "after call to getMaxTransferOpStencilWidth.\n");
   }
   d_coarsen_operators.addItem(coarsen_op);
}

void TransferOperatorRegistry::addRefineOperator(
   tbox::Pointer<RefineOperator> refine_op)
{
   if (d_max_op_stencil_width_req &&
       (refine_op->getStencilWidth() > getMaxTransferOpStencilWidth())) {
      TBOX_WARNING(
         "Adding refine operator " << refine_op->getOperatorName()
                                   << "\nwith stencil width greater than current maximum\n"
                                   << "after call to getMaxTransferOpStencilWidth.\n");
   }
   d_refine_operators.addItem(refine_op);
}

void TransferOperatorRegistry::addTimeInterpolateOperator(
   tbox::Pointer<TimeInterpolateOperator> time_op)
{
   d_time_operators.addItem(time_op);
}

/*
 *************************************************************************
 *                                                                       *
 * Search operator lists for operator matching request.                  *
 *                                                                       *
 *************************************************************************
 */

tbox::Pointer<CoarsenOperator>
TransferOperatorRegistry::lookupCoarsenOperator(
   const tbox::Pointer<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(!var.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   tbox::Pointer<CoarsenOperator> coarsen_op(NULL);
   bool found_op = false;

   if ((op_name == "NO_COARSEN") ||
       (op_name == "USER_DEFINED_COARSEN") ||
       (op_name.empty())) {
      found_op = true;
   } else {

      tbox::List<tbox::Pointer<CoarsenOperator> >::Iterator
         lop = d_coarsen_operators.listStart();

      while (coarsen_op.isNull() && lop) {
         if (lop()->findCoarsenOperator(var, op_name)) {
            found_op = true;
            coarsen_op = lop();
         }
         lop++;
      }
   }

   if (!found_op) {
      coarsen_op = buildCoarsenOperator(var, op_name);
   }

   return coarsen_op;
}

tbox::Pointer<RefineOperator>
TransferOperatorRegistry::lookupRefineOperator(
   const tbox::Pointer<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(!var.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   tbox::Pointer<RefineOperator> refine_op(NULL);
   bool found_op = false;

   if ((op_name == "NO_REFINE") ||
       (op_name == "USER_DEFINED_REFINE") ||
       (op_name.empty())) {
      found_op = true;
   } else {

      tbox::List<tbox::Pointer<RefineOperator> >::Iterator
         lop = d_refine_operators.listStart();

      while (refine_op.isNull() && lop) {
         if (lop()->findRefineOperator(var, op_name)) {
            found_op = true;
            refine_op = lop();
         }
         lop++;
      }
   }

   if (!found_op) {
      refine_op = buildRefineOperator(var, op_name);
   }

   return refine_op;
}

tbox::Pointer<TimeInterpolateOperator>
TransferOperatorRegistry::lookupTimeInterpolateOperator(
   const tbox::Pointer<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(!var.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   tbox::Pointer<TimeInterpolateOperator> time_op(NULL);
   bool found_op = false;

   if ((op_name == "NO_TIME_INTERPOLATE") ||
       (op_name.empty())) {
      found_op = true;
   } else {

      tbox::List<tbox::Pointer<TimeInterpolateOperator> >::Iterator
         lop = d_time_operators.listStart();

      while (time_op.isNull() && lop) {
         if (lop()->findTimeInterpolateOperator(var, op_name)) {
            found_op = true;
            time_op = lop();
         }
         lop++;
      }
   }

   if (!found_op) {
      time_op = buildTimeInterpolateOperator(var, op_name);
   }

   return time_op;
}

/*
 *************************************************************************
 * Compute the max operator stencil width from all constructed
 * refine and coarsen operators and the user-specified minimum value.
 *************************************************************************
 */
IntVector
TransferOperatorRegistry::getMaxTransferOpStencilWidth()
{
   IntVector max_width(d_min_stencil_width);
   max_width.max(RefineOperator::getMaxRefineOpStencilWidth(getDim()));
   max_width.max(CoarsenOperator::getMaxCoarsenOpStencilWidth(getDim()));
   d_max_op_stencil_width_req = true;
   return max_width;
}

/*
 *************************************************************************
 * Set the mininum value to be returned by
 * getMaxTransferOpStencilWidth().
 *************************************************************************
 */

void
TransferOperatorRegistry::setMinTransferOpStencilWidth(
   const IntVector& min_width)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_min_stencil_width, min_width);
   d_min_stencil_width = min_width;
}

/*
 *************************************************************************
 *******Get the dimension of the hier::GridGeometry object holding this singleton.
 *************************************************************************
 */
const tbox::Dimension&
TransferOperatorRegistry::getDim() const
{
   return d_dim;
}

/*
 *************************************************************************
 *                                                                       *
 * Print CartesianGridGeometry class data.                               *
 *                                                                       *
 *************************************************************************
 */

void
TransferOperatorRegistry::printClassData(
   std::ostream& os) const
{
   os << "printing TransferOperatorRegistry data..." << std::endl;
   os << "TransferOperatorRegistry: this = "
      << (TransferOperatorRegistry *)this << std::endl;

   os << "Coarsen operator list: " << std::endl;
   tbox::List<tbox::Pointer<CoarsenOperator> >::Iterator
      cop = d_coarsen_operators.listStart();
   while (cop) {
      os << (CoarsenOperator *)cop() << std::endl;
      cop++;
   }

   os << "Refine operator list: " << std::endl;
   tbox::List<tbox::Pointer<RefineOperator> >::Iterator
      rop = d_refine_operators.listStart();
   while (rop) {
      os << (RefineOperator *)rop() << std::endl;
      rop++;
   }

   os << "Time interpolate operator list: " << std::endl;
   tbox::List<tbox::Pointer<TimeInterpolateOperator> >::Iterator
      top = d_time_operators.listStart();
   while (top) {
      os << (TimeInterpolateOperator *)top() << std::endl;
      top++;
   }
}

}
}
#endif
