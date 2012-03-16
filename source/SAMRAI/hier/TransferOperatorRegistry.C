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

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * Constructor and destructor for TransferOperatorRegistry objects.
 *
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
 *
 * Add operator to appropriate lookup list.
 *
 *************************************************************************
 */

void
TransferOperatorRegistry::addCoarsenOperator(
   const boost::shared_ptr<CoarsenOperator>& coarsen_op)
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

void
TransferOperatorRegistry::addRefineOperator(
   const boost::shared_ptr<RefineOperator>& refine_op)
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

/*
 *************************************************************************
 *
 * Search operator lists for operator matching request.
 *
 *************************************************************************
 */

boost::shared_ptr<CoarsenOperator>
TransferOperatorRegistry::lookupCoarsenOperator(
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(var);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   boost::shared_ptr<CoarsenOperator> coarsen_op;
   bool found_op = false;

   if ((op_name == "NO_COARSEN") ||
       (op_name == "USER_DEFINED_COARSEN") ||
       (op_name.empty())) {
      found_op = true;
   } else {

      tbox::List<boost::shared_ptr<CoarsenOperator> >::Iterator lop =
         d_coarsen_operators.listStart();

      while (!coarsen_op && lop) {
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

boost::shared_ptr<RefineOperator>
TransferOperatorRegistry::lookupRefineOperator(
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(var);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   boost::shared_ptr<RefineOperator> refine_op;
   bool found_op = false;

   if ((op_name == "NO_REFINE") ||
       (op_name == "USER_DEFINED_REFINE") ||
       (op_name.empty())) {
      found_op = true;
   } else {

      tbox::List<boost::shared_ptr<RefineOperator> >::Iterator lop =
         d_refine_operators.listStart();

      while (!refine_op && lop) {
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

boost::shared_ptr<TimeInterpolateOperator>
TransferOperatorRegistry::lookupTimeInterpolateOperator(
   const boost::shared_ptr<Variable>& var,
   const std::string& op_name)
{
   TBOX_ASSERT(var);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   boost::shared_ptr<TimeInterpolateOperator> time_op;
   bool found_op = false;

   if ((op_name == "NO_TIME_INTERPOLATE") ||
       (op_name.empty())) {
      found_op = true;
   } else {

      tbox::List<boost::shared_ptr<TimeInterpolateOperator> >::Iterator lop =
         d_time_operators.listStart();

      while (!time_op && lop) {
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
 *
 * Print CartesianGridGeometry class data.
 *
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
   tbox::List<boost::shared_ptr<CoarsenOperator> >::Iterator cop =
      d_coarsen_operators.listStart();
   while (cop) {
      os << cop().get() << std::endl;
      cop++;
   }

   os << "Refine operator list: " << std::endl;
   tbox::List<boost::shared_ptr<RefineOperator> >::Iterator rop =
      d_refine_operators.listStart();
   while (rop) {
      os << rop().get() << std::endl;
      rop++;
   }

   os << "Time interpolate operator list: " << std::endl;
   tbox::List<boost::shared_ptr<TimeInterpolateOperator> >::Iterator top =
      d_time_operators.listStart();
   while (top) {
      os << top().get() << std::endl;
      top++;
   }
}

}
}
#endif
