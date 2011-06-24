/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Singleton registry for all tranfer operators.
 *
 ************************************************************************/

#ifndef included_hier_TransferOperatorRegistry
#define included_hier_TransferOperatorRegistry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/Pointer.h"

#include <string>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class TransferOperatorRegistry is intended to serve as the registry
 * for SAMRAI transfer operators.  It will be a singleton object held by class
 * hier::GridGeometry.
 *
 * This hier::TransferOperatorRegistry class provides a lookup mechanism to
 * search for time interpolation and spatial coarsening/refining operators.
 * That is, algorithms and applications that manage communication on an
 * AMR hierarchy may query the hier::TransferOperatorRegistry object for
 * operators that may be applied to specific variables.
 *
 * Typically, the operators are assigned to the hier::TrnasferOperatorRegistry
 * object in the constructor of the geometry object that defines the mesh
 * coordinate system.
 *
 * @note
 * <ol>
 *   <li> Operators are added to the heads of the operator lists, so the most
 *        recently added operator will be returned if more than one operator
 *        satisfies a given request.
 *   <li> Additional operators may be added to a transfer geometry object at
 *        any time during program execution. However, each operator must be
 *        added BEFORE it is requested or an unrecoverable assertion will be
 *        thrown and the program will abort.
 * </ol>
 *
 * See the time interpolation,spatial coarsening, and spatial refinement
 * operator base classes for more information about adding new operators for
 * either new patch data types or new operators for pre-existing patch data
 * types.
 *
 * @see hier::GridGeometry
 * @see hier::RefineOperator
 * @see hier::CoarsenOperator
 * @see hier::TimeInterpolateOperator
 */
class TransferOperatorRegistry
{
public:
   /*!
    * @brief Set the state of the hier::TransferOperatorRegistry members.
    *
    * @param[in]  dim
    */
  TransferOperatorRegistry(const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   virtual ~TransferOperatorRegistry();

   /*!
    * @brief Add a concrete spatial coarsening operator.
    *
    * @param[in]  coarsen_op The concrete coarsening operator to add to the
    *             lookup list.
    */
   void
   addCoarsenOperator(
      tbox::Pointer<CoarsenOperator> coarsen_op);

   /*!
    * @brief Add a concrete spatial refinement operator.
    *
    * @param[in]  refine_op The concrete refinement operator to add to the
    *             lookup list.
    */
   void
   addRefineOperator(
      tbox::Pointer<RefineOperator> refine_op);

   /*!
    * @brief Add a concrete time interpolation operator.
    *
    * @param[in]  time_op The concrete time interpolation operator to add
    *             to the lookup list.
    */
   void
   addTimeInterpolateOperator(
      tbox::Pointer<TimeInterpolateOperator> time_op);

   /*!
    * @brief Lookup function for coarsening operator.
    * 
    * Search list for the spatial coarsening operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    *
    * @param[in]     var The Variable for which the corresponding coarsening
    *                operator should match.
    * @param[in]     op_name The string identifier of the coarsening operator.
    */
   tbox::Pointer<CoarsenOperator>
      lookupCoarsenOperator(
         const tbox::Pointer<Variable>& var,
         const std::string& op_name);

   /*!
    * @brief Lookup function for refinement operator.
    * 
    * Search list for the spatial refinement operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    *
    * @param[in]     var The Variable for which the corresponding refinement
    *                operator should match.
    * @param[in]     op_name The string identifier of the refinement operator.
    */
   tbox::Pointer<RefineOperator>
      lookupRefineOperator(
         const tbox::Pointer<Variable>& var,
         const std::string& op_name);

   /*!
    * @brief Lookup function for time interpolation operator.
    * 
    * Search list for the time interpolation operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    *
    * @param[in]     var The Variable for which the corresponding time
    *                interpolation operator should match.
    * @param[in]     op_name The string identifier of the time interpolation
    *                operator.  \b Default: STD_LINEAR_TIME_INTERPOLATE
    */
   tbox::Pointer<TimeInterpolateOperator>
      lookupTimeInterpolateOperator(
         const tbox::Pointer<Variable>& var,
         const std::string& op_name =
            "STD_LINEAR_TIME_INTERPOLATE");

   /*!
    * @brief Get the max stencil width of all transfer operators.
    *
    * The max stencil width computed from all registered (constructed)
    * transfer operators.  This function is simply returns the max
    * from the registered refine, coarsen and time-refine operators.
    *
    * @return The max stencil width computed from all registered
    * operators.
    *
    * @see hier::GridGeometry::getMaxTransferOpStencilWidth().
    * @see hier::RefineOperator::getMaxRefineOpStencilWidth().
    * @see hier::CoarsenOperator::getMaxCoarsenOpStencilWidth().
    */
   IntVector
   getMaxTransferOpStencilWidth();

   /*!
    * @brief Set a minimum value on the value returned by
    * getMaxTransferOpStencilWidth().
    *
    * This method allows users to specify a minimum value returned by
    * getMaxTransferOpStencilWidth().  The default minimum is zero.
    * This value can be used as a substitute for data that is not yet
    * registered with the Geometry and therefore cannot be reflected
    * in getMaxTransferOpStencilWidth().
    *
    * @param[in]  min_value The minimum value to set.
    */
   void
   setMinTransferOpStencilWidth(
      const IntVector& min_value);

   /*!
    * @brief Get the dimension of the hier::GridGeometry holding this object.
    *
    * @return The dimension of the hier::GridGeometry holding this object.
    */
   const tbox::Dimension&
   getDim() const;

   /*!
    * @brief Print class data representation.
    *
    * @param[in]  os The std::ostream to print to.
    */
   void
   printClassData(
      std::ostream& os) const;

protected:
   /*!
    * @brief Builds the SAMRAI CoarsenOperator with name op_name associated
    * with the centering of variable var.  Can be specialized by applications
    * which introduce their own CoarsenOperators.  Specialized class will add
    * knowledge of how to build application's operators and will call
    * SAMRAITransferOperatorRegistry::buildCoarsenOperator to build a SAMRAI
    * libaray CoarsenOperator.
    */
   virtual tbox::Pointer<CoarsenOperator> buildCoarsenOperator(
      const tbox::Pointer<Variable>& var,
      const std::string& op_name) = 0;

   /*!
    * @brief Builds the SAMRAI RefineOperator with name op_name associated
    * with the centering of variable var.  Can be specialized by applications
    * which introduce their own RefineOperators.  Specialized class will add
    * knowledge of how to build application's operators and will call
    * SAMRAITransferOperatorRegistry::buildRefineOperator to build a SAMRAI
    * libaray RefineOperator.
    */
   virtual tbox::Pointer<RefineOperator> buildRefineOperator(
      const tbox::Pointer<Variable>& var,
      const std::string& op_name) = 0;

   /*!
    * @brief Builds the SAMRAI TimeInterpolateOperator with name op_name
    * associated with the centering of variable var.  Can be specialized by
    * applications which introduce their own TimeInterpolateOperators.
    * Specialized class will add knowledge of how to build application's
    * operators and will call
    * SAMRAITransferOperatorRegistry::buildTimeInterpolateOperator to build a
    * SAMRAI libaray TimeInterpolateOperator.
    */
   virtual tbox::Pointer<TimeInterpolateOperator> buildTimeInterpolateOperator(
      const tbox::Pointer<Variable>& var,
      const std::string& op_name) = 0;

private:
   /*
    * The list of spatial coarsening operators is maintained to lookup
    * operators for specific variables as requested by algorithms and/or
    * applications using the hier::GridGeometry holding this object.
    * Standard concrete coarsening operators can be found in the patchdata
    * package.  Additional operators may be added to this list at any time
    * (see addCoarsenOperator() function).
    */
   tbox::List<tbox::Pointer<CoarsenOperator> > d_coarsen_operators;

   /*
    * The list of spatial refinement operators is maintained to lookup
    * operators for specific variables as requested by algorithms and/or
    * applications using the hier::GridGeometry holding this object.
    * Standard concrete refinement operators can be found in the patchdata
    * package.  Additional operators may be added to this list at any time
    * (see addRefineOperator() function).
    */
   tbox::List<tbox::Pointer<RefineOperator> > d_refine_operators;

   /*
    * The list of time interpolation operators is maintained to lookup
    * operators for specific variables as requested by algorithms and/or
    * applications using the hier::GridGeometry holding this object.
    * Standard concrete time interpolation operators can be found in the
    * patchdata package.  Additional operators may be added to this list at
    * any time (see addTimeInterpolateOperator() function).
    */
   tbox::List<tbox::Pointer<TimeInterpolateOperator> > d_time_operators;

   /*!
    * @brief Value set by setMinTransferOpStencilWidth().
    */
   IntVector d_min_stencil_width;

   /*!
    * @brief The dimension of the GridGeometry holding this object.
    */
   tbox::Dimension d_dim;

   /*!
    * @brief true if a call to getMaxTransferOpStencilWidth has been made.
    */
   bool d_max_op_stencil_width_req;
   
};

}
}
#endif
