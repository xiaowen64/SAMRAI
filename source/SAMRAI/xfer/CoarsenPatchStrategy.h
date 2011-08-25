/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenPatchStrategy
#define included_xfer_CoarsenPatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/DescribedClass.h"

#include <set>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class for user-defined patch data coarsening operations.
 *
 * CoarsenPatchStrategy provides an interface for a user to supply
 * methods for application-specific coarsening of data between two levels in
 * an AMR patch hierarchy.  A concrete subclass must define three member
 * functions to perform the following tasks:
 *
 * <ul>
 *    <li> define maximum stencil width for user-defined coarsen operations
 *    <li> preprocess the coarsening
 *    <li> postprocess the coarsening
 * </ul>
 *
 * Note that the preprocess member function is called before standard data
 * coarsening using CoarsenOperators and the postprocess member function is
 * called afterwards.
 *
 * @see xfer::CoarsenAlgorithm
 * @see xfer::CoarsenSchedule
 */

class CoarsenPatchStrategy:
   public virtual tbox::DescribedClass
{
public:
   /*!
    * @brief Get the maximum stencil width over all CoarsenPatchStrategy objects
    * used in an application.
    *
    * @return The maximum of the return values of calls to
    * getCoarsenOpStencilWidth() for every CoarsenPatchStrategy of the
    * given Dimension used in an application.
    *
    * @param[in] dim   Only objects with this dimension will be used to
    *                  calculate the max.  If a CoarsenPatchStrategy with
    *                  another dimension is registered, it will be ignored.
    */
   static hier::IntVector
   getMaxCoarsenOpStencilWidth(
      const tbox::Dimension& dim);

   /*!
    * @brief Constructor.
    *
    * The constructor will register the constructed object with a static
    * set that manages all CoarsenPatchStrategy objects in an application.
    */
   explicit CoarsenPatchStrategy(
      const tbox::Dimension& dim);

   /*!
    * @brief Destructor.
    */
   virtual ~CoarsenPatchStrategy();

   /*!
    * @brief Return maximum stencil width needed for user-defined
    * data coarsening operations performed by this object.
    *
    * This is needed to determine the correct coarsening data dependencies and
    * to ensure that the data has a sufficient amount of ghost width.
    *
    * For any user-defined coarsening operations implemented in the
    * preprocess or postprocess methods, return the maximum stencil needed
    * on a fine patch to coarsen data to a coarse patch.
    */
   virtual hier::IntVector
   getCoarsenOpStencilWidth() const = 0;

   /*!
    * @brief Perform user-defined patch data coarsening operations.
    *
    * This member function is called before standard coarsening operations
    * (expressed using concrete subclasses of the CoarsenOperator base class).
    * The preprocess function should move data from the source components
    * on the fine patch into the source components on the coarse patch
    * in the specified coarse box region.  Recall that the source components
    * are specified in calls to the registerCoarsen() function in the
    * CoarsenAlgorithm class.
    *
    * @param coarse[out] Coarse patch that will receive coarsened data.
    * @param fine[in]    Fine patch containing source data.
    * @param coarse_box[in]  Box region on coarse patch into which data is
    *                        coarsened.
    * @param ratio[in]   Refinement ratio between coarse and fine patches.
    */
   virtual void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * @brief Perform user-defined patch data coarsening operations.
    *
    * This member function is called after standard coarsening operations
    * (expressed using concrete subclasses of the CoarsenOperator base class).
    * The postprocess function should move data from the source components on
    * the fine patch into the source components on the coarse patch in the
    * specified coarse box region.  Recall that the source components are
    * specified in calls to the registerCoarsen() function in the
    * CoarsenAlgorithm class.
    *
    * @param coarse[out]  Coarse patch that will receive coarsened data.
    * @param fine[in]     Fine patch containing source data.
    * @param coarse_box[in]  hier::Box region on coarse patch into which data
    *                        is coarsened.
    * @param ratio[in]    Refinement ratio between coarse and fine patches.
    */
   virtual void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

private:
   /*!
    * @brief Get the set of CoarsenPatchStrategy objects that have been
    * registered.
    */
   static std::set<CoarsenPatchStrategy *>&
   getCurrentObjects();

   /*!
    * @brief Dimension of the object.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Register the object with a set of all CoarsenPatchStrategy
    * objects used in an application.
    */
   void
   registerObject();

};

}
}
#endif
