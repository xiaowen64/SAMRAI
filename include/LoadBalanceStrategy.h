//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/load_balance/LoadBalanceStrategy.h $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Strategy interface for box load balancing routines.
//
 
#ifndef included_mesh_LoadBalanceStrategy
#define included_mesh_LoadBalanceStrategy
 
#include "SAMRAI_config.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "ProcessorMapping.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace mesh {

/*!
 * @brief Class LoadBalanceStrategy<DIM> is an abstract base class that defines 
 * a Strategy pattern interface for operations that load balance patches
 * on a single AMR patch hierarchy level.  Typically, such operations are 
 * invoked after the domain of a new hierarchy level is determined (e.g., via
 * some error estimation procedure) and is applied to the collection of 
 * boxes that describe the domain.  The load balancing process produces
 * a set of boxes from which patches on the new level are created and a 
 * processor mapping describing how the new patches are mapped to processors.
 * 
 * @see hier::PatchLevel
 * @see hier::ProcessorMapping
 */

template<int DIM> class LoadBalanceStrategy  : public tbox::DescribedClass
{
public:
   /*!
    * Construct load balance strategy object.
    */
   LoadBalanceStrategy();

   /*!
    * This virtual destructor does nothing interesting.
    */
   virtual ~LoadBalanceStrategy<DIM>();

   /*!
    * Indicate whether load balancing procedure for given level depends on
    * patch data on mesh.  This can be used to determine whether a level 
    * needs to be rebalanced although its box configuration is unchanged. 
    * 
    * @return Boolean value of true if load balance routines for level
    *         depend on patch data; false otherwise.
    *
    * @param level_number Integer level number.
    */
   virtual bool getLoadBalanceDependsOnPatchData(int level_number) const = 0;

   /*!
    * Given a list of boxes, the union of which represents the domain of
    * a specified level in the AMR hierarchy, generate an array of boxes 
    * and an associated processor mapping from which the patches for the
    * level may be generated.  This process typically involves chopping each
    * box in the original list that is "too large" (in a manner defined by
    * the concrete implementation of this subclass) into a set of boxes 
    * each smaller than some size.  Thus, the union of the boxes in the 
    * generated box array is the same as that of the original box list. 
    * A variety of constraints must be applied typically during the chopping
    * process.  These input arguments are described here.
    *
    * @param out_boxes        Output box array for generating patches on level.
    * @param mapping          Output processor mapping for patches on level.
    * @param in_boxes         Input box list representing union of patches on level.
    * @param hierarchy        Input patch hierarchy in which level will reside.
    * @param level_number     Input integer number of level in patch hierarchy.
    * @param physical_domain  Input box array representing physical extent
    *                         of the problem domain in the index space of the
    *                         level to be load balanced.
    * @param ratio_to_hierarchy_level_zero  Input hier::IntVector indicating
    *                         ratio between index space of level and coarsest 
    *                         hierarchy level (i.e., level zero).
    * @param min_size         Input hier::IntVector representing mimimum box size.
    * @param max_size         Input hier::IntVector representing maximum box size.
    * @param cut_factor       Input hier::IntVector indicating factor for chopping
    *                         each side of a box; i.e., after chopping a box,
    *                         the number of cells along each direction of each 
    *                         piece must be an integer multiple of the corresponding
    *                         entry in the cut factor vector.  For example, the 
    *                         cut factor may be related to the coarsen ratio between 
    *                         levels in the hierarchy in which case it may be used
    *                         to produce boxes that can be coarsened by a certain 
    *                         factor if needed.  See hier::BoxUtilities<DIM> header file 
    *                         for more information.
    * @param bad_interval     Input hier::IntVector indicating the length of an interval 
    *                         of cells along each side of the box where chopping 
    *                         the box may produce boxes with certain "bad" properties. 
    *                         For example, this is primiarily used to avoid generating
    *                         ghost regions for patches that intersect the domain
    *                         boundary in ways that may it difficult for a use to 
    *                         provide boundary values.  Thus, it is typically related 
    *                         to the maximum ghost cell width in the problem.  See 
    *                         hier::BoxUtilities<DIM> header file for more information.
    */
   virtual void loadBalanceBoxes(
      hier::BoxArray<DIM>& out_boxes,
      hier::ProcessorMapping& mapping,
      const hier::BoxList<DIM>& in_boxes,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      const int level_number,
      const hier::BoxArray<DIM>& physical_domain,
      const hier::IntVector<DIM>& ratio_to_hierarchy_level_zero,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& max_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval) const = 0;

private:
   // The following are not implemented:
   LoadBalanceStrategy(const LoadBalanceStrategy<DIM>&);
   void operator=(const LoadBalanceStrategy<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LoadBalanceStrategy.C"
#endif
