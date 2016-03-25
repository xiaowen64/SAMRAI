//
// File:        BalanceUtilities.h
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 601 $
// Modified:    $Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
// Description: utility routines useful for load balancing operations
//
 
#ifndef included_mesh_BalanceUtilities
#define included_mesh_BalanceUtilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_ProcessorMapping
#include "ProcessorMapping.h"
#endif
#ifndef included_math_PatchCellDataNormOpsReal
#include "PatchCellDataNormOpsReal.h"
#endif
#ifndef included_mesh_SpatialKey
#include "SpatialKey.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace mesh {

/*!
 * @brief Utility class BalanceUtilities<DIM> provides several functions 
 * useful in various load balancing operations.  These utilities include
 * bin packing operations, box chopping by recursive bisection, and
 * computation of effective processor layouts for boxes. 
 */

template<int DIM> struct BalanceUtilities
{
   /*!
    * Assign workloads to processors using a greedy algorithm that attempts
    * to distribute the sum of weights on each processor evenly across 
    * the given number of processors.
    *
    * @return         double-valued estimate of the load balance efficiency 
    *                 (ranges from zero to one hundred percent)
    *  
    * @param mapping  Output processor mapping.
    * @param weights  tbox::Array of double-valued weights to distribute.
    * @param nproc    Integer number of processors, must be > 0.
    */
   static double binPack(
      hier::ProcessorMapping& mapping,
      tbox::Array<double>& weights,
      int nproc);

   /*!
    * Assign boxes to processors so that boxes spatially near each 
    * other are likely to be assigned to processors near each other 
    * (assuming that processor ordering is reflected in processor rank) 
    * and so that the workload is approximately evenly distributed among 
    * the processors.  The routine uses a Morton space-filling curve
    * algorithm.
    *
    * Note that this routine potentially reorders the boxes passed in to
    * achieve the first goal.
    *
    * @return         Double-valued estimate of the load balance efficiency
    *                 (ranges from zero to one hundred percent)
    *
    * @param mapping  Output processor mapping.
    * @param weights  tbox::Array of double-valued box weights to distribute.
    * @param boxes    tbox::Array of boxes to distribute to processors.
    * @param nproc    Integer number of processors, must be > 0.
    * 
    * Note that the wight and box arrrays must be the same size.
    */
   static double spatialBinPack(
      hier::ProcessorMapping& mapping,
      tbox::Array<double>& weights,
      hier::BoxArray<DIM>& boxes,
      const int nproc);

   /*!
    * Recursively chop chops boxes in input boxlist until each box has
    * a workload less than the prescribed ideal workload or no more more 
    * chopping is allowed by the given constraints.   A spatially-uniform 
    * workload is assumed; i.e., all cells are weighted equally.  This routine
    * attempts to create as many boxes as possible with loads equal to or 
    * slightly less than the ideal workload value so that they can be 
    * mapped to processors effectively.  
    *
    * @param out_boxes       Output box list.
    * @param out_workloads   Output list of box workloads.
    * @param in_boxes        Input boxlist for chopping.
    * @param ideal_workload  Input double ideal box workload, must be > 0.
    * @param min_size        Input integer vector of minimum dimensions for
    *                        output boxes. All entries must be > 0.
    * @param cut_factor      Input integer vector used to create boxes with
    *                        correct dimensions.  The length of each box
    *                        dimension will be an integer multiple of the
    *                        corresponding cut factor vector entry.  All
    *                        vector entries must be > 0.  See hier::BoxUtilities
    *                        documentation for more details.
    * @param bad_interval    Input integer vector used to create boxes near
    *                        physical domain boundary with sufficient number
    *                        of cells.  No box face will be closer to the
    *                        boundary than the corresponding interval of cells
    *                        to the boundary (the corresponding value is given
    *                        by the normal direction of the box face) unless
    *                        the face coincides with the boundary itself.  The
    *                        point of this argument is to have no patch live
    *                        within a certain ghost cell width of the boundary
    *                        if its boundary does not coincide with that
    *                        boundary .  That is, all ghost cells along a face
    *                        will be either in the domain interior or outside
    *                        the domain.  All entries must be >= 0. See
    *                        hier::BoxUtilities documentation for more details.
    * @param physical_domain tbox::Array of boxes describing the physical extent of
    *                        the index space associated with the in_boxes.
    *                        This box array cannot be empty.
    */
   static void recursiveBisectionUniform(
      hier::BoxList<DIM>& out_boxes,
      tbox::List<double>& out_workloads,
      const hier::BoxList<DIM>& in_boxes,
      double ideal_workload,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval,
      const hier::BoxArray<DIM>& physical_domain);

   /*!
    * Recursively chops boxes given by patches on input patch level until each 
    * box has a workload less than the prescribed ideal workload or no more more
    * chopping is allowed by the given constraints.   A spatially-nonuniform
    * workload is assumed.  Cell weights must be given bydata defined by the 
    * given patch data id on the given patch level.  This routine attempts to 
    * create as many boxes as possible with loads equal to or slightly less 
    * than the ideal workload value so that they can be mapped to processors 
    * effectively.
    *
    * @param out_boxes       Output box list.
    * @param out_workloads   Output list of box workloads.
    * @param in_level        Input patch level whose patches describe input
    *                        box regions and whose patch data contain workload 
    *                        estimate for each cell.
    * @param work_id         Input integer patch data id for cell-centered 
    *                        double work estimate for each cell.
    * @param ideal_workload  Input double ideal box workload, must be > 0.
    * @param min_size        Input integer vector of minimum dimensions for
    *                        output boxes. All entries must be > 0.
    * @param cut_factor      Input integer vector used to create boxes with
    *                        correct dimensions.  The length of each box
    *                        dimension will be an integer multiple of the
    *                        corresponding cut factor vector entry.  All
    *                        vector entries must be > 0.  See hier::BoxUtilities
    *                        documentation for more details.
    * @param bad_interval    Input integer vector used to create boxes near
    *                        physical domain boundary with sufficient number
    *                        of cells.  No box face will be closer to the
    *                        boundary than the corresponding interval of cells
    *                        to the boundary (the corresponding value is given
    *                        by the normal direction of the box face) unless
    *                        the face coincides with the boundary itself.  The
    *                        point of this argument is to have no patch live
    *                        within a certain ghost cell width of the boundary
    *                        if its boundary does not coincide with that
    *                        boundary .  That is, all ghost cells along a face
    *                        will be either in the domain interior or outside
    *                        the domain.  All entries must be >= 0. See
    *                        hier::BoxUtilities documentation for more details.
    * @param physical_domain tbox::Array of boxes describing the physical extent of
    *                        the index space associated with the in_boxes.
    *                        This box array cannot be empty.
    */
   static void recursiveBisectionNonuniform(
      hier::BoxList<DIM>& out_boxes,
      tbox::List<double>& out_workloads,
      const tbox::Pointer< hier::PatchLevel<DIM> >& in_level,
      int work_id,
      double ideal_workload,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval,
      const hier::BoxArray<DIM>& physical_domain);

   /*!
    * Compute factorization of processors corresponding to 
    * dimensions of given box.
    *
    * @param proc_dist  Output number of processors for each 
    *                   coordinate direction.
    * @param num_procs  Input integer number of processors, must be > 0.
    * @param box        Input box to be distributed.
    */
   static void computeDomainDependentProcessorLayout(
      hier::IntVector<DIM>& proc_dist,
      int num_procs,
      const hier::Box<DIM>& box);

   /*!
    * Compute a factorization of processors that does NOT necessarily
    * correspond to the dimensions of the supplied box.  For example, the 
    * processor distribution in each direction may simply be a square root 
    * (cube root in 3D) of the number of processors.  The box information
    * is used simply to determine a maximum number of processors in each
    * coordinate direction. 
    * 
    * @param proc_dist  Output number of processors for each
    *                   coordinate direction.
    * @param num_procs  Input integer number of processors, must be > 0.
    * @param box        Input box to be distributed.
    */
   static void computeDomainIndependentProcessorLayout(
      hier::IntVector<DIM>& proc_dist,
      int num_procs,
      const hier::Box<DIM>& box);

   /*!
    * Sort box array in descending order of workload according to the
    * workload array.  Both the box array and the work array will be 
    * sorted on return.  
    *
    * Note that if you simply want to sort boxes based on their size,
    * see hier::BoxUtilities. 
    * 
    * @param boxes     Boxes to be sorted based on workload array. 
    * @param workload  Workloads to use for sorting boxes.
    * 
    * Note that both arrays must be the same size.
    */
   static void sortDescendingBoxWorkloads(
      hier::BoxArray<DIM>& boxes,
      tbox::Array<double>& workload);

   /*!
    * Compute total workload in region of argument box based on patch
    * data defined by given integer index.  The sum is computed on the
    * intersection of argument box and box over which data associated with 
    * workload is defined. 
    *
    * @return          Double-valued sum of workload values in box region.
    *
    * @param patch     Input patch on which workload data is defined.
    * @param wrk_indx  Input integer patch data identifier for work data.
    * @param box       Input box region 
    *
    * Note that wrk_indx must refer to a valid cell-centered patch data entry.
    */
   static double computeNonUniformWorkload(
      tbox::Pointer< hier::Patch<DIM> > patch,
      int wrk_indx,
      const hier::Box<DIM>& box);

   /*!
    * Compute and return load balance efficiency for a level.  
    * 
    * @return         Double-valued estimate of the load balance efficiency
    *                 (ranges from zero to one hundred percent)
    *
    * @param level            Input patch level to consider, can't be null.
    * @param os               Output stream for reporting load balance 
    *                         details.
    * @param workload_data_id (Optional) Input integer id for workload 
    *                         data on level.  If no value is given, the 
    *                         calculation assumes spatially-uniform load.
    */
   static double computeLoadBalanceEfficiency(
      const tbox::Pointer< hier::PatchLevel<DIM> >& level,
      ostream& os,
      int workload_data_id = -1);

private:
   static math::PatchCellDataNormOpsReal<DIM,double> s_norm_ops;

   static void privateHeapify(
      tbox::Array<int>& permutation,
      tbox::Array<double>& workload,
      const int index,
      const int heap_size);
   
   static void privateHeapify(
      tbox::Array<int>& permutation,
      tbox::Array<SpatialKey>& spatial_keys,
      const int index,
      const int heap_size);

   static void privateRecursiveProcAssign(
      const int wt_index_lo,
      const int wt_index_hi,
      tbox::Array<double>& weights,
      const int proc_index_lo,
      const int proc_index_hi,
      hier::ProcessorMapping& mapping,
      const double avg_weight);

   static void privatePrimeFactorization(
      const int N, 
      tbox::Array<int>& p);

   static void privateResetPrimesArray(
      tbox::Array<int>& p);

   static bool privateBadCutPointsExist(
      const hier::BoxArray<DIM>& physical_domain);
   

   static void privateInitializeBadCutPointsForBox(
      tbox::Array< tbox::Array<bool> >& bad_cut_points,
      hier::Box<DIM>& box,
      bool bad_domain_boundaries_exist,
      const hier::IntVector<DIM>& bad_interval,
      const hier::BoxArray<DIM>& physical_domain);

   static bool privateFindBestCutDimension(
      int& cut_dim_out,
      const hier::Box<DIM>& in_box,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& cut_factor,
      tbox::Array< tbox::Array<bool> >& bad_cut_points);


   static int privateFindCutPoint(
      double total_work,
      double ideal_workload,
      int mincut,
      int numcells,
      const tbox::Array<double>& work_in_slice,
      const tbox::Array<bool>& bad_cut_points);
   
   static void privateCutBoxesAndSetBadCutPoints(
      hier::Box<DIM>& box_lo,
      tbox::Array< tbox::Array<bool> >& bad_cut_points_for_boxlo,
      hier::Box<DIM>& box_hi,
      tbox::Array< tbox::Array<bool> >& bad_cut_points_for_boxhi,
      const hier::Box<DIM>& in_box,
      int cutdim,
      int cut_index,
      const tbox::Array< tbox::Array<bool> >& bad_cut_points);


   static void privateRecursiveBisectionUniformSingleBox(
      hier::BoxList<DIM>& out_boxes,
      tbox::List<double>& out_workloads,
      const hier::Box<DIM>& in_box,
      double in_box_workload,
      double ideal_workload,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& cut_factor,
      tbox::Array< tbox::Array<bool> >& bad_cut_points);

   static void privateRecursiveBisectionNonuniformSingleBox(
      hier::BoxList<DIM>& out_boxes,
      tbox::List<double>& out_workloads,
      const tbox::Pointer< hier::Patch<DIM> >& patch,
      const hier::Box<DIM>& in_box, 
      double in_box_workload,
      int work_data_index,
      double ideal_workload,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& cut_factor,
      tbox::Array< tbox::Array<bool> >& bad_cut_points);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BalanceUtilities.C"
#endif
