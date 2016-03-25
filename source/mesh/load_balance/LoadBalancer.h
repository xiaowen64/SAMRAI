//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/load_balance/LoadBalancer.h $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Load balance routines for uniform and non-uniform workloads.
//
 
#ifndef included_mesh_LoadBalancer
#define included_mesh_LoadBalancer
 
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_mesh_LoadBalanceStrategy
#include "LoadBalanceStrategy.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif

namespace SAMRAI {
    namespace mesh {

/*!
 * @brief Class LoadBalancer<DIM> provides load balancing routines for 
 * AMR hierarchy levels based on either uniform or non-uniform workload
 * estimates.  This class is derived from the abstract base class
 * LoadBalanceStrategy<DIM>; thus, it is a concrete implementation of the
 * load balance Strategy pattern interface.
 *
 * Load balancing operations, whether based on uniform or non-uniform
 * workloads, can be specified for each level in the hierarchy individually
 * or for the entire hierarchy.  Basic load balance parameters can be set
 * from an input file, while more complex behavior can be set at run-time
 * via member functions, including dynamic reconfiguration of balance
 * operations.
 *
 * Required input keys: NONE
 *
 * Optional input keys, data types, and defaults:
 * 
 *    - \b    processor_layout
 *       Integer array (length = DIM) indicating the way in which the
 *       domain should be chopped when a level can be described as a single
 *       parallelepiped region (i.e., a box). If no input value is provided,
 *       or if the product of these entries does not equal the number of
 *       processors, then the processor layout computed will be computed
 *       from the dimensions of the domain box and the number of processors
 *       in use if necessary.
 *
 *       NOTE: The largest patch size constraint specified in the input for
 *       the GriddingAlgorithm<DIM> object takes precedence over the processor
 *       layout specification.  That is, if the processor layout indicates
 *       that the resulting level patches would be larger than the largest
 *       patch size, the layout will be ignored and boxes obeying the patch
 *       size constrint will result.
 *
 *    - \b    bin_pack_method
 *       String value indicating the type of bin packing to use to map patches
 *       to processors.  Currently, two options are supported: "GREEDY",
 *       "SPATIAL". The "GREEDY" process simply maps each patch (box) to
 *       the first processor (bin), in ascending tbox::MPI process number, whose
 *       difference between the average workload and its current workload is
 *       less than the workload of the patch in question. The "SPATIAL"
 *       process first constructs an ordering of the patches (boxes) by
 *       passing a Morton-type curve through the center of each box.  Then,
 *       it attempts to map the patches to processors by assigning patches
 *       that are near each other on the curve to the same processor.  If no
 *       input value is specified, a default value of "SPATIAL" is used.  The
 *       input value will be used for all levels and will remain so until
 *       reset via the setBinPackMethod() member function below.
 *
 *    - \b   max_workload_factor
 *       Double used during the box-splitting phase to determine which boxes
 *       to split.  Specifically, boxes will be chopped if their estimated
 *       workload is greater than max_workload_factor * A, where A is the
 *       average workload (i.e., A = (total work)/(nume processors)).  The
 *       default value for this parameter is 1.0. It can be set to any value
 *       greater than zero, either in the input file or via the
 *       setMaxWorkloadFactor() member function below.   The default or input
 *       value will remain in use for all levels until reset by the member
 *       function.
 *
 *       NOTE: Setting this value > 1.0 increases the splitting threshold
 *       which effectively reduces the total number of boxes generated.
 *       Setting it less than 1.0 decreases the splitting threshold and will
 *       generally increase the total number of boxes.
 *
 *    - \b ignore_level_box_union_is_single_box
 *       Boolean flag to control chopping of level boxes when the union of 
 *       the input boxes passed to the loadBalanceBoxes() routine is a single box.
 *       The default value is false, which means that the domain will be chopped 
 *       to make patch boxes based on the (single box) union of the boxes describing
 *       the level regardless of the input boxes.  When the value is set to true,
 *       either via the setIgnoreLevelDomainIsSingleBox() function or an input file,
 *       the domain will be chopped by chopping each of the input boxes. 
 *
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 *
 *    processor_layout = 4 , 4 , 4    // number of processors is 64
 *    bin_pack = "GREEDY"
 *    max_workload_factor = 0.9
 *    ignore_level_box_union_is_single_box = TRUE
 *
 * \endverbatim 
 *
 * @see mesh::LoadBalanceStrategy
 */

template<int DIM> class LoadBalancer 
:
public LoadBalanceStrategy<DIM>
{
public:
   /*!
    * Construct load balancer object, including setting default object state 
    * and reading input data from the input data base, if required.
    *
    * @param name       User-defined string identifier used for error
    *                   reporting.  This string must be non-empty.
    * @param input_db   (optional) database pointer providing parameters from
    *                   input file.  This pointer may be null indicating no
    *                   input will be read.
    */
   LoadBalancer(const std::string& name,
                      tbox::Pointer<tbox::Database> input_db = (tbox::Database*)NULL);

   /*!
    * Construct load balancer object, including setting default object state
    * and reading input data from the input data base, if required.  The only
    * difference between this constructor and the previous one is the string 
    * identifier input.  If this constructor is used, the default object name
    * "LoadBalancer" applies.
    *
    * @param input_db   (optional) database pointer providing parameters from
    *                   input file.  This pointer may be null indicating no
    *                   input will be read.
    */
   LoadBalancer(tbox::Pointer<tbox::Database> input_db = (tbox::Database*)NULL);

   /*!
    * The virtual destructor releases all internal storage.
    */
   virtual ~LoadBalancer<DIM>();

   /*!
    * Set the max workload factor for either the specified level or all
    * hierarchy levels.  See discussion about inputs above for information
    * on how this value is used during load balancing operations.
    *
    * @param factor        Double value of multiplier for average workload
    *                      used in box chopping.  The default value is 1.0.
    * @param level_number  Optional integer number for level to which factor
    *                      is applied. If no value is given, the factor will
    *                      be used for all levels.
    */
   void setMaxWorkloadFactor(
      double factor,
      int level_number = -1);

   /*!
    * Configure the load balancer to use the data stored in the hierarchy at
    * the specified descriptor index for estimating the workload on each cell.
    *
    * @param data_id       Integer value of patch data identifier for workload
    *                      estimate on each cell.  An invalid value (i.e., < 0)
    *                      indicates that a spatially-uniform work estimate
    *                      will be used.  The default value is -1 (undefined)
    *                      implying the uniform work estimate.
    * @param level_number  Optional integer number for level on which data id
    *                      is used.  If no value is given, the data will be
    *                      used for all levels.
    */
   void setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1);

   /*!
    * Configure the load balancer to load balance boxes by assuming all cells
    * on the specified level or all hierarchy levels are weighted equally.
    *
    * @param level_number  Optional integer number for level on which uniform
    *                      workload estimate will be used.  If the level
    *                      number is not specified, a uniform workload
    *                      estimate will be used on all levels.
    */
   void setUniformWorkload(int level_number = -1);

   /*!
    * Configure the load balancer to use the bin-packing procedure for
    * mapping patches to processors indicated by the string.
    *
    * @param method        String value indicating bin-packing method to use.
    *                      See input file description above for valid options.
    *                      The default value is "GREEDY".
    * @param level_number  Optional integer number for level on which
    *                      bin-packing method will be used. If no value is
    *                      given, the prescribed methods will be used on all
    *                      levels.
    */
   void setBinPackMethod(
      const std::string& method,
      int level_number = -1);

   /*!
    * Set the boolean flag to control chopping of level boxes when the union of
    * the input boxes passed to the loadBalanceBoxes() routine is a single box.
    * The default value is false, which means that the domain will be chopped
    * to make patch boxes based on the (single box) union of the boxes describing
    * the level regardless of the input boxes.  When the value is set to true,
    * the domain will be chopped by chopping each of the input boxes.
    *
    * @param flag          Boolean value indicating whether to ignore the set of 
    *                      input boxes to the loadBalanceBoxes() routine when the
    *                      union of those boxes is a single box.
    */
   void setIgnoreLevelDomainIsSingleBox(bool flag);

   /*!
    * Return true if load balancing procedure for given level depends on
    * patch data on mesh; otherwise return false.  This can be used to
    * determine whether a level needs to be rebalanced although its box
    * configuration is unchanged.  This function is pure virtual in
    * the LoadBalanceStrategy<DIM> base class.
    *
    * @return Boolean value indicating whether load balance depends on 
    *         patch data.
    *
    * @param level_number  Integer patch level number.
    */
   bool getLoadBalanceDependsOnPatchData(int level_number) const;

   /*!
    * Given a list of boxes, representing the domain of a level in the AMR
    * hierarchy, generate an array of boxes and an associated processor
    * mapping from which the patches for the level will be generated and
    * assigned.  The resulting boxes and processor mapping will be determined
    * based on parameters set via input or member functions above.  This 
    * function is pure virtual in the LoadBalanceStrategy<DIM> base class.
    *
    * @param out_boxes       Output box array for generating patches on level.
    * @param mapping         Output processor mapping for patches on level.
    * @param in_boxes        Input box list representing union of patches on level.
    * @param hierarchy       Input patch hierarchy in which level will reside.
    * @param level_number    Input integer number of level in patch hierarchy.
    *                        This value must be >= 0.
    * @param physical_domain Array of boxes describing the physical extent of
    *                        the problem domain in the index space associated
    *                        with the level.  This box array cannot be empty.
    * @param ratio_to_hierarchy_level_zero  Input integer vector indicating
    *                        ratio between index space of level to load balance
    *                        and hierarchy level 0 (i.e., coarsest hierarchy level).
    * @param min_size        Input integer vector of minimum dimensions for
    *                        output boxes. All entries must be > 0.
    * @param max_size        Input integer vector of maximum dimensions for
    *                        output boxes. All entries must be >= min_size.
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
    */
   void loadBalanceBoxes(
      hier::BoxArray<DIM>& out_boxes,
      hier::ProcessorMapping& mapping,
      const hier::BoxList<DIM>& in_boxes,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      int level_number,
      const hier::BoxArray<DIM>& physical_domain,
      const hier::IntVector<DIM>& ratio_to_hierarchy_level_zero,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& max_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval) const;

   /*!
    * Print out all members of the class instance to given output stream.
    */
   virtual void printClassData(std::ostream& os) const; 

private:
   // The following are not implemented, but are provided here for
   // dumb compilers.
   LoadBalancer(const LoadBalancer<DIM>&);
   void operator=(const LoadBalancer<DIM>&);

   /*
    * Read parameters from input database.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   /*
    * Chop single box using uniform workload estimate.
    */
   void chopUniformSingleBox(
      hier::BoxArray<DIM>& out_boxes,
      tbox::Array<double>& out_workloads,
      const hier::Box<DIM>& in_box,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& max_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval,
      const hier::BoxArray<DIM>& physical_domain) const;

   /*
    * Chop boxes in list using uniform workload estimate.
    */
   void chopBoxesWithUniformWorkload(
      hier::BoxArray<DIM>& out_boxes,
      tbox::Array<double>& out_workloads,
      const hier::BoxList<DIM>& in_boxes,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      int level_number,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& max_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval,
      const hier::BoxArray<DIM>& physical_domain) const;

   /*
    * Chop boxes in list using non-uniform workload estimate.
    */
   void chopBoxesWithNonuniformWorkload(
      hier::BoxArray<DIM>& out_boxes,
      tbox::Array<double>& out_workloads,
      const hier::BoxList<DIM>& in_boxes,
      const tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
      int level_number,
      const hier::IntVector<DIM>& ratio_to_coarsest_hierarchy_level,
      int wrk_indx,
      const hier::IntVector<DIM>& min_size,
      const hier::IntVector<DIM>& max_size,
      const hier::IntVector<DIM>& cut_factor,
      const hier::IntVector<DIM>& bad_interval,
      const hier::BoxArray<DIM>& physical_domain) const;

   /*
    * Map boxes to processors using chosen bin pack method.
    */
   void binPackBoxes(
      hier::BoxArray<DIM>& boxes,
      hier::ProcessorMapping& mapping,
      tbox::Array<double>& workloads,
      const std::string& bin_pack_method) const;
 
   /*
    * Utility functions to determine parameter values for level.
    */
   int    getWorkloadDataId(int level_number) const;
   double getMaxWorkloadFactor(int level_number) const;
   std::string getBinPackMethod(int level_number) const;

   /*
    * String identifier for load balancer object.
    */
   std::string d_object_name;

   /*
    * Specification of processor layout.
    */
   bool d_processor_layout_specified;
   hier::IntVector<DIM> d_processor_layout;

   /*
    * Flag to control domain chopping when union of boxes for level 
    * is a single box.
    */
   bool d_ignore_level_box_union_is_single_box;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * that will be used for all levels unless specified for individual levels.
    */
   int    d_master_workload_data_id;
   double d_master_max_workload_factor;
   std::string d_master_bin_pack_method;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   tbox::Array<int>    d_workload_data_id;
   tbox::Array<double> d_max_workload_factor;
   tbox::Array<std::string> d_bin_pack_method;

   bool d_opt_for_single_box;

   /*
    * Performance timers.
    */
   tbox::Pointer<tbox::Timer> t_load_balance_boxes;
   tbox::Pointer<tbox::Timer> t_load_balance_boxes_remove_intersection;
   tbox::Pointer<tbox::Timer> t_bin_pack_boxes;
   tbox::Pointer<tbox::Timer> t_bin_pack_boxes_sort;
   tbox::Pointer<tbox::Timer> t_bin_pack_boxes_pack;
   tbox::Pointer<tbox::Timer> t_chop_boxes;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LoadBalancer.C"
#endif
