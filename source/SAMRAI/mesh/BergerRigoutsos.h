/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Asynchronous Berger-Rigoutsos clustering algorithm.
 *
 ************************************************************************/
#ifndef included_mesh_BergerRigoutsos
#define included_mesh_BergerRigoutsos

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/mesh/BoxGeneratorStrategy.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/shared_ptr.hpp"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Asynchronous Berger-Rigoutsos implementation.
 * This class is derived from the abstract base class
 * mesh::BoxGeneratorStrategy.  Thus, it serves as a concrete
 * implementation of the box generator Strategy pattern interface.
 *
 * This class uses the BergerRigoutsosNode class
 * to carry out the asynchronous Berger-Rigoutsos algorithm.
 * It handles aspects not central to that algorithm.  It:
 * - Implements the box generator Strategy pattern interface.
 * - Provides an interface with the input database for setting
 *   parameters influencing the implementation.
 * - Sorts the output data (if user requests).
 * - Performs some additional error checking.
 * For more details on the parallel implementation,
 * see mesh::BergerRigoutsosNode.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b max_box_size
 *       The maximum cluster size allowed.  This parameter is not critical to
 *       clustering but limiting the cluster size may improve performance of
 *       load balancing algorithms (due to the excessive work required by the
 *       owner of huge clusters).
 *
 *    - \b sort_output_nodes
 *       Whether to sort the output.  This makes the normally non-deterministic
 *       ordering deterministic and the results repeatable.
 *
 *    - \b check_min_box_size
 *       A flag to control how to resolve an initial box that violates the
 *       minimum box size.  Set to one of these strings: <br>
 *       \b "IGNORE" - violations will be quietly disregarded. <br>
 *       \b "WARN" - violations will cause a warning but the code will
 *       continue anyway. <br>
 *       \b "ERROR" - violations will cause an unrecoverable assertion.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>max_box_size</td>
 *     <td>int[]</td>
 *     <td>all valules max int</td>
 *     <td>all value > 0</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>sort_output_nodes</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>check_min_box_size</td>
 *     <td>string</td>
 *     <td>"WARN"</td>
 *     <td>"WARN", "IGNORE", "ERROR"</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 * </table>
 */
class BergerRigoutsos:public BoxGeneratorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   explicit BergerRigoutsos(
      const tbox::Dimension& dim,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>());

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   virtual ~BergerRigoutsos();

   /*!
    * @brief Set the MPI communication object.
    *
    * Duplicate the communicator in the given for private use.  A private
    * communicator isolates the complex communications used by the
    * asynchronous algorithm from other communications.  Duplicating
    * the communicator is expensive but should only be need once.  All
    * processes in the communicator must participate.  The duplicate
    * communicator is active until this object is destructed.
    * Using a duplicated MPI communicator is optional but recommended.
    * When a duplicate MPI communicator is in use, it must be congruent
    * with the communicator associated with the tag level.
    *
    * If the communicator is not set, the parallel clustering
    * algorithm uses the communicator of the input tag
    * box_level.  If it is set, then the algorithm only works
    * for input tag box_levels with a congruent communicator.
    *
    * If communicator is SAMRAI_MPI::commNull, it is the same as not
    * using a duplicate communicator.
    */
   void
   setMPI(
      const tbox::SAMRAI_MPI& mpi);

   /*!
    * @brief Implement the mesh::BoxGeneratorStrategy interface
    * method of the same name.
    *
    * Create a set of boxes that covers all integer tags on
    * the patch level that match the specified tag value.
    * Each box will be at least as large as the given minimum
    * size and the tolerances will be met.
    *
    * The efficiency tolerance is a threshold value for the percentage of
    * tagged cells in each box.  If this percentage is below the tolerance,
    * the box will continue to be split into smaller boxes.
    *
    * The combine tolerance is a threshold value for the sum of the volumes
    * of two boxes into which a box may be potentially split.  If ratio of
    * that sum and the volume of the original box, the box will not be split.
    */
   void
   findBoxesContainingTags(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::Box& bound_box,
      const hier::IntVector& min_box,
      const double efficiency_tol,
      const double combine_tol,
      const hier::IntVector& max_gcw,
      const hier::BlockId& block_id,
      const hier::LocalId& first_local_id) const;

protected:
   /*!
    * @brief Read parameters from input database.
    *
    * @param input_db Input Database.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& input_db);

private:
   const tbox::Dimension d_dim;

   void
   assertNoMessageForPrivateCommunicator() const;

   void
   sortOutputBoxes(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag) const;

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   //! @brief Communication object.
   tbox::SAMRAI_MPI d_mpi;

   //! @brief Max box size constraint used by BergerRigoutsosNode.
   hier::IntVector d_max_box_size;

   //! @brief Max distance from center for Laplace cut.
   double d_max_lap_cut_from_center;

   //! @brief Threshold for avoiding thinner directions for Laplace cut.
   double d_laplace_cut_threshold_ar;

   //! @brief Whether to log execution node allocation and deallocation.
   bool d_log_node_history;

   //! @brief Whether to briefly log cluster summary.
   bool d_log_cluster_summary;

   //! @brief Whether to log cluster summary.
   bool d_log_cluster;

   //! @brief How to select the owner of a node.
   std::string d_owner_mode;

   //! @brief Asynchronous mode for advancing algorithm.
   std::string d_algo_advance_mode;

   //! @brief Whether to sort results to make them deterministic.
   bool d_sort_output_nodes;

   //! @brief How to resolve initial boxes smaller than min box size.
   char d_check_min_box_size;

   //@{
   //! @name Used for evaluating performance;
   bool d_barrier_before;
   bool d_barrier_after;
   //@}

   static boost::shared_ptr<tbox::Timer> t_barrier_before;
   static boost::shared_ptr<tbox::Timer> t_barrier_after;
   static boost::shared_ptr<tbox::Timer> t_find_boxes_with_tags;
   static boost::shared_ptr<tbox::Timer> t_run_abr;
   static boost::shared_ptr<tbox::Timer> t_global_reductions;
   static boost::shared_ptr<tbox::Timer> t_sort_output_nodes;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif  // included_mesh_BergerRigoutsos
