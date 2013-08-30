/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Asynchronous Berger-Rigoutsos clustering algorithm.
 *
 ************************************************************************/
#ifndef included_mesh_TileClustering
#define included_mesh_TileClustering

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/mesh/BoxGeneratorStrategy.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"
#include "SAMRAI/tbox/Database.h"

#include "boost/shared_ptr.hpp"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Tiled patch clustering algorithm.
 * This is UNSUPPORTED, EXPERIMENTAL code not for general use.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *
 *   - \b box_size
 *   Box size in the index space of the tag level.
 *
 *   - \b coalesce_boxes
 *   Whether to coalesce boxes after clustering.  This can lead to
 *   clusters that are bigger than specified tile size.
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
 *     <td>box_size</td>
 *     <td>int[]</td>
 *     <td>all values are 8</td>
 *     <td>????????</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>coalesce_boxes</td>
 *     <td>bool</td>
 *     <td>true</td>
 *     <td>false/true</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @internal The following are developer inputs.
 * Defaults are listed in parenthesis:
 *
 * @internal DEV_print_steps (FALSE)
 * boolean
 */
class TileClustering:public BoxGeneratorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   explicit TileClustering(
      const tbox::Dimension& dim,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>());

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   virtual ~TileClustering();

   /*!
    * @brief Implement the mesh::BoxGeneratorStrategy interface
    * method of the same name.
    *
    * Create a set of boxes that covers all integer tags on
    * the patch level that match the specified tag value.
    * Each box will be at least as large as the given minimum
    * size and the tolerances will be met.
    */
   void
   findBoxesContainingTags(
      boost::shared_ptr<hier::BoxLevel>& new_box_level,
      boost::shared_ptr<hier::Connector>& tag_to_new,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::BoxContainer& bound_boxes,
      const hier::IntVector& min_box,
      const hier::IntVector& max_gcw);

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named
    * "mesh::BergerRigoutsosNode::*", where the third field is
    * the specific steps performed by the BergerRigoutsosNode.
    * You can override the first two fields with this method.
    * Conforming to the timer naming convention, timer_prefix should
    * have the form "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

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


   /*!
    * @brief Cluster, cutting off tiles at process boundaries.
    */
   void clusterWithinProcessBoundaries(
      hier::BoxLevel &new_box_level,
      hier::Connector &tag_to_new,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const hier::BoxContainer& bound_boxes,
      int tag_data_index,
      int tag_val);

   void
   clusterWithinLevelBoundaries(
      hier::BoxLevel &new_box_level,
      hier::Connector &tag_to_new,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const hier::BoxContainer& bound_boxes,
      int tag_data_index,
      int tag_val);

   /*!
    * @brief Create, populate and return a coarsened version of the
    * given tag data.
    *
    * The coarse cell values are set to tag_data if any corresponding
    * fine cell value is tag_data.  Otherwise, the coarse cell value
    * is set to zero.
    */
   boost::shared_ptr<pdat::CellData<int> >
   makeCoarsenedTagData(const pdat::CellData<int> &tag_data,
                        int tag_value) const;

   /*!
    * @brief Find tagged tiles in a single patch.
    */
   int
   findTilesContainingTags(
      hier::BoxContainer &tiles,
      const pdat::CellData<int> &tag_data,
      int tag_val,
      int first_tile_index);

   void
   clusterWholeTiles(
      hier::BoxLevel &new_box_level,
      boost::shared_ptr<hier::Connector> &tag_to_new,
      int &local_tiles_have_remote_extent,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const hier::BoxContainer& bound_boxes,
      int tag_data_index,
      int tag_val);

   void
   detectSemilocalEdges( boost::shared_ptr<hier::Connector> &tag_to_tile );

   void
   removeDuplicateTiles(
      hier::BoxLevel &tile_box_level,
      hier::Connector &tag_to_tiles);

   /*!
    * @brief Coalesce clusters (and update Connectors).
    */
   void
   coalesceClusters(
      hier::BoxLevel &tile_box_level,
      boost::shared_ptr<hier::Connector> &tag_to_tile,
      int tiles_have_remote_extent);

   void setTimers();

   const tbox::Dimension d_dim;

   //! @brief Box size constraint.
   hier::IntVector d_box_size;

   /*!
    * @brief Whether to allow tiles to have remote extents.
    *
    * If false, tiles will be cut at process boundaries, resulting in
    * completely local tiles.  If true, allow tiles to cross process
    * boundaries where, resulting in less tile fragmentation.
    */
   bool d_allow_remote_tile_extent;

   /*!
    * @brief Whether to coalesce tiled-boxes after clustering to
    * create boxes bigger than tile size.
    */
   bool d_coalesce_boxes;

   /*!
    * @brief Thread locker for modifying clustering outputs with multi-threads.
    */
   TBOX_omp_lock_t l_outputs;

   /*!
    * @brief Thread locker for modifying intermediate data with multi-threads.
    */
   TBOX_omp_lock_t l_interm;

   //@{
   //! @name Diagnostics and performance evaluation
   hier::OverlapConnectorAlgorithm d_oca;
   bool d_log_cluster_summary;
   bool d_log_cluster;
   bool d_barrier_and_time;
   bool d_print_steps;
   //@}


   //@{
   //! @name Performance timer data for this class.

   /*
    * @brief Structure of timers used by this class.
    *
    * Each object can set its own timer names through
    * setTimerPrefix().  This leads to many timer look-ups.  Because
    * it is expensive to look up timers, this class caches the timers
    * that has been looked up.  Each TimerStruct stores the timers
    * corresponding to a prefix.
    */
   struct TimerStruct {
      boost::shared_ptr<tbox::Timer> t_find_boxes_containing_tags;
      boost::shared_ptr<tbox::Timer> t_cluster;
      boost::shared_ptr<tbox::Timer> t_coalesce;
      boost::shared_ptr<tbox::Timer> t_coalesce_adjustment;
      boost::shared_ptr<tbox::Timer> t_global_reductions;
      boost::shared_ptr<tbox::Timer> t_cluster_setup;
      boost::shared_ptr<tbox::Timer> t_cluster_wrapup;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

   static int s_primary_mpi_tag;
   static int s_secondary_mpi_tag;
   static int s_first_data_length;

   /*!
    * @brief Structure of timers in s_static_timers, matching this
    * object's timer prefix.
    */
   TimerStruct* d_object_timers;

   //@}

};

}
}

#endif  // included_mesh_TileClustering
