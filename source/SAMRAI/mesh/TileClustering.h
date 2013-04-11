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
      boost::shared_ptr<hier::BoxLevel>& new_box_level,
      boost::shared_ptr<hier::Connector>& tag_to_new,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::BoxContainer& bound_boxes,
      const hier::IntVector& min_box,
      const double efficiency_tol,
      const double combine_tol,
      const hier::IntVector& max_gcw);

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
    * Create, populate and return a coarsened version of the given tag data.
    *
    * The coarse cell values are set to tag_data if any corresponding
    * fine cell value is tag_data.  Otherwise, the coarse cell value
    * is set to zero.
    */
   boost::shared_ptr<pdat::CellData<int> >
   makeCoarsenedTagData(const pdat::CellData<int> &tag_data,
                        int tag_value) const;

   void setTimers();

   const tbox::Dimension d_dim;

   //! @brief Box size constraint.
   hier::IntVector d_box_size;

   /*!
    * @brief Whether to coalesce tiled-boxes after clustering to
    * create boxes bigger than tile size.
    */
   bool d_coalesce_boxes;

   //@{
   //! @name Diagnostics and performance evaluation
   hier::OverlapConnectorAlgorithm d_oca;
   bool d_log_cluster_summary;
   bool d_log_cluster;
   bool d_barrier_and_time;
   bool d_print_steps;
   boost::shared_ptr<tbox::Timer> t_find_boxes_containing_tags;
   boost::shared_ptr<tbox::Timer> t_cluster;
   boost::shared_ptr<tbox::Timer> t_global_reductions;
   //@}

};

}
}

#endif  // included_mesh_TileClustering
