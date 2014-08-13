/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontTagger class implementation
 *
 ************************************************************************/
#include "SinusoidalFrontTagger.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iomanip>

using namespace SAMRAI;

SinusoidalFrontTagger::SinusoidalFrontTagger(
   const std::string& object_name,
   const tbox::Dimension& dim,
   tbox::Database* database):
   d_name(object_name),
   d_dim(dim),
   d_amplitude(0.2),
   d_ghost_cell_width(dim, 0),
   d_buffer_cells(dim, 1),
   d_allocate_data(true),
   d_time(0.5)
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
   TBOX_ASSERT(variable_db != 0);

   std::vector<double> init_disp;
   std::vector<double> velocity;
   std::vector<double> period;

   if (database != 0) {
      d_allocate_data =
         database->getBoolWithDefault("allocate_data",
            d_allocate_data);
      if (database->isInteger("buffer_cells")) {
         database->getIntegerArray("buffer_cells",
            &d_buffer_cells[0], d_dim.getValue());
      }
      for (int ln = 0; true; ++ln) {
         std::string name("buffer_space_");
         name = name + tbox::Utilities::intToString(ln);
         if (database->isDouble(name)) {
            d_buffer_space.resize(d_dim.getValue() * (ln + 1));
            database->getDoubleArray(name, &d_buffer_space[d_dim.getValue() * ln], d_dim.getValue());
         } else {
            break;
         }
      }
      if (database->isDouble("period")) {
         period = database->getDoubleVector("period");
      }
      if (database->isDouble("init_disp")) {
         init_disp = database->getDoubleVector("init_disp");
      }
      if (database->isDouble("velocity")) {
         velocity = database->getDoubleVector("velocity");
      }
      d_amplitude =
         database->getDoubleWithDefault("amplitude",
            d_amplitude);
      d_time =
         database->getDoubleWithDefault("time",
            d_time);

      if (database->isInteger("ghost_cell_width")) {
         database->getIntegerArray("ghost_cell_width",
            &d_ghost_cell_width[0], d_dim.getValue());
      }
   }

   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      d_init_disp[idim] =
         idim < static_cast<int>(init_disp.size()) ? init_disp[idim] : 0.0;
      d_velocity[idim] =
         idim < static_cast<int>(velocity.size()) ? velocity[idim] : 0.0;
      d_period[idim] =
         idim < static_cast<int>(period.size()) ? period[idim] : 1.0e20;
   }

   const std::string context_name = d_name + std::string(":context");
   d_context = variable_db->getContext(context_name);

   boost::shared_ptr<hier::Variable> dist_var(
      new pdat::NodeVariable<double>(dim, d_name + ":dist"));
   d_dist_id = variable_db->registerVariableAndContext(dist_var,
         d_context,
         d_ghost_cell_width);

   boost::shared_ptr<hier::Variable> tag_var(
      new pdat::CellVariable<int>(dim, d_name + ":tag"));
   d_tag_id = variable_db->registerVariableAndContext(tag_var,
         d_context,
         d_ghost_cell_width);

   t_setup = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::setup");
   t_node_pos = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::node_pos");
   t_distance = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::distance");
   t_uval = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::uval");
   t_tag_cells = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::tag_cells");
   t_copy = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::copy");
}

SinusoidalFrontTagger::~SinusoidalFrontTagger()
{
}

void SinusoidalFrontTagger::initializeLevelData(
   /*! Hierarchy to initialize */
   const boost::shared_ptr<hier::PatchHierarchy>& base_hierarchy,
   /*! Level to initialize */
   const int ln,
   const double init_data_time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
   const boost::shared_ptr<hier::PatchLevel>& old_base_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(old_base_level);

   TBOX_ASSERT(base_hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   boost::shared_ptr<hier::PatchLevel> level(
      base_hierarchy->getPatchLevel(ln));

   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      hier::Patch& patch = **pi;
      initializePatchData(patch,
         init_data_time,
         initial_time,
         allocate_data);
   }

#if 0
   if (d_allocate_data) {
      /*
       * If instructed, allocate all patch data on the level.
       * Allocate only persistent data.  Scratch data will
       * generally be allocated and deallocated as needed.
       */
      if (allocate_data) {
         level->allocatePatchData(d_dist_id);
         level->allocatePatchData(d_tag_id);
      }
      computeLevelData(base_hierarchy, ln, d_time /*init_data_time*/,
         d_dist_id, d_tag_id, old_base_level);
   }
#endif
}

void SinusoidalFrontTagger::initializePatchData(
   hier::Patch& patch,
   const double init_data_time,
   const bool initial_time,
   const bool allocate_data)
{
   NULL_USE(initial_time);

   if (d_allocate_data) {
      /*
       * If instructed, allocate all patch data on the level.
       * Allocate only persistent data.  Scratch data will
       * generally be allocated and deallocated as needed.
       */
      if (allocate_data) {
         if (!patch.checkAllocated(d_dist_id)) {
            patch.allocatePatchData(d_dist_id);
         }
         if (!patch.checkAllocated(d_tag_id)) {
            patch.allocatePatchData(d_tag_id);
         }
         boost::shared_ptr<pdat::NodeData<double> > dist_data(
            BOOST_CAST<pdat::NodeData<double>, hier::PatchData>(
               patch.getPatchData(d_dist_id)));
         boost::shared_ptr<pdat::CellData<int> > tag_data(
            BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
               patch.getPatchData(d_tag_id)));
         TBOX_ASSERT(dist_data);
         TBOX_ASSERT(tag_data);
         computePatchData(patch, init_data_time,
            dist_data.get(), 0, tag_data.get());
      }
   }
}

void SinusoidalFrontTagger::resetHierarchyConfiguration(
   /*! New hierarchy */ const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
   /*! Coarsest level */ int coarsest_level,
   /*! Finest level */ int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
   d_hierarchy = new_hierarchy;
   TBOX_ASSERT(d_hierarchy);
}

void SinusoidalFrontTagger::applyGradientDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& base_hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation);
   TBOX_ASSERT(base_hierarchy_);
   boost::shared_ptr<hier::PatchLevel> level_(
      base_hierarchy_->getPatchLevel(ln));
   TBOX_ASSERT(level_);

   hier::PatchLevel& level = *level_;

   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi) {
      hier::Patch& patch = **pi;

      boost::shared_ptr<pdat::CellData<int> > tag_cell_data_(
         BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
      TBOX_ASSERT(tag_cell_data_);

      if (d_allocate_data) {
         // Use internally stored data.
         boost::shared_ptr<hier::PatchData> saved_tag_data =
            patch.getPatchData(d_tag_id);
         TBOX_ASSERT(saved_tag_data);
         tag_cell_data_->copy(*saved_tag_data);
      } else {
         // Compute tag data for patch.
         computePatchData(patch,
            error_data_time,
            0,
            0,
            tag_cell_data_.get());
      }

   }
}

/*
 * Deallocate patch data allocated by this class.
 */

void SinusoidalFrontTagger::deallocatePatchData(
   hier::PatchHierarchy& hierarchy)
{
   int ln;
   for (ln = 0; ln < hierarchy.getNumberOfLevels(); ++ln) {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy.getPatchLevel(ln));
      deallocatePatchData(*level);
   }
}

/*
 * Deallocate patch data allocated by this class.
 */

void SinusoidalFrontTagger::deallocatePatchData(
   hier::PatchLevel& level)
{
   level.deallocatePatchData(d_dist_id);
   level.deallocatePatchData(d_tag_id);
}

/*
 * Deallocate patch data allocated by this class.
 */
void SinusoidalFrontTagger::computeHierarchyData(
   hier::PatchHierarchy& hierarchy,
   double time)
{
   d_time = time;
   if (!d_allocate_data) return;

   for (int ln = 0; ln < hierarchy.getNumberOfLevels(); ++ln) {
      computeLevelData(hierarchy, ln, time, d_dist_id, d_tag_id);
   }
}

/*
 * Compute the solution data for a level.
 * Can copy data from old level (if any) to support
 * initializeLevelData().
 */

void SinusoidalFrontTagger::computeLevelData(
   const hier::PatchHierarchy& hierarchy,
   const int ln,
   const double time,
   const int dist_id,
   const int tag_id,
   const boost::shared_ptr<hier::PatchLevel>& old_level) const
{
   NULL_USE(old_level);

   const boost::shared_ptr<hier::PatchLevel> level(
      hierarchy.getPatchLevel(ln));

   /*
    * Initialize data in all patches in the level.
    */
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      hier::Patch& patch = **pi;
      boost::shared_ptr<pdat::NodeData<double> > dist_data;
      if (dist_id >= 0) {
         dist_data =
            boost::dynamic_pointer_cast<pdat::NodeData<double>, hier::PatchData>(
               patch.getPatchData(dist_id));
      }
      boost::shared_ptr<pdat::CellData<int> > tag_data;
      if (tag_id >= 0) {
         tag_data =
            boost::dynamic_pointer_cast<pdat::CellData<int>, hier::PatchData>(
               patch.getPatchData(tag_id));
      }
      computePatchData(patch, time,
         dist_data.get(),
         0,
         tag_data.get());
   }
}

/*
 * Compute various data for a patch.
 */

void SinusoidalFrontTagger::computePatchData(
   const hier::Patch& patch,
   const double time,
   pdat::NodeData<double>* dist_data,
   pdat::CellData<double>* uval_data,
   pdat::CellData<int>* tag_data) const
{

   t_setup->start();

   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT(patch.inHierarchy());

   const int ln = patch.getPatchLevelNumber();
   const boost::shared_ptr<hier::PatchLevel> level(
      d_hierarchy->getPatchLevel(ln));
   const hier::IntVector& ratio(level->getRatioToLevelZero());

   boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));

   TBOX_ASSERT(patch_geom);

   const double* xlo = patch_geom->getXLower();

   const double* dx = patch_geom->getDx();

   /*
    * Compute the size of buffer to tag around cells crossing front.
    * They should be at least d_buffer_cells in reference index space
    * and at least getBufferSpace(ln) in physical space.
    */
   hier::IntVector buffer(d_buffer_cells);
   for (int i = 0; i < d_dim.getValue(); ++i) {
      const double* buffer_space = getBufferSpace(ln);
      if (buffer_space != 0) {
         int space_based_buffer =
            int(d_buffer_space[ln * d_dim.getValue() + i] / dx[i] + 0.5);
         if (space_based_buffer > buffer(i)) {
            buffer(i) = space_based_buffer;
         }
      }
   }
   buffer *= ratio;

   t_setup->stop();

   /*
    * Determine the max ghost data we need to compute.
    */
   hier::Box fill_box = patch.getBox();
   if (tag_data) fill_box += tag_data->getGhostBox();
   if (uval_data) fill_box += uval_data->getGhostBox();
   if (dist_data) fill_box += dist_data->getGhostBox();

   computeFrontsData(dist_data, uval_data, tag_data,
      fill_box, buffer, xlo, dx, time);
}

/*
 * Compute various data due to the fronts.
 */
void SinusoidalFrontTagger::computeFrontsData(
   pdat::NodeData<double>* dist_data,
   pdat::CellData<double>* uval_data,
   pdat::CellData<int>* tag_data,
   const hier::Box& fill_box,
   const hier::IntVector& tag_buffer,
   const double xlo[],
   const double dx[],
   const double time) const
{
   TBOX_ASSERT(dist_data != 0 || uval_data != 0 || tag_data != 0);
   if (dist_data != 0 && tag_data != 0) {
      TBOX_ASSERT(dist_data->getBox().isSpatiallyEqual(tag_data->getBox()));
   }
   if (dist_data != 0 && uval_data != 0) {
      TBOX_ASSERT(dist_data->getBox().isSpatiallyEqual(uval_data->getBox()));
   }
   if (tag_data != 0 && uval_data != 0) {
      TBOX_ASSERT(tag_data->getBox().isSpatiallyEqual(uval_data->getBox()));
   }

   hier::Box pbox(d_dim);
   if (tag_data != 0) {
      pbox = tag_data->getBox();
   } else if (uval_data != 0) {
      pbox = uval_data->getBox();
   } else {
      pbox = dist_data->getBox();
   }

   /*
    * Initialize node x-distances from front.
    */

   double wave_number[SAMRAI::MAX_DIM_VAL];
   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      wave_number[idim] = 2 * 3.141592654 / d_period[idim];
   }

   /*
    * Compute box to contain the nodes' front position.  The box
    * should be big enough for computing data in fill_box.
    * Squash the front box to a single i-plane because front position
    * is indepenent of i.
    */
   hier::Box front_box = fill_box;
   front_box.grow(tag_buffer);
   front_box.growUpper(hier::IntVector(d_dim, 1));
   front_box.setUpper(0, pbox.lower(0));
   front_box.setLower(0, pbox.lower(0));

   const int ifront = front_box.lower(0);

   /*
    * Compute front_x: the front's x-coordinates, as a function of j-
    * and k- indices.
    */
   t_node_pos->start();

   pdat::ArrayData<double> front_x(front_box, 1);
   pdat::ArrayData<int>::iterator aiend(front_x.getBox(), false);
   for (pdat::ArrayData<int>::iterator ai(front_x.getBox(), true);
        ai != aiend; ++ai) {

      const hier::Index& index = *ai;
      double y = 0.0, siny = 0.0, z = 0.0, sinz = 1.0;

      y = xlo[1] + dx[1] * (index(1) - pbox.lower()[1]);
      siny = sin(wave_number[1] * (y + d_init_disp[1] - d_velocity[1] * time));

      if (d_dim.getValue() > 2) {
         z = xlo[2] + dx[2] * (index(2) - pbox.lower()[2]);
         sinz = sin(wave_number[2] * (z + d_init_disp[2] - d_velocity[2] * time));
      }

      front_x(index, 0) = d_velocity[0] * time + d_init_disp[0]
         + d_amplitude * siny * sinz;
   }
   t_node_pos->stop();

   if (tag_data != 0) {

      t_setup->start();

      /*
       * We need at least tag_buffer ghost cells outside fill_box to
       * compute the that come from buffering tags outside fill_box.
       * Create temporary patch data with the required tag_buffer for
       * computing tag values.  (We could give the real data the
       * required ghost cells, but that may affect the regridding
       * algorithm I'm testing.)
       */
      pdat::CellData<int> tmp_tag(fill_box, 1, tag_buffer);

      t_setup->stop();

      /*
       * Initialize tmp_tag to zero then tag specific cells.
       */
      hier::Box tag_fill_box = tmp_tag.getGhostBox() * fill_box;
      tmp_tag.fill(0, tag_fill_box);
      hier::BlockId blk0(0);
      pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(tag_fill_box));
      for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(tag_fill_box));
           ci != ciend; ++ci) {

         const pdat::CellIndex& cell_index = *ci;
         const hier::Box cell_box(cell_index, cell_index, blk0);

         /*
          * Compute the cell's min and max (signed) distance from
          * front, which are the min/max of its corners from the
          * front.  If the min and max have different signs, then the
          * cell overlaps the front surface and are tagged.
          */
         double min_distance_to_front = tbox::MathUtilities<double>::getMax();
         double max_distance_to_front = -tbox::MathUtilities<double>::getMax();
         pdat::NodeIterator niend(pdat::NodeGeometry::end(cell_box));
         for (pdat::NodeIterator ni(pdat::NodeGeometry::begin(cell_box));
              ni != niend; ++ni) {

            const pdat::NodeIndex& node_index = *ni;
            hier::Index front_index = node_index;
            front_index(0) = pbox.lower(0);

            double node_x = xlo[0] + dx[0] * (node_index(0) - pbox.lower() (0));

            double distance_to_front = node_x - front_x(front_index, 0);
            min_distance_to_front = tbox::MathUtilities<double>::Min(min_distance_to_front,
                  distance_to_front);
            max_distance_to_front = tbox::MathUtilities<double>::Max(max_distance_to_front,
                  distance_to_front);

         }
         while (min_distance_to_front < -0.5 * d_period[0]) {
            min_distance_to_front += d_period[0];
            max_distance_to_front += d_period[0];
         }
         while (max_distance_to_front > 0.5 * d_period[0]) {
            min_distance_to_front -= d_period[0];
            max_distance_to_front -= d_period[0];
         }
         if (min_distance_to_front < 0 && max_distance_to_front > 0) {
            // This cell has nodes on both sides of the front.  Tag it and the tag_buffer around it.
            hier::Box cell_and_buffer(cell_index, cell_index, blk0);
            cell_and_buffer.grow(tag_buffer);
            tmp_tag.fill(1, cell_and_buffer);
         }

      }

      t_copy->start();
      tag_data->copy(tmp_tag);
      t_copy->stop();

   }

   /*
    * Initialize U-value data.
    * The exact value of U increases by 1 across each front.
    */
   if (uval_data != 0) {
      t_uval->start();

      pdat::CellData<double>& uval(*uval_data);
      hier::Box uval_fill_box = uval.getGhostBox() * fill_box;
      uval.fill(0.0, uval_fill_box);
      const pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(uval_fill_box));
      for (pdat::CellData<double>::iterator ci = pdat::CellGeometry::begin(uval_fill_box);
           ci != ciend; ++ci) {
         const pdat::CellIndex& cindex = *ci;
         pdat::CellIndex squashed_cindex = cindex;
         squashed_cindex(0) = front_box.lower(0);
         double cellx = xlo[0] + dx[0] * (cindex(0) - pbox.lower(0) + 0.5);
         // Approximate cell's distance to front as average of its node distances.
         double dist_from_front = 0.0;
         if (d_dim == tbox::Dimension(2)) {
            dist_from_front = cellx - 0.5 * (
                  front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LowerLeft), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::UpperLeft), 0));
         } else if (d_dim == tbox::Dimension(3)) {
            dist_from_front = cellx - 0.25 * (
                  front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LLL), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LUL), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LLU), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LUU), 0));
         }
         while (dist_from_front > 0) {
            dist_from_front -= d_period[0];
            uval(cindex) += 1.0;
         }
      }

      t_uval->stop();
   }

   /*
    * Initialize distance data.
    */
   if (dist_data != 0) {
      t_distance->start();

      pdat::NodeData<double>& dist(*dist_data);
      hier::Box dist_fill_box = dist.getGhostBox() * fill_box;
      pdat::NodeData<double>::iterator ni(pdat::NodeGeometry::begin(dist_fill_box));
      pdat::NodeData<double>::iterator niend(pdat::NodeGeometry::end(dist_fill_box));
      for ( ; ni != niend; ++ni) {
         const pdat::NodeIndex& index = *ni;
         pdat::NodeIndex front_index(index);
         front_index(0) = ifront;
         dist(index) = xlo[0] + (index(0) - pbox.lower(0)) * dx[0]
            - front_x(front_index, 0);
      }

      t_distance->stop();
   }

}

#ifdef HAVE_HDF5
int SinusoidalFrontTagger::registerVariablesWithPlotter(
   appu::VisItDataWriter& writer)
{
   /*
    * Register variables with plotter.
    */
   if (d_allocate_data) {
      writer.registerPlotQuantity("Distance to front", "SCALAR", d_dist_id);
      writer.registerPlotQuantity("Tag value", "SCALAR", d_tag_id);
   } else {
      writer.registerDerivedPlotQuantity("Distance to front", "SCALAR", this,
         // hier::IntVector(0),
         1.0,
         "NODE");
      writer.registerDerivedPlotQuantity("Tag value", "SCALAR", this);
   }
   return 0;
}
#endif

bool SinusoidalFrontTagger::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index) const
{
   NULL_USE(region);
   NULL_USE(depth_index);

   TBOX_ASSERT(d_allocate_data == false);
   if (variable_name == "Distance to front") {
      pdat::NodeData<double> dist_data(patch.getBox(), 1, hier::IntVector(d_dim,
                                          0));
      computePatchData(patch, d_time, &dist_data, 0, 0);
      pdat::NodeData<double>::iterator ciend(pdat::NodeGeometry::end(patch.getBox()));
      for (pdat::NodeData<double>::iterator ci(pdat::NodeGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = dist_data(*ci);
      }
   } else if (variable_name == "Tag value") {
      pdat::CellData<int> tag_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      computePatchData(patch, d_time, 0, 0, &tag_data);
      pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(patch.getBox()));
      for (pdat::CellData<double>::iterator ci(pdat::CellGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = tag_data(*ci);
      }
   } else {
      TBOX_ERROR("Unrecognized name " << variable_name);
   }
   return true;
}

void SinusoidalFrontTagger::setTime(
   double time)
{
   d_time = time;
}

const double *SinusoidalFrontTagger::getBufferSpace(int ln) const
{
   if (static_cast<int>(d_buffer_space.size()) > ln * d_dim.getValue()) {
      return &d_buffer_space[ln * d_dim.getValue()];
   }
   return 0;
}
