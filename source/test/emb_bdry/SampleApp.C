/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Class to manage functions for QM calculations. 
 *
 ************************************************************************/

#include "SampleApp.h"

/*
 * Header files for SAMRAI library classes
 */
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"

/*
 * Fortran function(s) used, if fortran is desired.
 */
extern "C" {
#if (NDIM == 2)
void F77_FUNC(tagcells, TAGCELLS) (const int &, const int &, // ifirst0, ilast0
   const int &, const int &,               // ifirst1, ilast1
   const int &, const int &,               // fgh0, fgh1
   const int &, const int &,               // tgh0, tgh1
   const int &,                            // flag value to tag
   const int &,                            // tag_value
   const int *,                            // flag
   int *);                                 // tags
#endif
#if (NDIM == 3)
void F77_FUNC(tagcells, TAGCELLS) (const int &, const int &, // ifirst0, ilast0
   const int &, const int &,               // ifirst1, ilast1
   const int &, const int &,               // ifirst2, ilast2
   const int &, const int &, const int &,             // fgh0, fgh1, fgh2
   const int &, const int &, const int &,             // tgh0, tgh1, tgh2
   const int &,                            // cut_cell_tag
   const int &,                            // tag_value
   const int *,                            // flag
   int *);                                 // tags
#endif
}

using namespace pdat;
using namespace geom;
using namespace xfer;

/*************************************************************************
 *
 * Constructor and Destructor for SampleApp class.
 *
 ************************************************************************/

SampleApp::SampleApp(
   const string& object_name,
   Pointer<Database> input_db,
   Pointer<CartesianGridGeometry> grid_geom,
   Pointer<EmbeddedBoundaryGeometry> eb_geom,
   Pointer<VisItDataWriter> viz_writer)
{
   d_object_name = object_name;
   d_eb_geom = eb_geom;

   /*
    * Read in input values
    */
   getFromInput(input_db);

   /*
    * Create variable(s)
    */
   VariableDatabase* variable_db = VariableDatabase::getDatabase();

   Pointer<CellVariable<int> > mark_var =
      new CellVariable<int>("mark", 1);

   Pointer<VariableContext> cur_cxt = variable_db->getContext("CURRENT");
   IntVector<NDIM> zero_ghosts(0);
   d_mark_id = variable_db->registerVariableAndContext(mark_var,
         cur_cxt,
         zero_ghosts);

   viz_writer->registerPlotQuantity("MARK", "SCALAR", d_mark_id);

#ifdef USE_NONUNIFORM_LB
   Pointer<CellVariable<double> > weight_var =
      new CellVariable<double>("weight", 1);
   IntVector<NDIM> one_ghost(0);
   d_weight_id = variable_db->registerVariableAndContext(weight_var,
         cur_cxt,
         one_ghost);

   viz_writer->registerPlotQuantity("WEIGHT", "SCALAR", d_weight_id);
#else
   d_weight_id = -1;
#endif

   /*
    * Create RefineAlgorithm used for regridding
    */
   d_fill_new_level = new RefineAlgorithm();

   Pointer<RefineOperator<NDIM> > int_refine_op =
      grid_geom->lookupRefineOperator(mark_var, "CONSTANT_REFINE");

   d_fill_new_level->registerRefine(d_mark_id, // dst
      d_mark_id,                               // src
      d_mark_id,                               // scratch
      int_refine_op);

}

SampleApp::~SampleApp()
{
}

int SampleApp::getWorkloadIndex() const
{
   return d_weight_id;
}

void SampleApp::setMarkerOnPatch(
   const Pointer<Patch>& patch,
   const double time,
   const bool initial_time)
{
   Pointer<CellData<int> > mark = patch->getPatchData(d_mark_id);

   if (initial_time) {
      mark->fillAll(0);
   }

   // Patch lower, upper indices, dx
   const Pointer<CartesianPatchGeometry> patch_geom =
      patch->getPatchGeometry();
   const double* dx = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const Index ifirst = patch->getBox().lower();
   const Index ilast = patch->getBox().upper();

   double xyz[NDIM];
   int loc_i, i;
   double centroid[NDIM], radius, cell_rad_squared;
   for (i = 0; i < NDIM; i++) {
      centroid[i] = 0.;
      xyz[i] = 0.;
   }

   Box pbox = patch->getBox();
   for (CellIterator ib(pbox); ib; ib++) {

      // determine centroid of cell
      for (i = 0; i < NDIM; i++) {
         loc_i = ib() (i) - ifirst(i);
         xyz[i] = xlo[i] + ((double)(loc_i) + 0.5) * dx[i];
      }

      // Determine centroid/radius of moving,growing sphere
      radius = d_radius + time * d_radius_growth_rate;
      for (i = 0; i < NDIM; i++) {
         centroid[i] = d_centroid[i] + time * d_centroid_velocity[i];
      }

      // see if cell_centroid falls in sphere
      cell_rad_squared = 0.;
      for (i = 0; i < NDIM; i++) {
         cell_rad_squared += (xyz[i] - centroid[i]) * (xyz[i] - centroid[i]);
      }
      if (((radius * radius) - cell_rad_squared) > 0.) {
         (*mark)(ib()) = 1;
      }
   }

}

void SampleApp::setWeightOnPatch(
   const Pointer<Patch>& patch,
   const double time,
   const bool initial_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_weight_id >= 0);
#endif

   if (!patch->checkAllocated(d_weight_id)) {
      patch->allocatePatchData(d_weight_id);
   }
   Pointer<CellData<double> > weight = patch->getPatchData(d_weight_id);

   // Initiallly, all cells have weight = 1
   weight->fillAll(1.0);

   // Loop through cut cells on patch and assign weights.  Based on timings,
   // cut cells take about 116.7 times longer to compute than regular cells,
   // so thats the weight we use.
   int eb_index = d_eb_geom->getIndexCutCellDataId();
   if (patch->checkAllocated(eb_index)) {
      Pointer<IndexData<NDIM, CutCell, pdat::CellGeometry> > eboundary =
         patch->getPatchData(eb_index);

      for (IndexData<NDIM, CutCell, pdat::CellGeometry>::Iterator cc(*eboundary);
           cc;
           cc++) {
         Index ind = cc().getIndex();
         CellIndex cind(ind);
         Box box(ind, ind);
         if (!box.intersects(weight->getGhostBox())) {
            TBOX_ERROR("ouch");
         }
         (*weight)(cind) = 116.7;
         //(*weight)(cind) = 1.;
      }
   }

}

void SampleApp::tagCellsOnPatch(
   const Pointer<Patch>& patch,
   const double time,
   const int tag_index,
   const bool initial_time)
{
   int flag_index = d_eb_geom->getCellFlagDataId();
   Pointer<CellData<int> > tags = patch->getPatchData(tag_index);
   Pointer<CellData<int> > flag = patch->getPatchData(flag_index);

   tags->fillAll(0);
   Box pbox = patch->getBox();

   const Index ifirst = patch->getBox().lower();
   const Index ilast = patch->getBox().upper();

   IntVector<NDIM> tgh = tags->getGhostCellWidth();
   IntVector<NDIM> fgh = flag->getGhostCellWidth();

   /****************************************************************
    * Tag cells that are CUT
    ***************************************************************/

   if (d_tag_cut_cells) {

      int cut_cell_tag = EmbeddedBoundaryDefines::CUT;
      int tag_value = 1;

      F77_FUNC(tagcells, TAGCELLS) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         fgh(0), fgh(1),
#if (NDIM == 3)
         fgh(2),
#endif
         tgh(0), tgh(1),
#if (NDIM == 3)
         tgh(2),
#endif
         cut_cell_tag,
         tag_value,
         flag->getPointer(),
         tags->getPointer());
   }

   /****************************************************************
    * If desired, tag around directed sphere
    ***************************************************************/

   if (d_tag_growing_sphere) {

      Pointer<CellData<int> > mark = patch->getPatchData(d_mark_id);

      IntVector<NDIM> mgh = mark->getGhostCellWidth();

      int mark_tag = 1;
      int tag_value = 1;

      F77_FUNC(tagcells, TAGCELLS) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         mgh(0), mgh(1),
#if (NDIM == 3)
         mgh(2),
#endif
         tgh(0), tgh(1),
#if (NDIM == 3)
         tgh(2),
#endif
         mark_tag,
         tag_value,
         mark->getPointer(),
         tags->getPointer());

   } // if d_tag_growing_sphere

}

void SampleApp::printBoundaryNodeData(
   const Pointer<PatchLevel>& level,
   const Pointer<EmbeddedBoundaryGeometry>& eb_geom)
{

   // loop over patches on level
   for (PatchLevel::Iterator ip(level); ip; ip++) {
      Pointer<Patch> patch = level->getPatch(ip());

      // Pointers to eb_geometry cell, node flag data.
      Pointer<CellData<int> > cell_flag =
         patch->getPatchData(eb_geom->getCellFlagDataId());
      Pointer<NodeData<int> > node_flag =
         patch->getPatchData(eb_geom->getNodeInsideOutsideDataId());

      Pointer<IndexData<NDIM, CutCell, pdat::CellGeometry> > eboundary =
         patch->getPatchData(eb_geom->getIndexCutCellDataId());

      Index ij;
      for (IndexData<NDIM, CutCell, pdat::CellGeometry>::Iterator cc(*eboundary);
           cc;
           cc++) {

         ij = cc().getIndex();
         // print boundary node information for this cut cell
         if (appu::CutCell::boundaryNodesEnabled()) {
            cc().printBoundaryNodes(pout);

#if 0
            int num_bdry_nodes = cc().getNumberOfBoundaryNodes();
            const double* front_centroid = cc().getFrontCentroid();
            const double exposed_cut_surf_area = cc().getFrontArea();

            // Loop thru boundary nodes of cut cell.
            for (j = 0; j < num_bdry_nodes; j++) {
               BoundaryNode bn = cc().getBoundaryNode(j);

               NodeIndex bdry_node = bn.getIndex();
            }
#endif
         }

      } // loop over cut cells

   } // loop over patches
}

#ifdef SGS
void SampleApp::initializeLevelData(
   const Pointer<PatchHierarchy> hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const Pointer<PatchLevel> old_level,
   const bool allocate_data)
#endif

void SampleApp::initializeLevelData(
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer<hier::PatchLevel> old_level,
   const tbox::Pointer<hier::Connector>& to_old_mapped_box_level,
   const tbox::Pointer<hier::Connector>& from_old_mapped_box_level,
   const bool allocate_data)
{
   Pointer<PatchHierarchy> reg_hierarchy = hierarchy;

   Pointer<PatchLevel> level =
      hierarchy->getPatchLevel(level_number);

   /*
    * Allocate storage for local vars.
    */
   if (allocate_data) {
      level->allocatePatchData(d_mark_id, init_data_time);
   } else {
      level->setTime(init_data_time, d_mark_id);
   }

   /*
    * Create schedules and fill marker and weight data on new level
    */
   if ((level_number > 0) || (!old_level.isNull())) {

#ifdef SGS
      Pointer<RefineSchedule> sched =
         d_fill_new_level->createSchedule(level,
            old_level,
            level_number - 1,
            hierarchy,
            NULL);
#endif

      Pointer<RefineSchedule> sched =
         d_fill_new_level->createSchedule(level,
            old_level,
            level_number - 1,
            hierarchy,
            NULL);

      sched->fillData(init_data_time);
   }

   /*
    * Set the marker on the new level
    */
   if (d_tag_growing_sphere) {

      for (PatchLevel::Iterator ip(level); ip; ip++) {
         Pointer<Patch> patch = level->getPatch(ip());
         setMarkerOnPatch(patch, init_data_time, initial_time);
      }

   }

#ifdef USE_NONUNIFORM_LB
   /*
    * Set the weight on the new level
    */
   for (PatchLevel::Iterator ip(level); ip; ip++) {
      Pointer<Patch> patch = level->getPatch(ip());
      setWeightOnPatch(patch, init_data_time, initial_time);
   }
#endif

}

void SampleApp::resetHierarchyConfiguration(
   const Pointer<PatchHierarchy> hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   (void)hierarchy;
   (void)coarsest_level;
   (void)finest_level;
}

void SampleApp::applyGradientDetector(
   const Pointer<PatchHierarchy> hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   Pointer<PatchLevel> level = hierarchy->getPatchLevel(level_number);
   for (PatchLevel::Iterator ip(level); ip; ip++) {
      Pointer<Patch> patch = level->getPatch(ip());

      tagCellsOnPatch(patch, time, tag_index, initial_time);

   } // loop over patches
}

void
SampleApp::getFromInput(
   tbox::Pointer<tbox::Database> db)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   d_tag_cut_cells =
      db->getBoolWithDefault("tag_cut_cells", false);

   d_tag_growing_sphere =
      db->getBoolWithDefault("tag_growing_sphere", false);

   if (d_tag_growing_sphere) {
      d_centroid = db->getDoubleArray("centroid");
      d_centroid_velocity = db->getDoubleArray("centroid_velocity");
      d_radius = db->getDouble("radius");
      d_radius_growth_rate = db->getDouble("radius_growth_rate");
   }

}
