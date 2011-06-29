/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Class to manage functions for QM calculations. 
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

/*
 * Header files for SAMRAI classes referenced in this class.
 */
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/appu/EmbeddedBoundaryGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"

using namespace std;
using namespace SAMRAI;
using namespace appu;
using namespace tbox;
using namespace mesh;
using namespace geom;
using namespace xfer;
using namespace hier;

/**
 * This class provides basic level initialization and boundary condtions
 * for a simple application code.  It is intended to mimic a real app code
 * for tests with the EmbeddedBoundaryGeometry classes.
 */

class SampleApp:
   public StandardTagAndInitStrategy
{
public:
   /**
    * Default constructor for SampleApp.
    */
   SampleApp(
      const string& object_name,
      Pointer<Database> input_db,
      Pointer<CartesianGridGeometry> grid_geom,
      Pointer<EmbeddedBoundaryGeometry> eb_geom,
      Pointer<VisItDataWriter> viz_writer);

   /**
    * Empty destructor for SampleApp.
    */
   virtual ~SampleApp();

   /**
    * Returns workload index for nonuniform load balance.
    */
   int
   getWorkloadIndex() const;

   /*
    * Outputs boundary node information.
    */
   void
   printBoundaryNodeData(
      const Pointer<PatchLevel>& level,
      const Pointer<EmbeddedBoundaryGeometry>& eb_geom);

   /**
    * Initialize data on level.
    */
   virtual void
   initializeLevelData(
      const tbox::Pointer<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const tbox::Pointer<hier::PatchLevel> old_level =
         tbox::Pointer<hier::PatchLevel>(NULL),
      const tbox::Pointer<hier::Connector>& to_old_mapped_box_level =
         tbox::Pointer<hier::Connector>(NULL),
      const tbox::Pointer<hier::Connector>& from_old_mapped_box_level =
         tbox::Pointer<hier::Connector>(NULL),
      const bool allocate_data = true);

   /**
    * Reset hierarchy after regrid.
    */
   virtual void
   resetHierarchyConfiguration(
      const Pointer<PatchHierarchy> hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Tag cells for refinement using gradient detector.
    */
   void
   applyGradientDetector(
      const Pointer<PatchHierarchy> hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /**
    * Sets the "mark" variable used to identify where to place refinement
    * (called by initializeLevelData).
    */
   void
   setMarkerOnPatch(
      const Pointer<Patch>& patch,
      const double time,
      const bool initial_time);

   /**
    * Sets the "weight" variable used to assign workload weights on patch
    * cells (called by initializeLevelData).
    */
   void
   setWeightOnPatch(
      const Pointer<Patch>& patch,
      const double time,
      const bool initial_time);

   /**
    * Tags cells on a patch (called by gradient detector).
    */
   void
   tagCellsOnPatch(
      const Pointer<Patch>& patch,
      const double time,
      const int tag_index,
      const bool initial_time);
private:
   /*
    * Reads in input parameters from the specified database.
    */
   void
   getFromInput(
      Pointer<Database> database);

   /*
    * Object name
    */
   string d_object_name;

   /*
    * RefineAlgorithm used in regridding
    */
   Pointer<RefineAlgorithm> d_fill_new_level;

   /*
    * Pointer to EmbeddedBoundaryGeometry object.
    */
   Pointer<EmbeddedBoundaryGeometry> d_eb_geom;

   /*
    * Local data
    */
   int d_mark_id;
   int d_weight_id;

   /*
    * Tagging criteria
    */
   bool d_tag_growing_sphere;
   bool d_tag_cut_cells;
   Array<double> d_centroid;
   Array<double> d_centroid_velocity;
   double d_radius;
   double d_radius_growth_rate;

};
