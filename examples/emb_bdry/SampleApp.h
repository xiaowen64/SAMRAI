//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/emb_bdry/SampleApp.h $
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Class to manage functions for QM calculations.
//
 
#include "SAMRAI_config.h"

/*
 * Header files for SAMRAI classes referenced in this class.
 */
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "EmbeddedBoundaryGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "tbox/Pointer.h"
#include "tbox/IOStream.h"
#include "StandardTagAndInitStrategy.h"

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

class SampleApp 
   : public StandardTagAndInitStrategy<NDIM>
{
public:
   /**
    * Default constructor for SampleApp.
    */     
   SampleApp(const string& object_name,
             Pointer<Database> input_db,
             Pointer<CartesianGridGeometry<NDIM> > grid_geom,
             Pointer<EmbeddedBoundaryGeometry<NDIM> > eb_geom,
             Pointer<VisItDataWriter<NDIM> > viz_writer);

   /**
    * Empty destructor for SampleApp.
    */
   virtual ~SampleApp();

   /**
    * Returns workload index for nonuniform load balance.
    */
   int getWorkloadIndex() const;
   
   /*
    * Outputs boundary node information.
    */
   void printBoundaryNodeData(const Pointer<PatchLevel<NDIM> >& level,
                              const Pointer<EmbeddedBoundaryGeometry<NDIM> >& eb_geom);

   /**
    * Initialize data on level.
    */  
   void initializeLevelData(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int level_number,
      const double time,
      const bool can_be_refined,
      const bool initial_time,
      const Pointer<BasePatchLevel<NDIM> > old_level,
      const bool allocate_data);
   
   /**
    * Reset hierarchy after regrid.
    */
   virtual void resetHierarchyConfiguration(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Tag cells for refinement using gradient detector.
    */
   void applyGradientDetector(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /**
    * Sets the "mark" variable used to identify where to place refinement
    * (called by initializeLevelData).
    */
   void setMarkerOnPatch(const Pointer<Patch<NDIM> >& patch,
                         const double time,
                         const bool initial_time);

   /**
    * Sets the "weight" variable used to assign workload weights on patch
    * cells (called by initializeLevelData).
    */
   void setWeightOnPatch(const Pointer<Patch<NDIM> >& patch,
                         const double time,
                         const bool initial_time);
  
   /**
    * Tags cells on a patch (called by gradient detector).
    */
   void tagCellsOnPatch(const Pointer<Patch<NDIM> >& patch,
                        const double time,
                        const int tag_index,
                        const bool initial_time);
private:
 
   /*
    * Reads in input parameters from the specified database.
    */
   void getFromInput(Pointer<Database> database);
  

   /*
    * Object name
    */
   string d_object_name;

   /*
    * RefineAlgorithm used in regridding
    */
   Pointer<RefineAlgorithm<NDIM> > d_fill_new_level;

   /*
    * Pointer to EmbeddedBoundaryGeometry object.
    */
   Pointer<EmbeddedBoundaryGeometry<NDIM> > d_eb_geom;

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
