//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/mblkcomm/CellMultiblockTest.h $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: AMR communication tests for cell-centered patch data
//

#ifndef included_CellMultiblockTest
#define included_CellMultiblockTest

#include "SAMRAI_config.h"

#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "Box.h"
#include "CellData.h"
#include "BlockGridGeometry.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchMultiblockTestStrategy.h"

#include "BlockPatchGeometry.h"
#include "Variable.h"

using namespace SAMRAI;

/**
 * Class CellMultiblockTest provides routines to test communication operations
 * for cell-centered patch data on an AMR patch hierarchy.
 *
 * See PatchMultiblockTestStrategy header file comments for variable and
 * refinement input data description.
 */

class CellMultiblockTest : public PatchMultiblockTestStrategy
{
public:
  /**
   * The constructor initializes variable data arrays to zero length.
   */
   CellMultiblockTest(const string& object_name,
                tbox::Pointer<tbox::Database> main_input_db,
                bool do_refine,
                bool do_coarsen,
                const string& refine_option);

   /**
    * Virtual destructor for CellMultiblockTest.
    */
   virtual ~CellMultiblockTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   void setPhysicalBoundaryConditions(
      hier::Patch<NDIM>& patch,
      const double time,
      const hier::IntVector<NDIM>& gcw_to_fill) const;

   void fillSingularityBoundaryConditions(
      hier::Patch<NDIM>& patch,
      tbox::List<xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
         sing_patches,
      const hier::Box<NDIM>& fill_box,
      const hier::BoundaryBox<NDIM>& boundary_box);

   /**
    * This function is called from the MultiblockTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the MultiblockTester object, which
    * manages the variable storage.
    */
   void registerVariables(MultiblockTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void initializeDataOnPatch(
      hier::Patch<NDIM>& patch,
      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
      int level_number,
      int block_number,
      char src_or_dst);

   /**
    * Function for tagging cells on each patch to refine.
    */
   void tagCellsToRefine(
      hier::Patch<NDIM>& patch,
      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
      int level_number,
      int tag_index);
   
   /**
    * Function for checking results of communication operations.
    */
   bool verifyResults(
      hier::Patch<NDIM>& patch,
      const tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> > hierarchy,
      int level_number,
      int block_number);

   ///
   void postprocessRefine(hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const tbox::Pointer<hier::VariableContext>& context,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio) const;

private:
   /**
    * Function for reading test data from input file.
    */
   void readTestInput(tbox::Pointer<tbox::Database> db);

   /*
    * Object string identifier for error reporting
    */
   string d_object_name;

   /*
    * Data members specific to this cell data test.
    */
   tbox::Array< tbox::Pointer< hier::GridGeometry<NDIM> > > d_skel_grid_geometry;

   string d_refine_option;
   int d_finest_level_number;

   tbox::Array< tbox::Pointer< hier::Variable<NDIM> > > d_variables;

};

#endif
