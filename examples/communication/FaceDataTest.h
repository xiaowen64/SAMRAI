//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/communication/FaceDataTest.h $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: AMR communication tests for face-centered patch data
//

#ifndef included_pdat_FaceDataTest
#define included_pdat_FaceDataTest

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_geom_CartesianGridGeometry
#include "CartesianGridGeometry.h"
#endif
#ifndef included_geom_CartesianPatchGeometry
#include "CartesianPatchGeometry.h"
#endif
#ifndef included_pdat_CellIndex
#include "CellIndex.h"
#endif
#ifndef included_pdat_FaceData
#include "FaceData.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchDataTestStrategy
#include "PatchDataTestStrategy.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif

namespace SAMRAI {

class CommTester;

/**
 * Class FaceDataTest provides routines to test communication operations
 * for face-centered patch data on an AMR patch hierarchy.
 *
 * Required input keys and data types for test:
 *
 *   NONE...
 *
 * See PatchDataTestStrategy header file comments for variable and
 * refinement input data description.  Additionally, there are two
 * optional input parameters for each face variable.  These are:
 *
 * 


 *
 *    - \b  test_direction   face directions to test 
 *                             (default = -1 ie, all directions)
 *    - \b  use_fine_value_at_interface   which values to use at coarse-
 *                                          fine interface (default = TRUE)
 *
 * 


 * 
 */

class FaceDataTest : public PatchDataTestStrategy
{
public:
  /**
   * The constructor initializes variable data arrays to zero length.
   */
   FaceDataTest(const string& object_name,
                tbox::Pointer<tbox::Database> main_input_db,
                bool do_refine,
                bool do_coarsen,
                const string& refine_option);

   /**
    * Virtual destructor for FaceDataTest.
    */
   ~FaceDataTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   virtual void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                              const double time,
                                              const hier::IntVector<NDIM>& gcw) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   void registerVariables(CommTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void initializeDataOnPatch(hier::Patch<NDIM>& patch,
                                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
				      char src_or_dst);

   /**
    * Function for checking results of communication operations.
    */
   bool verifyResults(hier::Patch<NDIM>& patch,
                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                      int level_number);

private:
   /*
    * Function for reading test data from input file.
    */
   void readTestInput(tbox::Pointer<tbox::Database> db);

   void setLinearData(tbox::Pointer< pdat::FaceData<NDIM,double> > data,
                      const hier::Box<NDIM>& box,
                      hier::Patch<NDIM>& patch) const;

   /*
    * Set conservative linear function data for testing coarsening
    */
   void setConservativeData(tbox::Pointer< pdat::FaceData<NDIM,double> > data,
                            const hier::Box<NDIM>& box,
                            hier::Patch<NDIM>& patch, 
                            const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
                            int level_number) const;

   void checkPatchInteriorData(
      const tbox::Pointer< pdat::FaceData<NDIM,double> >& data,
      const hier::Box<NDIM>& interior,
      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const;

   /*
    * Object string identifier for error reporting
    */
   string d_object_name;

   /*
    * Data members specific to this face data test.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_cart_grid_geometry;

   /*
    * Data members specific to this face data test.
    */
   tbox::Array<int> d_test_direction;
   tbox::Array<bool> d_use_fine_value_at_interface;

   double d_Acoef;
   double d_Bcoef;
   double d_Ccoef;
   double d_Dcoef;

   bool d_do_refine;
   bool d_do_coarsen;
   string d_refine_option;
   int d_finest_level_number;

   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > > d_variables;

};

}
#endif
