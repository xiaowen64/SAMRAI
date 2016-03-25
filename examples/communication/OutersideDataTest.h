//
// File:        OutersideDataTest.h
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 317 $
// Modified:    $Date: 2005-04-27 21:26:06 -0700 (Wed, 27 Apr 2005) $
// Description: AMR communication tests for outerside-centered patch data
//

#ifndef included_OutersideDataTest
#define included_OutersideDataTest

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_Array
#include "tbox/Array.h"
#endif
#ifndef included_BoundaryBox
#include "BoundaryBox.h"
#endif
#ifndef included_Box
#include "Box.h"
#endif
#ifndef included_CartesianGridGeometry
#include "CartesianGridGeometry.h"
#endif
#ifndef included_CartesianPatchGeometry
#include "CartesianPatchGeometry.h"
#endif
#ifndef included_CellIndex
#include "CellIndex.h"
#endif
#ifndef included_SideData
#include "SideData.h"
#endif
#ifndef included_OutersideData
#include "OutersideData.h"
#endif
#ifndef included_Database
#include "tbox/Database.h"
#endif
#ifndef included_IntVector
#include "IntVector.h"
#endif
#ifndef included_Patch
#include "Patch.h"
#endif
#ifndef included_PatchDataTestStrategy
#include "PatchDataTestStrategy.h"
#endif
#ifndef included_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_Variable
#include "Variable.h"
#endif

namespace SAMRAI {

class CommTester;

/**
 * Class OutersideDataTest provides routines to test communication operations
 * for outerside-centered patch data on an AMR patch hierarchy.
 *
 * Required input keys and data types for test:
 *
 *   NONE...
 *
 * See PatchDataTestStrategy header file comments for variable and
 * refinement input data description.  Additionally, there are two
 * optional input parameters for each side variable.  These are:
 *
 * 


 *
 *    - \b  test_direction   side directions to test 
 *                             (default = -1 ie, all directions)
 *    - \b  use_fine_value_at_interface   which values to use at coarse-
 *                                          fine interface (default = TRUE)
 *
 * 


 * 
 */

class OutersideDataTest : public PatchDataTestStrategy
{
public:
  /**
   * The constructor initializes variable data arrays to zero length.
   */
   OutersideDataTest(const string& object_name,
                     tbox::Pointer<tbox::Database> main_input_db,
                     bool do_refine,
                     bool do_coarsen,
                     const string& refine_option);

   /**
    * Virtual destructor for SideDataTest.
    */
   ~OutersideDataTest();

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
   virtual void initializeDataOnPatch(
      hier::Patch<NDIM>& patch,
      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
      int level_number,
      char src_or_dst);

   /**
    * Function for checking results of communication operations.
    */
   void verifyResults(
      hier::Patch<NDIM>& patch,
      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
      int level_number);

private:
   /*
    * Function for reading test data from input file.
    */
   void readTestInput(tbox::Pointer<tbox::Database> db);

   void setLinearData(tbox::Pointer< pdat::OutersideData<NDIM,double> > data,
                      const hier::Box<NDIM>& box,
                      hier::Patch<NDIM>& patch) const;

   void setLinearData(tbox::Pointer< pdat::SideData<NDIM,double> > data,
                      const hier::Box<NDIM>& box,
                      hier::Patch<NDIM>& patch) const;


   void checkPatchInteriorData(
      const tbox::Pointer< pdat::OutersideData<NDIM,double> >& data,
      const hier::Box<NDIM>& interior,
      const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const;

   /*
    * Object string identifier for error reporting
    */
   string d_object_name;

   /*
    * Data members specific to this outerside data test.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_cart_grid_geometry;

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

   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > > d_variables_src;
   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > > d_variables_dst;

};

}
#endif
