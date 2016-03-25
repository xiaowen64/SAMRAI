//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/patch_configuration/PatchConfigurationTester.h $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Class to test patch configuration utility
//
 
#include "SAMRAI_config.h"

/*
 * Header files for SAMRAI classes referenced in this class.
 */
#include "Box.h"
#include "BoxList.h"
#include "tbox/InputDatabase.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/IOStream.h"

#include "PatchConfigurationUtilities.h"

using namespace std;
using namespace SAMRAI;

/**
 * The PatchConfigurationTester class tests the functionality of the 
 * hier::PatchConfigurationUtilities class.  Specifically, it builds a 
 * patch hierarchy, sets up a utilities object, and then tests the utilities
 * object data against the hierarchy configuration.  PatchConfigurationTester 
 * input data must follow this format:
 *
 * PatchConfigurationTester {
 *
 * Required 'test_to_run' input string:
 * 
 * test_to_run = <test option> 
 *    test option must be either "HIERARCHY_AT_ONCE", 
 *                               "LEVELS_COARSE_TO_FINE", or
 *                               "LEVELS_FINE_TO_COARSE"
 *   
 * Required GriddingParameters input describes AMR patch hierarchy configuration:
 *
 *    GriddingParameters { 
 *       // number of hierarchy levels [optional; default = 2]
 *       num_levels = 2
 *
 *       // integer vectors (length = NDIM) specifying ratio between index 
 *       // space of patch level to next coarser level 
 *       // [required if num_levels > 1] 
 *       ratio_to_coarser {
 *          level_1 = 2, 2, 2
 *          // all finer levels use level_1 values unless specified...
 *       }
 *    }
 *
 * }  
 */

class PatchConfigurationTester
{
public:

   /**
    * Constructor for PatchConfigurationTester.
    */     
   PatchConfigurationTester(const string& object_name,
                            tbox::Pointer<tbox::InputDatabase> input_db,
                            tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy);

   /**
    * Boring destructor for PatchConfigurationTester.
    */
   virtual ~PatchConfigurationTester();

   /**
    * Construct patch hierarchy based on input file information.
    */
   void buildPatchHierarchy();

   /**
    * Setup information about patch configuration.
    */
   void setupPatchConfiguration();

   /**
    * Check accuracy of patch configuration data against actual hierarchy
    * and return integer number of test failures.
    */
   int checkPatchConfiguration(ostream &os) const;

   /**
    * Dump the hierarchy configuration to the specified output stream.
    */
   void printHierarchyData(ostream &os) const;

private:
   /*
    * Read input parameters for patch hierarchy construction.
    */
   void getGriddingParametersFromInput(
      tbox::Pointer<tbox::InputDatabase> gridding_db);

   /*
    * Private member functions to check patch configuration data.
    */
   bool checkPatchNeighbors(int ip, int ln) const;
   bool comparePatchNeighbors(
      int ip,
      int ln,
      int codim,
      const tbox::Array< hier::PatchConfigurationUtilities<NDIM>::NeighborPatchInfo >&
         util_neighbors,
      const tbox::Array< hier::PatchConfigurationUtilities<NDIM>::NeighborPatchInfo >&
         test_neighbors) const;
   bool checkFinerLevelPatchOverlap(int ip, int ln) const;
   bool checkCoarserLevelPatchOverlap(int ip, int ln) const;

   /*
    * Object string name identifier for error reporting, etc.
    */
   string d_object_name;

   string d_test_to_run;
  
   /*
    * Objects storing AMR patch hierarchy information. 
    */ 
   tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;
   int d_num_levels;
   tbox::Array< hier::IntVector<NDIM> > d_ratio_to_coarser;
   hier::IntVector<NDIM> d_npatches_on_coarsest;

   tbox::Pointer< hier::PatchConfigurationUtilities<NDIM> > d_patch_config_util;

};
