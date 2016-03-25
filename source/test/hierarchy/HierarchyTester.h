// 
// File:        HierarchyTester.h
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Manager class for patch hierarchy refine/coarsen tests.
//

#ifndef included_HierarchyTester
#define included_HierarchyTester

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_noprefix_Database
#include "tbox/Database.h"
#endif
#ifndef included_noprefix_IntVector
#include "IntVector.h"
#endif
#ifndef included_noprefix_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_noprefix_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_noprefix_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_noprefix_StandardTagAndInitStrategy
#include "StandardTagAndInitStrategy.h"
#endif
#ifndef included_tbox_String
#include <string>
using namespace std;
#define included_String
#endif

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace mesh;

namespace SAMRAI {
/**
 * Class HierarchyTester tests patch hierarchy coarsen/refine operations.
 * It sets up an initial patch hierarchy based on input file specifications.
 * Then, depending on how the test is specified in the input file, a 
 * new patch hierarchy will be created which is either a coarsened version
 * or a refined version of the initial hierarchy.  Then, all aspects of
 * of the new hierarchy are compared against the initial hierarchy to
 * check whether the hierarchy coarsenin or refining were performed 
 * correctly.
 * 
 * The functions in this class called from main() are:
 * \begin{enumerate}
 *    - [setupInitialHierarchy(...)] 
 *                           
 *    - [runHierarchyTestAndVerify(...)]
 * \end{enumerate}
 */

class HierarchyTester : public StandardTagAndInitStrategy<NDIM>
{
public:

   /**
    * Constructor initailizes test operations based on input database.
    */
   HierarchyTester(const string& object_name,
                   Pointer<Database> hier_test_db);

   /**
    * Destructor deallocates internal storage.
    */
   ~HierarchyTester();

   /**
    * Set up initial hierarchy used to test coarsen/refine operations.
    */
   void setupInitialHierarchy(Pointer<Database> main_input_db);

   /**
    * After hierarchy refine/coarsen operations are performed, check results.
    */
   void runHierarchyTestAndVerify();

   /**
    * The following two functions are declared pure virtual in the
    * StandardTagAndInitStrategy<NDIM> base class.  Although they do nothing,
    * they must be defined here for proper operation of the 
    * GriddingAlgorithm<NDIM> class which creates the initial patch hierarchy.
    */
   void initializeLevelData(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const Pointer<BasePatchLevel<NDIM> > old_level = Pointer<BasePatchLevel<NDIM> >(NULL),
      const bool allocate_data = true) 
   {
      NULL_USE(level_number);
      NULL_USE(init_data_time);
      NULL_USE(can_be_refined);
      NULL_USE(initial_time);
      NULL_USE(old_level);
      NULL_USE(allocate_data);
   }

   void resetHierarchyConfiguration(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int coarsest_level,
      const int finest_level)
   {
      NULL_USE(hierarchy);
      NULL_USE(coarsest_level);
      NULL_USE(finest_level);
   }

private:
   /*
    * Object name for error reporting.
    */
   string d_object_name;

   /*
    * Booleans to indicate whether refine or coarsen is operation to test.
    * They are initialized to false and set via input values.
    */
   bool d_do_refine_test;
   bool d_do_coarsen_test;

   /*
    * Ratio to coarsen/refine hierarchy during test.
    */
   IntVector<NDIM> d_ratio; 

   /*
    * Initial patch hierarchy set up based on input data and second hierarchy
    * generated during coarsen/refine operations.
    */
   Pointer<PatchHierarchy<NDIM> > d_initial_patch_hierarchy;
   Pointer<PatchHierarchy<NDIM> > d_test_patch_hierarchy;

};

}
#endif
