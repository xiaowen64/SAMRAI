//
// File:        AutoTester.h
// Package:     SAMRAI data transfer
// Copyright:   (c) 1997-2002 The Regents of the University of California
// Revision:    $Revision: 7 $
// Modified:    $Date: 2004-11-30 13:18:17 -0800 (Tue, 30 Nov 2004) $
// Description: Simple class used for autotesting. 
//

#ifndef included_AutoTesterXD
#define included_AutoTesterXD

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif


#include "tbox/Array.h"
#include "BoxArray.h"
#include "BoxIOUtility.h"
#include "tbox/Database.h"
#include "GriddingAlgorithm.h"
#include "HyperbolicLevelIntegrator.h"
#include "MethodOfLinesIntegrator.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#include "TimeRefinementIntegrator.h"


using namespace SAMRAI;

/**
 *  Class AutoTester sets and verifies certain testing information
 *  used for autotesting the applications codes. 
 *
 *  The following input parameters may be specified:
 *  


 *     - \b test_fluxes (bool) specifies whether we will do 
 *                 Riemann test or test on timesteps.                   
 *     - \b test_iter_num (int) iteration to carry out test.
 *     - \b correct_result (double array) holds correct result    
 *     - \b output_correct (bool) specifies whether we will write
 *                 correct result (useful if changing problems 
 *                 and want to set "correct" array).  
 *  

         
 */

class AutoTester
{
public:

   /**
    * Default constructor for AutoTester
    */
   AutoTester(const string& object_name,
              tbox::Pointer<tbox::Database> input_db);

   /**
    * Virtual destructor for AutoTester
    */
   virtual ~AutoTester();


   /**
    * Checks result for applications using TimeRefinementIntegrator
    * and HyperbolicLevelIntegrator.
    */
   void evalTestData(int iter,
                     const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                     const tbox::Pointer<algs::TimeRefinementIntegrator<NDIM> > tri,
                     const tbox::Pointer<algs::HyperbolicLevelIntegrator<NDIM> > hli,
                     const tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > ga);

   /**
    * Checks result for applications using MethodOfLinesIntegrator
    */
   void evalTestData(int iter,
                     const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                     const double time,
                     const tbox::Pointer<algs::MethodOfLinesIntegrator<NDIM> > mol,
                     const tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > ga);

   /**
    * Check boxes on specified hierarchy level against test boxes.
    */
   void checkHierarchyBoxes(const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                            int ln,
                            const hier::BoxArray<NDIM>& test_boxes,
                            int iter); 

private:
   /*
    *  Sets the parameters in the struct, based 
    *  on data read from input.
    */
   void getFromInput(tbox::Pointer<tbox::Database> input_db);

   string d_object_name;

   bool   d_test_fluxes;  //  specifies whether to check Riemann problem
   bool   d_output_correct; // output  result?
   int    d_test_iter_num;  // iteration number to check result.

   tbox::Array<double> d_correct_result;  // array to hold correct values

   bool d_write_patch_boxes;
   bool d_read_patch_boxes;
   tbox::Array<int> d_test_patch_boxes_at_steps;
   string d_test_patch_boxes_filename;
   bool d_simplify_test_boxes;

   hier::BoxIOUtility<NDIM>* d_io_box_utility;
   int d_test_patch_boxes_step_count;
   
};

#endif

