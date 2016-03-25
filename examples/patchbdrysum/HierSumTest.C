//
// File:        Hierarchy Sum test
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:     $Name$
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: SAMRAI interface class for hierarchy node and edge sum test
//

#include "HierSumTest.h"


#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "EdgeData.h"
#include "Index.h"
#include "NodeData.h"
#include "NodeIndex.h"
#include "NodeIterator.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"
#include "PatchBoundaryNodeSum.h"
#include "PatchBoundaryEdgeSum.h"
#include "CoarseFineBoundary.h" 

extern "C" {
#if (NDIM == 2)
   void setedges2d_(const int&, const int&, 
                    const int&, const int&,
                    const int&, const int&,
                    const int&, const int&,
                    const double*,
                    double*,
                    double*);

   void checkedges2d_(const int&, const int&, 
                      const int&, const int&,
                      const int&, const int&,
                      const double&, 
                      int&, 
                      const double*,
                      const double*);
#endif
#if (NDIM == 3)
   void setedges3d_(const int&, const int&, const int &, 
                    const int&, const int&, const int &,
                    const int&, const int&, const int &,
                    const int&, const int&, const int &,
                    const double*,
                    double*,
                    double*,
                    double*);

   void checkedges3d_(const int&, const int&, const int &,
                      const int&, const int&, const int &,
                      const int&, const int&, const int &,
                      const double&, 
                      int&, 
                      const double*,
                      const double*,
                      const double*);
#endif
}

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace pdat;
using namespace geom;
using namespace algs;
using namespace appu;

/*************************************************************************
 *
 * Constructor and Destructor.
 *
 ************************************************************************/

HierSumTest::HierSumTest(
   const string& object_name,
   Pointer<Database> input_db,
   Pointer<VisItDataWriter<NDIM> > viz_writer)
{

   d_object_name = object_name;

   d_depth = 1;
   
   /*
    * Initialize object with data read from given input databases.
    */
   getFromInput(input_db);
   
   /* 
    * Set up variables and contexts. 
    *
    * Vars:
    *    ucell - cell-centered variable U
    *    unode - node-centered variable U
    *    uedge - edge-centered variable U
    *    
    * Contexts:
    *
    *    NOGHOST  - zero ghosts (unode)
    *    ONEGHOST - one ghost (ucell)
    *
    * What the test case does:
    *   1. Set ucell = 1.0 on cells of all levels
    *   2. Set ucell = 0.0 on cells of L < LN that are
    *      covered by refined cells.
    *   3. Set node weight on patch interiors = sum(cell weights)
    *   3. Do a hier sum transaction
    *   4. Correct result - all nodes on all levels = 2^NDIM 
    *
    *
    * Below we construct u variable and its contexts.
    */
   VariableDatabase<NDIM>* variable_db = VariableDatabase<NDIM>::getDatabase();

   d_ucell_var = new CellVariable<NDIM,double>("ucell",d_depth);
   d_unode_var = new NodeVariable<NDIM,double>("unode",d_depth);
   d_uedge_var = new EdgeVariable<NDIM,double>("uedge",d_depth);
   
   Pointer<VariableContext> cxt1 = variable_db->getContext("CONTEXT1");
   Pointer<VariableContext> cxt2 = variable_db->getContext("CONTEXT2");

   IntVector<NDIM> one_ghost(1);

   d_ucell_node_id = 
      variable_db->registerVariableAndContext(d_ucell_var,
                                              cxt1,
                                              one_ghost);

   d_unode_id = variable_db->registerVariableAndContext(d_unode_var,
                                                        cxt1,
                                                        d_node_ghosts);
   d_ucell_edge_id = 
      variable_db->registerVariableAndContext(d_ucell_var,
                                              cxt2,
                                              one_ghost);

   d_uedge_id = variable_db->registerVariableAndContext(d_uedge_var,
                                                        cxt1,
                                                        d_edge_ghosts);

   /*
    * Register u values to be written by the viz writer.
    */
   viz_writer->registerPlotQuantity("ucell::node", "SCALAR", 
                                    d_ucell_node_id, 0);
   viz_writer->registerPlotQuantity("ucell::edge", "SCALAR", 
                                    d_ucell_edge_id, 0);
   viz_writer->registerPlotQuantity("unode", "SCALAR", d_unode_id);


}

HierSumTest::~HierSumTest()
{   
}

/*************************************************************************
 *
 * Set initial node values based on sum of surrounding cell weights
 *
 ************************************************************************/
int
HierSumTest::setInitialNodeValues(
   const Pointer<PatchHierarchy<NDIM> > hierarchy) 
{
   int fail_count = 0;

   /*
    * Set node weight on patch interiors = sum(cell weights)
    */
   // loop over hierarchy levels   
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      // loop over patches on level
      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(ip());  

         Pointer< NodeData<NDIM,double> > unode = patch->getPatchData(d_unode_id);
         Pointer< CellData<NDIM,double> > ucell = 
            patch->getPatchData(d_ucell_node_id);

         // output initial cell values
         int level_number = level->getLevelNumber();
         tbox::plog << "INITIAL Cell values for NODE - Level: " << level_number
              << "\tPatch: " << patch->getPatchNumber()  << endl;
         ucell->print(ucell->getGhostBox(),plog);

         // loop over nodes of patch
         Box<NDIM> pbox = patch->getBox();
         double cell_val, node_val;
         int d;
         for (NodeIterator<NDIM> ni(pbox); ni; ni++) {
            NodeIndex<NDIM> node = ni();
            for (d = 0; d < ucell->getDepth(); d++) {
               
               (*unode)(node,d) = 0.;

               /*
                * Sum contributions from surrounding cells for each node
                * value.
                */
#if (NDIM == 3)
               for (int k = 0; k <= 1; k++) {
#endif
                  for (int j = 0; j <= 1; j++) {
                     for (int i = 0; i <= 1; i++) {
                        CellIndex<NDIM> cell = ni();
                        cell(0) -= i;
                        cell(1) -= j;
#if (NDIM == 3)
                        cell(2) -= k;
#endif
                        cell_val = (*ucell)(cell,d);
                        node_val = (*unode)(node,d);
                        (*unode)(node,d) += cell_val;
                     }
                  }
#if (NDIM == 3)
               }
#endif
            } // loop over depth
         } // loop over nodes
      } // loop over patches

      /*
       * Any nodes that are *inside* the complement region (inside
       * meaning at least one cell away from the coarse-fine boundary) 
       * should not contribute to the nodal sum on finer levels.  Here
       * we verify this is the case by resetting all nodes inside this
       * region to -999.  
       * 
       * The so-called fine_overlap_shrunk is computed by finding the
       * region of the level that is not overlapped by fine
       * patches (complement), growing it, and removing intersection
       * with the level boxes.  This is effectively a shrunken coarse-fine
       * overlap region, on which the nodes shouldn't participate in any
       * communication. 
       */
      BoxArray<NDIM> level_boxes = level->getBoxes();
      BoxList<NDIM> fine_overlap_shrunk(level_boxes);
      BoxList<NDIM> complement(level_boxes);
      if (level->getLevelNumber() != hierarchy->getFinestLevelNumber()) {
         Pointer<PatchLevel<NDIM> > fine_level = hierarchy->getPatchLevel(ln+1);
         BoxArray<NDIM> fine_level_boxes = fine_level->getBoxes();
         fine_level_boxes.coarsen(fine_level->getRatioToCoarserLevel());
         BoxList<NDIM> fine_level_bl(fine_level_boxes);
         complement.removeIntersections(fine_level_bl);
         complement.grow(IntVector<NDIM>(1));
         fine_overlap_shrunk.removeIntersections(complement);

         for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
            Pointer<Patch<NDIM> > patch = level->getPatch(ip());  

            Pointer< NodeData<NDIM,double> > unode = 
               patch->getPatchData(d_unode_id);

            for (BoxList<NDIM>::Iterator b(fine_overlap_shrunk); b; b++) {
               Box<NDIM> fine_overlap = b();
               Box<NDIM> patch_interior = patch->getBox();
               Box<NDIM> data_box = fine_overlap * patch_interior;
               for (NodeIterator<NDIM> ni(data_box); ni; ni++) {
                  NodeIndex<NDIM> node = ni();
                  for (int d = 0; d < unode->getDepth(); d++) {
                     double node_val =  (*unode)(node,d);
                     if (tbox::MathUtilities<double>::equalEps(node_val, 0.0)) {
                        (*unode)(node,d) = -999.;
                     }
                  }
               }         
            } // loop over complement boxes
         } // loop over patches
      } // if a finer level exists

      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(ip());

         Pointer< NodeData<NDIM,double> > unode = 
            patch->getPatchData(d_unode_id);

         // output initial node values
         tbox::plog << "INITIAL Node values - Level: " << level->getLevelNumber()
              << "\tPatch: " << patch->getPatchNumber() << endl;
         unode->print(unode->getGhostBox(),plog);
         
      } // loop over patches
   } // loop over levels

   return(fail_count);
}

/*************************************************************************
 *
 * Set initial edge values based on sum of surrounding cell weights
 *
 ************************************************************************/
int
HierSumTest::setInitialEdgeValues(
   const Pointer<PatchLevel<NDIM> > level) 
{
   int fail_count = 0;

   double correct_val = 1.;
   for (int i = 0; i < NDIM-1; i++) {
      correct_val *= 2.;
   }
   
   /*
    * Set edge weight on patch interiors = sum(cell weights)
    */
   for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
      Pointer<Patch<NDIM> > patch = level->getPatch(ip());  
      
      Pointer< EdgeData<NDIM,double> > uedge = patch->getPatchData(d_uedge_id);
      Pointer< CellData<NDIM,double> > ucell = 
         patch->getPatchData(d_ucell_edge_id);

      // output initial cell values
      int level_number = level->getLevelNumber();
      tbox::plog << "INITIAL Cell values for EDGE - Level: " << level_number
           << "\tPatch: " << patch->getPatchNumber()  << endl;
      ucell->print(ucell->getGhostBox(),plog);

      const Index<NDIM> ifirst=patch->getBox().lower();
      const Index<NDIM> ilast =patch->getBox().upper();

      IntVector<NDIM> cellg = ucell->getGhostCellWidth();
      IntVector<NDIM> edgeg = uedge->getGhostCellWidth();

      for (int d = 0; d < uedge->getDepth(); d++) {      

#if (NDIM == 2)
         setedges2d_(ifirst(0),ifirst(1),
                     ilast(0),ilast(1),
                     cellg(0),cellg(1),
                     edgeg(0),edgeg(1),
                     ucell->getPointer(d),
                     uedge->getPointer(0,d),
                     uedge->getPointer(1,d));
#endif  
#if (NDIM == 3)
         setedges3d_(ifirst(0),ifirst(1),ifirst(2),
                     ilast(0),ilast(1),ilast(2),
                     cellg(0),cellg(1),cellg(2),
                     edgeg(0),edgeg(1),edgeg(2),
                     ucell->getPointer(d),
                     uedge->getPointer(0,d),
                     uedge->getPointer(1,d),
                     uedge->getPointer(2,d));
#endif  

         /*
          * If you want to check edges BEFORE communication (to make sure 
          * the communication is actually doing something) then do the 
          * check here.  Be forwarned that it dumps a lot of data because,
          * if things are working right, the data before communication will
          * have many errors, that the communication will fix.
          */
         
         if (d_check_data_before_communication) {
            
            int fort_all_correct = 1;
            
#if (NDIM == 2)
            checkedges2d_(ifirst(0),ifirst(1),
                          ilast(0),ilast(1),
                          edgeg(0),edgeg(1),
                          correct_val,
                          fort_all_correct,
                          uedge->getPointer(0,d),
                          uedge->getPointer(1,d));
#endif
#if (NDIM == 3)
            checkedges3d_(ifirst(0),ifirst(1),ifirst(2),
                          ilast(0),ilast(1),ilast(2),
                          edgeg(0),edgeg(1),edgeg(2),
                          correct_val,
                          fort_all_correct,
                          uedge->getPointer(0,d),
                          uedge->getPointer(1,d),
                          uedge->getPointer(2,d));
#endif
            
            if (fort_all_correct == 0) {
               fail_count++; 
               tbox::perr << "PatchBdrySum Edge test FAILED:  Errors on Level: " 
                          << level->getLevelNumber()
                          << "\t Patch: " << patch->getPatchNumber()
                          << "\t" << patch->getBox()
                          << "\nAll edges are not correct value." << endl;
            } else {
#if (TESTING == 1)
               tbox::plog 
#else
               tbox::pout
#endif 
                  << "All edges on Level: " << level->getLevelNumber()
                  << "\t Patch: " << patch->getPatchNumber()
                  << "\tare correct." << endl;
            }
         }

      } // loop over depth

#if (TESTING == 1)
      tbox::plog << "INITIAL Edge values - Level: " << level->getLevelNumber()
           << "\tPatch: " << patch->getPatchNumber() << endl;
      uedge->print(uedge->getGhostBox(),plog);
#endif
      
   } // loop over patches

   tbox::SAMRAI_MPI::barrier();

   return(fail_count);
   
}

/*************************************************************************
 *
 * Setup outer node sum.
 *
 ************************************************************************/
void
HierSumTest::setupOuternodeSum(
   const Pointer<PatchHierarchy<NDIM> > hierarchy) 
{
   d_node_sum_util = new PatchBoundaryNodeSum<NDIM>("Node Sum Util");

   d_node_sum_util->registerSum(d_unode_id);

   int num_levels = hierarchy->getNumberLevels();
   if (num_levels > 1) {

      int coarsest_level_number = 0;
      d_node_sum_util->setupSum(hierarchy, 
                                coarsest_level_number,
                                hierarchy->getFinestLevelNumber());
   } else {
      d_node_sum_util->setupSum(hierarchy->getPatchLevel(0));
   }
      
}

/*************************************************************************
 *
 * Perform outer node sum.
 *
 ************************************************************************/
void
HierSumTest::doOuternodeSum() 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_node_sum_util.isNull());
#endif
   
   bool fill_hanging_nodes = true;
   d_node_sum_util->computeSum(fill_hanging_nodes);
      
}

/*************************************************************************
 *
 * Setup Outeredge sum.
 *
 ************************************************************************/
void
HierSumTest::setupOuteredgeSum(
   const Pointer<PatchHierarchy<NDIM> > hierarchy,
   const int level_num)
{
   if (level_num >= d_edge_sum_util.getSize()) {
      d_edge_sum_util.resizeArray(level_num+1);
   }
   
   d_edge_sum_util[level_num] = 
      new PatchBoundaryEdgeSum<NDIM>("Level Edge Sum Util");

   d_edge_sum_util[level_num]->registerSum(d_uedge_id);

   d_edge_sum_util[level_num]->setupSum(hierarchy->getPatchLevel(level_num));
}


/*************************************************************************
 *
 * Setup and perform Outeredge sum.
 *
 ************************************************************************/
void
HierSumTest::doOuteredgeSum(const int level_num) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level_num < d_edge_sum_util.getSize());
   TBOX_ASSERT(!d_edge_sum_util[level_num].isNull());
#endif

   d_edge_sum_util[level_num]->computeSum();
      
}


/*************************************************************************
 *
 * Check correctness of hierarchy node sum operation
 *
 ************************************************************************/

int HierSumTest::checkNodeResult(
   const Pointer<PatchHierarchy<NDIM> > hierarchy)
{

   int fail_count = 0;

   /*
    * After the communication the sum on every node of every level should be
    * 2^NDIM.  Check this here...
    */

   double correct_val = 1.;
   for (int i = 0; i < NDIM; i++) {
      correct_val *= 2.;
   }

   // loop over hierarchy levels   
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      
      BoxArray<NDIM> level_boxes = level->getBoxes();
      BoxList<NDIM> level_boxes_complement(level_boxes);
      
      // If a finer level exists, remove overlap boxes by computing complement
      if (level->getLevelNumber() != hierarchy->getFinestLevelNumber()) {
         Pointer<PatchLevel<NDIM> > fine_level = hierarchy->getPatchLevel(ln+1);
         BoxArray<NDIM> fine_level_boxes = fine_level->getBoxes();
         fine_level_boxes.coarsen(fine_level->getRatioToCoarserLevel());
         BoxList<NDIM> fine_level_bl(fine_level_boxes);
         level_boxes_complement.removeIntersections(fine_level_bl);
      }
         
      // loop over patches on level 
      bool all_correct = true;
      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(ip());   

         Pointer< NodeData<NDIM,double> > unode = 
            patch->getPatchData(d_unode_id);

         // loop over Level complement boxlist
         for (BoxList<NDIM>::Iterator b(level_boxes_complement); b; b++) {
            Box<NDIM> complement = b();

            // intersect patch box with level box complement
            Box<NDIM> patch_interior = patch->getBox();
            Box<NDIM> data_box = complement * patch_interior;
            
            /*
             * Iterate over nodes and check correctness of result.
             */
            for (NodeIterator<NDIM> i(data_box); i; i++) {
               NodeIndex<NDIM> node = i();  // i,j
               for (int d = 0; d < unode->getDepth(); d++) {
                  
                  bool node_correct = false;
                  double node_val = (*unode)(node,d);
               
                  if (node_val == correct_val) node_correct = true;
                  if (!node_correct) {
                     tbox::pout << "BAD NODE = " << node_val << " at index " << i() 
                          << " depth = " << d << endl;
                     all_correct = false;
                     break;
                  }
               }
            }
         } // loop over complement boxes

         if (!all_correct) {
            fail_count++;
            tbox::perr << "PatchBdrySum Node test FAILED:  Errors on Level: " 
                       << level->getLevelNumber()
                       << "\t Patch: " << patch->getPatchNumber()
                       << "\t" << patch->getBox()
                       << "\nAll nodes are not correct value." << endl;
         } else {
#if (TESTING == 1)
            tbox::plog 
#else
               tbox::pout
#endif 
                  << "All nodes on Level: " << level->getLevelNumber()
                  << "\t Patch: " << patch->getPatchNumber()
                  << "\tare correct." << endl;
         }

         Pointer< CellData<NDIM,double> > ucell_node = 
            patch->getPatchData(d_ucell_node_id);
         
#if (TESTING == 1)
         tbox::plog << "FINAL Cell values for NODE - Level: " 
                    << level->getLevelNumber()
                    << "\tPatch: " << patch->getPatchNumber() << endl;
         ucell_node->print(ucell_node->getGhostBox(),plog);
         
         tbox::plog << "FINAL Node values - Level: " << level->getLevelNumber()
                    << "\tPatch " << patch->getPatchNumber() 
                    << " : " << patch->getBox() << endl;
         unode->print(unode->getBox(),plog);
#endif
         
      } // loop over patches
      
   } // loop over levels

   return(fail_count);
   
}


/*************************************************************************
 *
 * Check correctness of level edge sum operation
 *
 ************************************************************************/

int HierSumTest::checkEdgeResult(
   const Pointer<PatchLevel<NDIM> > level)
{

   int fail_count = 0;

   /*
    * After the communication the sum on every edge of every level should be
    * 2^(NDIM-1).  Check this here...
    */

   double correct_val = 1.;
   for (int i = 0; i < NDIM-1; i++) {
      correct_val *= 2.;
   }

   // loop over patches on level
   for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
      Pointer<Patch<NDIM> > patch = level->getPatch(ip());   
      
      Pointer< EdgeData<NDIM,double> > uedge = 
         patch->getPatchData(d_uedge_id);

      const Index<NDIM> ifirst=patch->getBox().lower();
      const Index<NDIM> ilast =patch->getBox().upper();

      IntVector<NDIM> edgeg = uedge->getGhostCellWidth();

      for (int d = 0; d < uedge->getDepth(); d++) {
         
         /*
          * In the fortran, set "fort_all_correct" to 0 if we
          * detect differences between uedge data and "correct_val".
          */
         int fort_all_correct = 1;
      
#if (NDIM == 2)
         checkedges2d_(ifirst(0),ifirst(1),
                       ilast(0),ilast(1),
                       edgeg(0),edgeg(1),
                       correct_val,
                       fort_all_correct,
                       uedge->getPointer(0,d),
                       uedge->getPointer(1,d));
#endif
#if (NDIM == 3)
         checkedges3d_(ifirst(0),ifirst(1),ifirst(2),
                       ilast(0),ilast(1),ilast(2),
                       edgeg(0),edgeg(1),edgeg(2),
                       correct_val,
                       fort_all_correct,
                       uedge->getPointer(0,d),
                       uedge->getPointer(1,d),
                       uedge->getPointer(2,d));
#endif


         if (fort_all_correct == 0) {
            fail_count++;
            tbox::perr << "PatchBdrySum Edge test FAILED:  Errors on Level: " 
                       << level->getLevelNumber()
                       << "\t Patch: " << patch->getPatchNumber()
                       << "\t" << patch->getBox()
                       << "\nAll edges are not correct value." << endl;
         } else {
#if (TESTING == 1)
            tbox::plog 
#else
            tbox::pout
#endif 
               << "All edges on Level: " << level->getLevelNumber()
               << "\t Patch: " << patch->getPatchNumber()
               << "\t" << patch->getBox()
               << "\tare correct." << endl;
         }
         
      } // loop over depth 
      
      Pointer< CellData<NDIM,double> > ucell_edge = 
         patch->getPatchData(d_ucell_edge_id);
      
#if (TESTING == 1)
      tbox::plog << "FINAL Cell values for EDGE - Level: " 
                 <<  level->getLevelNumber()
                 << "\tPatch: " << patch->getPatchNumber()  << endl;
      ucell_edge->print(ucell_edge->getGhostBox(),plog);      
      
      tbox::plog << "FINAL Edge values - Level: " 
                 << level->getLevelNumber()
                 << "\tPatch: " << patch->getPatchNumber() << endl;
      uedge->print(uedge->getGhostBox(),plog);
#endif
      
   } // loop over patches

   return(fail_count);

}

 
 

      

/*************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy<NDIM>.
 *
 ************************************************************************/


/*
 * Allocate storage and initialize data on the level.
 */
void
HierSumTest::initializeLevelData(
   const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const Pointer<BasePatchLevel<NDIM> > old_level,
   const bool allocate_data)
{

   Pointer<PatchHierarchy<NDIM> > local_hierarchy = hierarchy;

   /*
    * Set initial data on hierarchy level. 
    *   1. Set ucell = 1.0 on cells of all levels
    *   2. For NODE data only, set ucell = 0.0 on cells of L < LN that are
    *      covered by refined cells.
    */
   Pointer<PatchLevel<NDIM> > level =
      hierarchy->getPatchLevel(level_number);
   

   /*
    * Allocate storage for cell and node data.  
    */
   if (allocate_data) {
      level->allocatePatchData(d_ucell_node_id, time);
      level->allocatePatchData(d_ucell_edge_id, time);
      level->allocatePatchData(d_unode_id, time);
      level->allocatePatchData(d_uedge_id, time);
   }
   
   /*
    * Set edge/node values to zero initially.
    */
   for (PatchLevel<NDIM>::Iterator p0(level); p0; p0++) {
      Pointer<Patch<NDIM> > patch = level->getPatch(p0());

      Pointer< NodeData<NDIM,double> > unode = 
         patch->getPatchData(d_unode_id);
      Pointer< EdgeData<NDIM,double> > uedge = 
         patch->getPatchData(d_uedge_id);
      unode->fillAll(0.0);
      uedge->fillAll(0.0);
   }
   
      


   /*
    * Set cell weights.  We want interior cells set to 1.0 and 
    * ghost cells set to 0.0.  (Eventually, we will reset the ghosts
    * based on overlap with neighboring patches but this will be
    * in the next step).
    */
   for (PatchLevel<NDIM>::Iterator p0(level); p0; p0++) {
      Pointer<Patch<NDIM> > patch = level->getPatch(p0());

      Pointer< CellData<NDIM,double> > ucell_node = 
         patch->getPatchData(d_ucell_node_id);
      Pointer< CellData<NDIM,double> > ucell_edge = 
         patch->getPatchData(d_ucell_edge_id);

      ucell_node->fillAll(0.0, ucell_node->getGhostBox()); // ghost box
      ucell_node->fillAll(1.0, patch->getBox());          // interior patch box

      ucell_edge->fillAll(0.0, ucell_edge->getGhostBox()); // ghost box
      ucell_edge->fillAll(1.0, patch->getBox());          // interior patch box

      // set cell values at physical boundary
      const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = 
         patch->getPatchGeometry();
      const tbox::Array<BoundaryBox<NDIM> > node_bdry =
         patch_geom->getCodimensionBoundary(NDIM);
      const tbox::Array<BoundaryBox<NDIM> > edge_bdry =
         patch_geom->getCodimensionBoundary(NDIM-1);
#if (NDIM == 3)
      const tbox::Array<BoundaryBox<NDIM> > face_bdry =
         patch_geom->getCodimensionBoundary(1);
#else
      const tbox::Array<BoundaryBox<NDIM> > face_bdry;
#endif
      // node cell values
      setBoundaryConditions(*patch,
                            node_bdry,
                            edge_bdry,
                            face_bdry,
                            d_ucell_node_id);

      // edge cell values
      setBoundaryConditions(*patch,
                            node_bdry,
                            edge_bdry,
                            face_bdry,
                            d_ucell_edge_id);

   }
      

   if (level_number > 0) {

      /*
       * For node data, set the cell weights to zero on coarser level 
       * where there is overlap with fine level patches.
       */
      Pointer<PatchLevel<NDIM> > coarser_level = 
         hierarchy->getPatchLevel(level_number-1);
      BoxArray<NDIM> fine_level_boxes = level->getBoxes();


      IntVector<NDIM> ratio = level->getRatioToCoarserLevel();
      fine_level_boxes.coarsen(ratio);
      
      for (PatchLevel<NDIM>::Iterator p1(coarser_level); p1; p1++) {
         Pointer<Patch<NDIM> > cpatch = coarser_level->getPatch(p1());
         
         Pointer< CellData<NDIM,double> > ucell_node = 
            cpatch->getPatchData(d_ucell_node_id);

         Box<NDIM> cpbox = cpatch->getBox();
         for (int n = 0; n < fine_level_boxes.getNumberOfBoxes(); n++) {
            Box<NDIM> setbox = cpbox * fine_level_boxes(n);
            if (!setbox.empty()) {
               ucell_node->fillAll(0.0, setbox);
            }
         }

         // zero out cells on the boundary that lie at coarse-fine
         // interface.
         zeroOutPhysicalBoundaryCellsAtCoarseFineBoundary(*cpatch,
                                                          d_ucell_node_id);
         
      } // loop over coarser level patches

      /*
       * For edge data, set the ghosts of the cell data equal to 1.0 at the 
       * coarse-fine boundaries.
       */
      IntVector<NDIM> max_ghosts(1);
      CoarseFineBoundary<NDIM> cfbdry(*(Pointer<PatchHierarchy<NDIM> >)hierarchy, 
                                      level_number, 
                                      max_ghosts);

      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
      for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(p());
         Box<NDIM> pbox = patch->getBox();
         int pn = patch->getPatchNumber();

         const tbox::Array<BoundaryBox<NDIM> > node_bdry = 
            cfbdry.getNodeBoundaries(pn);
         const tbox::Array<BoundaryBox<NDIM> > edge_bdry = 
            cfbdry.getEdgeBoundaries(pn);
#if (NDIM == 3)
         const tbox::Array<BoundaryBox<NDIM> > face_bdry = 
            cfbdry.getFaceBoundaries(pn);
#else
         const tbox::Array<BoundaryBox<NDIM> > face_bdry;
#endif

         setBoundaryConditions(*patch,
                               node_bdry,
                               edge_bdry,
                               face_bdry,
                               d_ucell_edge_id);                               

      } // loop over level patches
      
   } // if level_number > 0
   
}

/*
 * Perform operations necessary when grid changes for dynamic grid
 * calculations (to be added later...)
 */
void 
HierSumTest::resetHierarchyConfiguration(
   const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   (void) hierarchy;
   (void) coarsest_level;
   (void) finest_level;
}

/*
 * Routine to do cell tagging.  Add later...
 */
void
HierSumTest::applyGradientDetector(
   const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   (void) hierarchy;
   (void) level_number;
   (void) time;
   (void) tag_index;
   (void) initial_time;
   (void) uses_richardson_extrapolation_too;
}

/*
 * Set boundary conditions by shifting patch appropriately and 
 * finding intersection with boundary fill box.
 */
void 
HierSumTest::setBoundaryConditions(
   Patch<NDIM>& patch,
   const tbox::Array<BoundaryBox<NDIM> >& node_bdry,
   const tbox::Array<BoundaryBox<NDIM> >& edge_bdry,
   const tbox::Array<BoundaryBox<NDIM> >& face_bdry,
   const int cell_data_id)
{

   const int num_node_bdry_boxes = node_bdry.getSize();
   const int num_edge_bdry_boxes = edge_bdry.getSize();
   const int num_face_bdry_boxes = face_bdry.getSize();

   const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = 
      patch.getPatchGeometry();

   /*
    * Pointer to data in ghost regions.
    */
   Pointer< CellData<NDIM,double> > ucell = patch.getPatchData(cell_data_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!ucell.isNull());
#endif
   IntVector<NDIM> ghost_cells = ucell->getGhostCellWidth();
   const Box<NDIM> pbox(patch.getBox());

   int i,d;
   IntVector<NDIM> shift(0);
   Box<NDIM> shifted_pbox = pbox;
   

#if (NDIM ==3) 
   /*
    * Set cell weights to 1.0 on FACES of patch.
    */
   for ( i = 0; i < num_face_bdry_boxes; i++ ) {
      Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(face_bdry[i],
                                                    pbox,
                                                    ghost_cells);
      /*
       * location index:
       *    0,1 - X lower,upper
       *    2,3 - Y lower,upper
       *    4,5 - Z lower,upper
       */   
      int loc_indx = face_bdry[i].getLocationIndex();
      shift = IntVector<NDIM>(0);
      shifted_pbox = pbox;
      if (loc_indx == 0) {
         shift(0) = -1;
      } else if (loc_indx == 1) {
         shift(0) = 1;
      } else if (loc_indx == 2) {
         shift(1) = -1;
      } else if (loc_indx == 3) {
         shift(1) = 1;
      } else if (loc_indx == 4) {
         shift(2) = -1;
      } else if (loc_indx == 5) {
         shift(2) = 1;
      }
      shifted_pbox.shift(shift);
      fill_box = fill_box * shifted_pbox;
         
      for (CellIterator<NDIM> i(fill_box); i; i++) {
         CellIndex<NDIM> cell = i();  
         for (d = 0; d < ucell->getDepth(); d++) {
            (*ucell)(cell,d) = 1.0;
         }
      }
   }
#endif   

   /*
    * Set cell weights to 1.0 on EDGES of patch.
    */
   for ( i = 0; i < num_edge_bdry_boxes; i++ ) {
      Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(edge_bdry[i],
                                                    pbox,
                                                    ghost_cells);


      int loc_indx = edge_bdry[i].getLocationIndex();
      shift = IntVector<NDIM>(0);
      shifted_pbox = pbox;

#if (NDIM == 3)
      /*
       * location index:
       * 3D:
       *    0,1,2,3 - ZloYlo, ZloYhi, ZhiYlo, ZhiYhi
       *    4,5,6,7 - XloZlo, XloZhi, XhiZlo, XhiZhi
       *    8,9,10,11 - YloXlo, YloXhi, YhiXlo, YhiXhi
       */
      if (loc_indx == 0) {
         shift(2) = -1;
         shift(1) = -1;
      } else if (loc_indx == 1) {
         shift(2) = -1;
         shift(1) = 1;
      } else if (loc_indx == 2) {
         shift(2) = 1;
         shift(1) = -1;
      } else if (loc_indx == 3) {
         shift(2) = 1;
         shift(1) = 1;
      } else if (loc_indx == 4) {
         shift(0) = -1;
         shift(2) = -1;
      } else if (loc_indx == 5) {
         shift(0) = -1;
         shift(2) = 1;
      } else if (loc_indx == 6) {
         shift(0) = 1;
         shift(2) = -1;
      } else if (loc_indx == 7) {
         shift(0) = 1;
         shift(2) = 1;
      } else if (loc_indx == 8) {
         shift(0) = -1;
         shift(1) = -1;
      } else if (loc_indx == 9) {
         shift(0) = 1;
         shift(1) = -1;
      } else if (loc_indx == 10) {
         shift(0) = -1;
         shift(1) = 1;
      } else if (loc_indx == 11) {
         shift(0) = 1;
         shift(1) = 1;
      }
#endif
#if (NDIM == 2)
      // to avoid warnings of unused vars
      (void) num_face_bdry_boxes;

      /*
       * location index:
       * 2D:
       *    0,1 - Xlo, Xhi
       *    2,3 - Ylo, Yhi
       */   
     if (loc_indx == 0) {
         shift(0) = -1;
      } else if (loc_indx == 1) {
         shift(0) = 1;
      } else if (loc_indx == 2) {
         shift(1) = -1;
      } else if (loc_indx == 3) {
         shift(1) = 1;
      }
#endif
      shifted_pbox.shift(shift);
      fill_box = fill_box * shifted_pbox;
         
      for (CellIterator<NDIM> i(fill_box); i; i++) {
         CellIndex<NDIM> cell = i();  
         for (d = 0; d < ucell->getDepth(); d++) {
            (*ucell)(cell,d) = 1.0;
         }
      }
   }

   /*
    * Set cell weights to 1.0 on NODES of patch.
    */
   
   for ( i = 0; i < num_node_bdry_boxes; i++ ) {
      Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(node_bdry[i],
                                                    pbox,
                                                    ghost_cells);
      for (CellIterator<NDIM> i(fill_box); i; i++) {
         CellIndex<NDIM> cell = i();  //i,j

         for (d = 0; d < ucell->getDepth(); d++) {
            (*ucell)(cell,d) = 1.0;
         }
      }
   }

}

/*
 * Zero out the cells on the physical boundary that lie at the
 * coarse-fine boundary.  This is needed for cases in which the
 * fine patch intersects the physical boundary.
 */
void HierSumTest::zeroOutPhysicalBoundaryCellsAtCoarseFineBoundary(
   Patch<NDIM>& cpatch,
   const int cell_data_id)
{

   const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = cpatch.getPatchGeometry();

   /*
    * Get node and edge boundary boxes.
    */

   const tbox::Array<BoundaryBox<NDIM> > edge_bdry =
      patch_geom->getCodimensionBoundary(NDIM-1);
   const int num_edge_bdry_boxes = edge_bdry.getSize();

#if (NDIM == 3)
   const tbox::Array<BoundaryBox<NDIM> > face_bdry =
      patch_geom->getCodimensionBoundary(1);
   const int num_face_bdry_boxes = face_bdry.getSize();
#endif

   /*
    * Pointer to data in ghost regions.
    */
   Pointer< CellData<NDIM,double> > ucell = cpatch.getPatchData(cell_data_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!ucell.isNull());
#endif
   IntVector<NDIM> ghost_cells = ucell->getGhostCellWidth();

   int i,d;
   Box<NDIM> cpbox = cpatch.getBox();
   double interior_value;
   

#if (NDIM ==3) 
   /*
    * Zero out values on the FACE boundary when they border a value of
    * zero on the interior.
    */
   for ( i = 0; i < num_face_bdry_boxes; i++ ) {
      Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(face_bdry[i],
                                                    cpbox,
                                                    ghost_cells);
      /*
       * location index:
       *    0,1 - X lower,upper
       *    2,3 - Y lower,upper
       *    4,5 - Z lower,upper
       */
      int loc_indx = face_bdry[i].getLocationIndex();
      for (CellIterator<NDIM> i(fill_box); i; i++) {
         CellIndex<NDIM> boundary_cell = i();  
         CellIndex<NDIM> interior = boundary_cell;
         if (loc_indx == 0) {
            interior(0) += 1;
         } else if (loc_indx == 1) {
            interior(0) -= 1;
         } else if (loc_indx == 2) {
            interior(1) += 1;
         } else if (loc_indx == 3) {
            interior(1) -= 1;
         } else if (loc_indx == 4) {
            interior(2) += 1;
         } else if (loc_indx == 5) {
            interior(2) -= 1;
         }

         for (d = 0; d < ucell->getDepth(); d++) {
            interior_value = (*ucell)(interior,d);
         
            if (tbox::MathUtilities<double>::equalEps(interior_value, 0.0)) {
               (*ucell)(boundary_cell,d) = 0.0;
            }
         }
      }
   }
#endif   
   /*
    * Zero out values on the EDGE boundary when they border a value of
    * zero on the interior.
    */
   for ( i = 0; i < num_edge_bdry_boxes; i++ ) {
      Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(edge_bdry[i],
                                                    cpbox,
                                                    ghost_cells);

      /*
       * location index:
       * 3D:
       *    0,1,2,3 - ZloYlo, ZloYhi, ZhiYlo, ZhiYhi
       *    4,5,6,7 - XloZlo, XloZhi, XhiZlo, XhiZhi
       *    8,9,10,11 - YloXlo, YloXhi, YhiXlo, YhiXhi
       * 2D:
       *    0,1 - Xlo, Xhi
       *    2,3 - Ylo, Yhi
       */
      int loc_indx = edge_bdry[i].getLocationIndex();
      for (CellIterator<NDIM> i(fill_box); i; i++) {
         CellIndex<NDIM> boundary_cell = i();  
         CellIndex<NDIM> interior = boundary_cell;
#if (NDIM == 3)
         if (loc_indx == 0) {
            interior(2) += 1;
            interior(1) += 1;
         } else if (loc_indx == 1) {
            interior(2) += 1;
            interior(1) -= 1;
         } else if (loc_indx == 2) {
            interior(2) -= 1;
            interior(1) += 1;
         } else if (loc_indx == 3) {
            interior(2) -= 1;
            interior(1) -= 1;
         } else if (loc_indx == 4) {
            interior(0) += 1;
            interior(2) += 1;
         } else if (loc_indx == 5) {
            interior(0) += 1;
            interior(2) -= 1;
         } else if (loc_indx == 6) {
            interior(0) -= 1;
            interior(2) += 1;
         } else if (loc_indx == 7) {
            interior(0) -= 1;
            interior(2) -= 1;
         } else if (loc_indx == 8) {
            interior(0) += 1;
            interior(1) += 1;
         } else if (loc_indx == 9) {
            interior(0) -= 1;
            interior(1) += 1;
         } else if (loc_indx == 10) {
            interior(0) += 1;
            interior(1) -= 1;
         } else if (loc_indx == 11) {
            interior(0) -= 1;
            interior(1) -= 1;
         }
#endif
#if (NDIM == 2)
         if (loc_indx == 0) {
            interior(0) += 1;
         } else if (loc_indx == 1) {
            interior(0) -= 1;
         } else if (loc_indx == 2) {
            interior(1) += 1;
         } else if (loc_indx == 3) {
            interior(1) -= 1;
         }
#endif
         for (d = 0; d < ucell->getDepth(); d++) {
            interior_value = (*ucell)(interior,d);
            
            if (tbox::MathUtilities<double>::equalEps(interior_value, 0.0)) {
               (*ucell)(boundary_cell,d) = 0.0;
            }
         }
      }
   }

   /*
    * Since we never use cell values at corners (nodes), no need to set 
    * them to anything.
    */
}




/*
*************************************************************************
*
* Get data from input database.                                         *
*                                                                       *
*************************************************************************
*/ 
void 
HierSumTest::getFromInput(Pointer<Database> input_db)
{

   int i;

   /*
    * Set number of ghosts for node and edge data.
    */ 
   tbox::Array<int> tmp_array;
   if (input_db->keyExists("node_ghosts")) {
      tmp_array = input_db->getIntegerArray("node_ghosts");
      if (tmp_array.getSize() != NDIM) {
         TBOX_ERROR("HierSumTest::getFromInput()" 
                    << "invalid 'node_ghosts' entry - must be integer"
                     << "array of size NDIM" << endl);
       }
   } else {
      tmp_array.resizeArray(NDIM);
      for (i = 0; i < NDIM; i++) {
         tmp_array[i] = 0;
      }
   }
   
   for (i = 0; i < NDIM; i++) {
      d_node_ghosts(i) = tmp_array[i];
   }
   
   if (input_db->keyExists("edge_ghosts")) {
      tmp_array = input_db->getIntegerArray("edge_ghosts");
      if (tmp_array.getSize() != NDIM) {
         TBOX_ERROR("HierSumTest::getFromInput()" 
                    << "invalid 'edge_ghosts' entry - must be integer"
                     << "array of size NDIM" << endl);
       }
   } else {
      tmp_array.resizeArray(NDIM);
      for (i = 0; i < NDIM; i++) {
         tmp_array[i] = 0;
      }
   }
   
   for (i = 0; i < NDIM; i++) {
      d_edge_ghosts(i) = tmp_array[i];
   }

   if (input_db->keyExists("var_depth")) {
      d_depth = input_db->getInteger("var_depth");
   }
   
   /*
    * See if we want to check data before communication
    */
   d_check_data_before_communication = false;
   if (input_db->keyExists("check_data_before_communication")) {
      d_check_data_before_communication = 
         input_db->getBool("check_data_before_communication");
   }
   
   
}


