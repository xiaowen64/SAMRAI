//
// File:        Hierarchy Sum test
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1973 $
// Modified:    $LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
// Description: SAMRAI interface class for hierarchy node and edge sum test
//

#include "SAMRAI_config.h"

/*
 * Header file for base classes.
 */
#include "StandardTagAndInitStrategy.h"

/*
 * Header file for SAMRAI classes referenced in this class.
 */
#include "Box.h"
#include "BoundaryBox.h"
#include "tbox/Database.h"
#include "CellVariable.h"
#include "IntVector.h"
#include "EdgeVariable.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchBoundaryNodeSum.h"
#include "PatchBoundaryEdgeSum.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "VariableContext.h"
#include "VisItDataWriter.h"

using namespace std;
using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace mesh;
using namespace pdat;
using namespace algs;
using namespace appu;


#ifndef NULL
#define NULL (0)
#endif

/*!
 * This code tests the hierarchy sum operations in SAMRAI.  Here's a brief
 * summary of what it does:
 *
 *   1. Create a cell-centered variable (ucell) which has ghost cell width
 *      equal 1.  Crate a node-centered variable (unode) which has zero
 *      ghosts. 
 *        - Set ucell = 1.0 on cells of patch INTERIORS
 *        - Set ucell = 0.0 on cells of patch GHOSTS
 *   2. Set ucell = 0.0 on cells of L < LN that are
 *      covered by refined cells.
 *   3. Set node values unode = sum(surrounding cells)
 *   3. Do a hier sum transaction
 *   4. Correct result - all nodes on all levels = 2^NDIM 
 * 
 * The class provides interfaces to problem dependent operations that
 * are expected by gridding operations performed in SAMRAI.  For example, 
 * initializing data on a level and resetting hierarchy after regrid.
 *
 * Other methods may be added as needed.
 * 
 * Input Parameters:
 * 
 *
 * A sample input entry might look like:
 *  
 *    HierSumTest {
 *
 *    }
 */

class HierSumTest : 
   public StandardTagAndInitStrategy<NDIM>
{
public:
   /*!
    * Default constructor for HierSumTest.
    */     
   HierSumTest(const string& object_name,
               Pointer<Database> input_db,
               Pointer<VisItDataWriter<NDIM> > viz_writer);

   /*!
    * Empty destructor for HierSumTest.
    */
   virtual ~HierSumTest();

/*************************************************************************
 *
 * Methods particular to HierSumTest class.
 *
 ************************************************************************/

   /*!
    * Set node values before the hierarchy sum operation and return integer 
    * number of failures.
    */
   int setInitialNodeValues(const Pointer<PatchHierarchy<NDIM> > hierarchy);

   /*!
    * Set edge values before the level sum operation and return integer
    * number of failures.
    */
   int setInitialEdgeValues(const Pointer<PatchLevel<NDIM> > level);

   /*!
    * Setup the node hierarchy sum.
    */
   void setupOuternodeSum(const Pointer<PatchHierarchy<NDIM> > hierarchy);

   /*!
    * Invoke the node hierarchy sum communication.
    */
   void doOuternodeSum();

   /*!
    * Setup the edge level sum.
    */
   void setupOuteredgeSum(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                          const int level_num);
   
   /*!
    * Invoke the edge level sum communication.
    */
   void doOuteredgeSum(const int level_num);
   
   /*!
    * Check node result after hierarchy sum operation and return integer number of
    * test failures.
    */
   int checkNodeResult(const Pointer<PatchHierarchy<NDIM> > hierarchy);

   /*!
    * Check edge result after level sum operation and return integer number of
    * test failures.
    */
   int checkEdgeResult(const Pointer<PatchLevel<NDIM> > level);



/***************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy<NDIM>.
 *
 ************************************************************************/
 
   /*!
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.
    *
    * Generally, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the old level pointer in the argument list is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.   In this case, the level number must match that of the old 
    * level.  The specific operations that occur when initializing level 
    * data are determined by the particular solution methods in use; i.e.,
    * in the subclass of this abstract base class.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    */
   virtual void
   initializeLevelData(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                       const int level_number,
                       const double time,
                       const bool can_be_refined,
                       const bool initial_time,
                       const Pointer<BasePatchLevel<NDIM> > old_level = 
                             Pointer<BasePatchLevel<NDIM> >(NULL),
                       const bool allocate_data = true);

   /*!
    * After hierarchy levels have changed and data has been initialized on 
    * the new levels, this routine can be used to reset any information 
    * needed by the solution method that is particular to the hierarchy 
    * configuration.  For example, the solution procedure may cache 
    * communication schedules to amortize the cost of data movement on the 
    * AMR patch hierarchy.  This function will be called by the gridding 
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must 
    * make the solution consistent across multiple levels after the hierarchy 
    * is changed, this process may be invoked by this routine.  Of course the 
    * details of these processes are determined by the particular solution 
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void resetHierarchyConfiguration(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    * Set tags to the specified tag value where refinement of the given
    * level should occur using the user-supplied gradient detector.  The 
    * value "tag_index" is the index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument indicates 
    * whether cells are being tagged on the level for the first time; 
    * i.e., when the hierarchy is initially constructed.  If it is false, 
    * it should be assumed that cells are being tagged at some later time 
    * after the patch hierarchy was initially constructed.  This information
    * is provided since the application of the error estimator may be 
    * different in each of those circumstances.
    */
   virtual void 
   applyGradientDetector(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                         const int level_number,
                         const double time,
                         const int tag_index,
                         const bool initial_time,
                         const bool uses_richardson_extrapolation_too);


private:
   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    *
    * An assertion results if the database pointer is null.
    */
   virtual void getFromInput(tbox::Pointer<tbox::Database> db);

   /*
    * Set boundary conditions at physical boundaries and coarse-fine 
    * boundaries.  Note that this is not called by the refine
    * schedule because we do not use a refine schedule in this test.  
    * It is called during initializeLevel().
    */
   void setBoundaryConditions(
      Patch<NDIM>& patch,
      const tbox::Array<BoundaryBox<NDIM> >& node_bdry,
      const tbox::Array<BoundaryBox<NDIM> >& edge_bdry,
      const tbox::Array<BoundaryBox<NDIM> >& face_bdry,
      const int cell_data_id);

   /*
    * Zero out the coarse cell boundary values where fine patch intersects
    * the physical boundary.
    */
   void zeroOutPhysicalBoundaryCellsAtCoarseFineBoundary(
      Patch<NDIM>& cpatch,
      const int cell_data_id);   


   
   /*
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   string d_object_name;

   /*
    * Variable<NDIM> - u
    */
   Pointer< CellVariable<NDIM,double> > d_ucell_var;
   Pointer< NodeVariable<NDIM,double> > d_unode_var;
   Pointer< EdgeVariable<NDIM,double> > d_uedge_var;

   /*
    * Ghost vectors
    */
   IntVector<NDIM> d_node_ghosts;
   IntVector<NDIM> d_edge_ghosts;   

   /*
    * Patch<NDIM> Data ids - used to access data off patch
    */
   int d_ucell_node_id;
   int d_ucell_edge_id;
   int d_unode_id;
   int d_uedge_id;
   
   /*
    * Node, edge, and cell variable depths - all are equal.
    */
   int d_depth;

   /*
    * Node and edge sum utilities.
    */
   Pointer<PatchBoundaryNodeSum<NDIM> > d_node_sum_util;
   Array<Pointer<PatchBoundaryEdgeSum<NDIM> > > d_edge_sum_util;
   
   /*
    * Flag to tell whether to check data before communication.  Usually,
    * this will be false, but sometimes you may want to see if there are
    * actually errors in the data before the communication.
    */
   bool d_check_data_before_communication;
   
   

};

