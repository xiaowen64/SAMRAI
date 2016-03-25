//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/algorithm/femutils/standard/PatchBoundaryNodeSum.C $
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description:	Routines for summing node data at patch boundaries
//

#ifndef included_algs_PatchBoundaryNodeSum_C
#define included_algs_PatchBoundaryNodeSum_C

#include "PatchBoundaryNodeSum.h"

#include "VariableDatabase.h"
#include "NodeData.h"
#include "NodeDataFactory.h"
#include "NodeGeometry.h"
#include "OuternodeData.h"
#include "OuternodeDoubleConstantCoarsen.h"
#include "OuternodeSumTransactionFactory.h"
#include "CoarsenAlgorithm.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"
#include "RefineOperator.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN 77 routines used to sum node and    *
* outernode data.                                                       *
*                                                                       *
*************************************************************************
*/

extern "C" {
// in algs_nodeouternodeops2d.f:
   void nodeouternodesum2d_( 
      const int&, const int&,  // fine patch lo
      const int&, const int&,  // fine patch hi
      const int&, const int&,  // coarse patch lo
      const int&, const int&,  // coarse patch hi
      const int*,              // ratio vector
      const int&,              // data depth
      const int&, const int&,  // node data gcw
      double*,                 // node array
      const double*, const double*,  // onode arrays
      const double*, const double*);

   void nodehangnodeinterp2d_( 
      const int&, const int&, // fine patch lo 
      const int&, const int&, // fine patch hi
      const int&, const int&, // coarse patch lo
      const int&, const int&, // coarse patch hi
      const int&, const int&, // bbox lo
      const int&, const int&, // bbox hi
      const int&,             // bbox location
      const int*,             // ratio vector
      const int&,             // data depth
      const int&, const int&, // node data gcw
      double*);               // node array

// in algs_nodeouternodeops3d.f:
   void nodeouternodesum3d_( 
      const int&, const int&, const int&, // fine patch lo
      const int&, const int&, const int&, // fine patch hi
      const int&, const int&, const int&, // coarse patch lo
      const int&, const int&, const int&, // coarse patch hi
      const int*,                         // ratio vector
      const int&,                         // data depth
      const int&, const int&, const int&, // node data gcw
      double*,                            // node array
      const double*, const double*,       // onode arrays
      const double*, const double*,
      const double*, const double*);

   void nodehangnodeinterp3d_( 
      const int&, const int&, const int&, // fine patch lo
      const int&, const int&, const int&, // fine patch hi
      const int&, const int&, const int&, // coarse patch lo
      const int&, const int&, const int&, // coarse patch hi
      const int&, const int&, const int&, // bbox lo
      const int&, const int&, const int&, // bbox hi
      const int&,                         // bbox location
      const int*,                         // ratio vector
      const int&,                         // data depth
      const int&, const int&, const int&, // node data gcw
      double*);                           // node array
}

namespace SAMRAI {
   namespace algs {

#ifndef NULL
#define NULL (0)
#endif

/*
*************************************************************************
*                                                                       *
* Initialize the static data members.                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> int PatchBoundaryNodeSum<DIM>::s_instance_counter = 0;

template<int DIM> tbox::Array< tbox::Array<int> > 
   PatchBoundaryNodeSum<DIM>::s_onode_src_id_array =
      tbox::Array< tbox::Array<int> >(0);
template<int DIM> tbox::Array< tbox::Array<int> > 
   PatchBoundaryNodeSum<DIM>::s_onode_dst_id_array =
      tbox::Array< tbox::Array<int> >(0);

/*
*************************************************************************
*                                                                       *
* Static functions to determine number of patch data slots needed       *
* for PatchBoundaryNodeSum objects.                                     *
*                                                                       *
*************************************************************************
*/
 
template<int DIM>
int
PatchBoundaryNodeSum<DIM>::getNumSharedPatchDataSlots(
   int max_variables_to_register)
{
   // node boundary sum requires two internal outernode variables
   // (source and destination) for each registered variable.
 
   return( 2 * max_variables_to_register );
}
 
template<int DIM>
int
PatchBoundaryNodeSum<DIM>::getNumUniquePatchDataSlots(
   int max_variables_to_register)
{
   NULL_USE(max_variables_to_register);
   // all patch data slots used by node boundary sum are static
   // and shared among all objects.
 
   return( 0 );
}

/*
*************************************************************************
*                                                                       *
* Constructor patch boundary node sum objects initializes data members  *
* to default (undefined) states.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> PatchBoundaryNodeSum<DIM>::PatchBoundaryNodeSum(
   const std::string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
#endif

   d_object_name = object_name;
   d_setup_called = false;

   d_num_reg_sum = 0;

   d_level.setNull();
   d_hierarchy.setNull();
   d_coarsest_level = -1;
   d_finest_level = -1;

   d_level_setup_called     = false;
   d_hierarchy_setup_called = false;

   d_sum_transaction_factory = new OuternodeSumTransactionFactory<DIM>();

   s_instance_counter++;   
}

/*
*************************************************************************
*                                                                       *
* Destructor removes temporary outernode patch data ids from            *
* variable database, if defined.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> PatchBoundaryNodeSum<DIM>::~PatchBoundaryNodeSum()
{

   s_instance_counter--;
   if (s_instance_counter == 0) {
      const int arr_length_depth = s_onode_src_id_array.size();

      for (int id = 0; id < arr_length_depth; id++) {
         const int arr_length_nvar = s_onode_src_id_array[id].size(); 

         for (int iv = 0; iv < arr_length_nvar; iv++) {
          
            if (s_onode_src_id_array[id][iv] >= 0) {
               hier::VariableDatabase<DIM>::getDatabase()->
                  removeInternalSAMRAIVariablePatchDataIndex(
                     s_onode_src_id_array[id][iv]);
            }
            if (s_onode_dst_id_array[id][iv] >= 0) {
               hier::VariableDatabase<DIM>::getDatabase()->
                  removeInternalSAMRAIVariablePatchDataIndex(
                     s_onode_dst_id_array[id][iv]);
            }

            s_onode_src_id_array[id].resizeArray(0); 
            s_onode_dst_id_array[id].resizeArray(0); 

         }

      }

      s_onode_src_id_array.resizeArray(0);
      s_onode_dst_id_array.resizeArray(0);

   }

}

/*
*************************************************************************
*                                                                       *
* Register node patch data index for summation.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::registerSum(
   int node_data_id)
{

   if (d_setup_called) {
      TBOX_ERROR("PatchBoundaryNodeSum<DIM>::register error..."
                 << "\nobject named " << d_object_name
                 << "\nCannot call registerSum with this PatchBoundaryNodeSum"
                 << "\nobject since it has already been used to create communication"
                 << "\nschedules; i.e., setupSum() has been called."
                 << std::endl);
   }

   if (node_data_id < 0) {
      TBOX_ERROR("PatchBoundaryNodeSum<DIM> register error..."
                 << "\nobject named " << d_object_name
                 << "\n node_data_id = " << node_data_id
                 << " is an invalid patch data identifier." << std::endl);
   }

   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   tbox::Pointer< pdat::NodeDataFactory<DIM,double> > node_factory =
      var_db->getPatchDescriptor()->getPatchDataFactory(node_data_id);

   if (node_factory.isNull()) {

      TBOX_ERROR("PatchBoundaryNodeSum<DIM> register error..."
                 << "\nobject named " << d_object_name
                 << "\n node_data_id = " << node_data_id
                 << " does not correspond to node data of type double." << std::endl);

   } else {

      static std::string tmp_onode_src_variable_name("PatchBoundaryNodeSum__internal-onode-src");
      static std::string tmp_onode_dst_variable_name("PatchBoundaryNodeSum__internal-onode-dst");

      const int reg_sum_id = d_num_reg_sum;

      d_num_reg_sum++;

      d_user_node_data_id.resizeArray(d_num_reg_sum);
         d_user_node_data_id[reg_sum_id] = ID_UNDEFINED;
      d_user_node_depth.resizeArray(d_num_reg_sum);
         d_user_node_depth[reg_sum_id] = ID_UNDEFINED;
      d_tmp_onode_src_variable.resizeArray(d_num_reg_sum);
      d_tmp_onode_dst_variable.resizeArray(d_num_reg_sum);
      d_onode_src_id.resizeArray(d_num_reg_sum);
         d_onode_src_id[reg_sum_id] = ID_UNDEFINED;
      d_onode_dst_id.resizeArray(d_num_reg_sum);
         d_onode_dst_id[reg_sum_id] = ID_UNDEFINED;

      const int data_depth = node_factory->getDefaultDepth();
      const int array_by_depth_size = data_depth + 1;

      if (d_num_registered_data_by_depth.size() < array_by_depth_size) {
         const int old_size = d_num_registered_data_by_depth.size();
         const int new_size = array_by_depth_size;
       
         d_num_registered_data_by_depth.resizeArray(new_size); 
         for (int i = old_size; i < new_size; i++) {
            d_num_registered_data_by_depth[i] = 0;
         } 
      }

      const int data_depth_id = d_num_registered_data_by_depth[data_depth];
      const int num_data_at_depth = data_depth_id + 1;

      if (s_onode_src_id_array.size() < array_by_depth_size) {
         s_onode_src_id_array.resizeArray(array_by_depth_size);
         s_onode_dst_id_array.resizeArray(array_by_depth_size);
      }

      if (s_onode_src_id_array[data_depth].size() < num_data_at_depth) {
         const int old_size = s_onode_src_id_array[data_depth].size();
         const int new_size = num_data_at_depth;

         s_onode_src_id_array[data_depth].resizeArray(new_size);
         s_onode_dst_id_array[data_depth].resizeArray(new_size);
         for (int i = old_size; i < new_size; i++) {
            s_onode_src_id_array[data_depth][i] = ID_UNDEFINED;
            s_onode_dst_id_array[data_depth][i] = ID_UNDEFINED;
         }
      }

      char var_suffix[17];
      sprintf(var_suffix, "%04d__depth=%04d", data_depth_id, data_depth);

      std::string tonode_src_var_name = tmp_onode_src_variable_name + var_suffix;
      d_tmp_onode_src_variable[reg_sum_id] = var_db->getVariable(tonode_src_var_name);
      if (d_tmp_onode_src_variable[reg_sum_id].isNull()) {
         d_tmp_onode_src_variable[reg_sum_id] =
            new pdat::OuternodeVariable<DIM,double>(tonode_src_var_name, data_depth);
      }

      std::string tonode_dst_var_name = tmp_onode_dst_variable_name + var_suffix;
      d_tmp_onode_dst_variable[reg_sum_id] = var_db->getVariable(tonode_dst_var_name);
      if (d_tmp_onode_dst_variable[reg_sum_id].isNull()) {
         d_tmp_onode_dst_variable[reg_sum_id] =
            new pdat::OuternodeVariable<DIM,double>(tonode_dst_var_name, data_depth);
      }

      if ( s_onode_src_id_array[data_depth][data_depth_id] < 0 ) {
         s_onode_src_id_array[data_depth][data_depth_id] =
            var_db->registerInternalSAMRAIVariable(
                    d_tmp_onode_src_variable[reg_sum_id], 
                    hier::IntVector<DIM>(0));
      }
      if ( s_onode_dst_id_array[data_depth][data_depth_id] < 0) {
         s_onode_dst_id_array[data_depth][data_depth_id] =
            var_db->registerInternalSAMRAIVariable(
                    d_tmp_onode_dst_variable[reg_sum_id], 
                    hier::IntVector<DIM>(0));
      }

      d_user_node_data_id[reg_sum_id] = node_data_id;
      d_user_node_depth[reg_sum_id] = data_depth;

      d_num_registered_data_by_depth[data_depth] = num_data_at_depth;

      d_onode_src_id[reg_sum_id] = s_onode_src_id_array[data_depth][data_depth_id];
      d_onode_dst_id[reg_sum_id] = s_onode_dst_id_array[data_depth][data_depth_id];

      d_onode_src_data_set.setFlag(d_onode_src_id[reg_sum_id]);
      d_onode_dst_data_set.setFlag(d_onode_dst_id[reg_sum_id]);

   }

}

/*
*************************************************************************
*                                                                       *
* Set up schedule to sum node data around patch boundaries              *
* on a single level.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::setupSum(
   tbox::Pointer<hier::PatchLevel<DIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   if (d_hierarchy_setup_called) {

      TBOX_ERROR("PatchBoundaryNodeSum<DIM>::setupSum error...\n"
         << " object named " << d_object_name
         << " already initialized via setupSum() for hierarchy" << std::endl);

   } else {

      d_setup_called = true;
      d_level_setup_called = true;

      d_level = level;

      d_single_level_sum_schedule.resizeArray(1);

      // Communication algorithm for summing outernode values on a level
      xfer::RefineAlgorithm<DIM> single_level_sum_algorithm;

      for (int i = 0; i < d_num_reg_sum; i++) {
         single_level_sum_algorithm.registerRefine(d_onode_dst_id[i],  // dst data
                                                   d_onode_src_id[i],  // src data
                                                   d_onode_dst_id[i],  // scratch data
                                                   (xfer::RefineOperator<DIM>*)NULL);
      }

      d_single_level_sum_schedule[0] =
         single_level_sum_algorithm.createSchedule(
                                    d_level,
                                    (xfer::RefinePatchStrategy<DIM>*)NULL,
                                    d_sum_transaction_factory);

   } 

}

/*
*************************************************************************
*                                                                       *
* Set up schedule to sum node data around patch boundaries              *
* for set of consecutive hierarchy levels.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::setupSum(
   tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT(   (coarsest_level >= 0)
          && (finest_level >= coarsest_level)
          && (finest_level <= hierarchy->getFinestLevelNumber()) );
#endif

   if (d_level_setup_called) {

      TBOX_ERROR("PatchBoundaryNodeSum<DIM>::setupSum error...\n"
         << " object named " << d_object_name
         << " already initialized via setupSum() for single level" << std::endl);

   } else {

      d_setup_called = true; 
      d_hierarchy_setup_called = true;

      d_hierarchy = hierarchy;
      d_coarsest_level = coarsest_level;
      d_finest_level = finest_level;

      const int num_levels = d_finest_level+1;

      d_single_level_sum_schedule.resizeArray(num_levels);
      d_cfbdry_copy_schedule.resizeArray(num_levels);
      d_sync_coarsen_schedule.resizeArray(num_levels);
      d_cfbdry_tmp_level.resizeArray(num_levels);

      d_coarse_fine_boundary.resizeArray(num_levels);

      // Communication algorithm for summing outernode values on each level
      xfer::RefineAlgorithm<DIM> single_level_sum_algorithm;

      // Communication algorithm for copying node values on each coarser 
      // level to outernode values on coarsened version of patches on 
      // next finer level
      xfer::RefineAlgorithm<DIM> cfbdry_copy_algorithm;

      // Communication algorithm for coarsening outernode values on 
      // each finer level to node data on next coarser level
      xfer::CoarsenAlgorithm<DIM> sync_coarsen_algorithm(false);
      tbox::Pointer<pdat::OuternodeDoubleConstantCoarsen<DIM> > coarsen_op =
         new pdat::OuternodeDoubleConstantCoarsen<DIM>();

      for (int i = 0; i < d_num_reg_sum; i++) {
         single_level_sum_algorithm.registerRefine(d_onode_dst_id[i],  // dst data
                                                   d_onode_src_id[i],  // src data
                                                   d_onode_dst_id[i],  // scratch data
                                                   (xfer::RefineOperator<DIM>*)NULL);

         cfbdry_copy_algorithm.registerRefine(d_onode_dst_id[i],      // dst data
                                              d_user_node_data_id[i], // src data
                                              d_onode_dst_id[i],      // scratch data
                                              (xfer::RefineOperator<DIM>*)NULL);

         sync_coarsen_algorithm.registerCoarsen(d_user_node_data_id[i], // dst data
                                                d_onode_dst_id[i],      // src data
                                                coarsen_op);
      }

      tbox::Pointer< hier::PatchLevel<DIM> > coarsest_level_loop = 
         d_hierarchy->getPatchLevel(d_coarsest_level);

      tbox::Pointer<hier::PatchLevel<DIM> > hier_level = 
         d_hierarchy->getPatchLevel(0);

      d_single_level_sum_schedule[d_coarsest_level] =
         single_level_sum_algorithm.createSchedule(
                                    coarsest_level_loop,
                                    (xfer::RefinePatchStrategy<DIM>*)NULL,
                                    d_sum_transaction_factory);

      for (int ln = d_coarsest_level+1; ln <= d_finest_level; ln++) {

         const int crse_level_num = ln-1;
         const int fine_level_num = ln;

         tbox::Pointer< hier::PatchLevel<DIM> > crse_level = 
            d_hierarchy->getPatchLevel(crse_level_num);
         tbox::Pointer< hier::PatchLevel<DIM> > fine_level = 
            d_hierarchy->getPatchLevel(fine_level_num);

         d_single_level_sum_schedule[fine_level_num] =
            single_level_sum_algorithm.createSchedule(
                                       fine_level,
                                       (xfer::RefinePatchStrategy<DIM>*)NULL,
                                       d_sum_transaction_factory);

         d_cfbdry_tmp_level[fine_level_num] = new hier::PatchLevel<DIM>();
         d_cfbdry_tmp_level[fine_level_num]->
            setCoarsenedPatchLevel(fine_level,
                                   fine_level->getRatioToCoarserLevel()); 
  
         d_cfbdry_copy_schedule[fine_level_num] = 
            cfbdry_copy_algorithm.createSchedule(
               d_cfbdry_tmp_level[fine_level_num],     // dst level 
               crse_level);                            // src level

         d_sync_coarsen_schedule[fine_level_num] = 
            sync_coarsen_algorithm.createSchedule(   
               crse_level,                             // dst level
               fine_level);                            // src level

         d_coarse_fine_boundary[fine_level_num].
            computeFromLevel(*(d_cfbdry_tmp_level[fine_level_num]), 
                             *hier_level,
                             hier::IntVector<DIM>(1));

      }

   }

}

/*
*************************************************************************
*                                                                       *
* Perform patch boundary node sum across single level or multiple       *
* hierarchy levels depending on how object was initialized.  In the     *
* single level case, values at shared nodes are summed.  In the         *
* multiple-level case, the algorithm involves the following steps:      *
*                                                                       *
*    1) Sum values at shared nodes on each level.                       *
*    2) Set node values at coarse-fine boundary on each finer level     *
*       to sum of fine level values and coarse level values at all      *
*       nodes that are shared between the coarse and fine level.        *
*                                                                       *
*       2a) Copy coarser level node values to finer level (outer)node   *
*           values at nodes on boundary of patches on a temporary       *
*           level that represents the finer level coarsened to the      *
*           index space of the coarser level.                           *
*       2b) Sum (outer)node values at patch boundaries on finer level   *
*           and (outer)node values at patch boundaries on coarsened     *
*           finer level and set values on finer level to sum.  Note     *
*           that the hanging nodes on the finer level may be treated    *
*           at this point if specified to do so by the user.            *
*                                                                       *
*    3) Inject (outer)node values around each finer level patch         *
*       boundary to corresponding node values on each coarser level.    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::computeSum(
   const bool fill_hanging_nodes) const
{

   if (d_level_setup_called) {

      d_level->allocatePatchData(d_onode_src_data_set);
      d_level->allocatePatchData(d_onode_dst_data_set);

      doLevelSum(d_level);

      d_level->deallocatePatchData(d_onode_src_data_set);
      d_level->deallocatePatchData(d_onode_dst_data_set);

   } else {  // assume d_hierarchy_setup_called

      int ln;

      for (ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         
        tbox::Pointer< hier::PatchLevel<DIM> > level =
            d_hierarchy->getPatchLevel(ln);
   
        level->allocatePatchData(d_onode_src_data_set);
        level->allocatePatchData(d_onode_dst_data_set);

        doLevelSum(level);

      }

      for (ln = d_coarsest_level+1; ln <= d_finest_level; ln++) {

         tbox::Pointer< hier::PatchLevel<DIM> > level =
            d_hierarchy->getPatchLevel(ln);

         d_cfbdry_tmp_level[ln]->allocatePatchData(d_onode_dst_data_set);

         d_cfbdry_copy_schedule[ln]->fillData(0.0, false);

         doLocalCoarseFineBoundarySum(level,
                                      d_cfbdry_tmp_level[ln],
                                      d_user_node_data_id,
                                      d_onode_dst_id,
                                      fill_hanging_nodes);

         d_cfbdry_tmp_level[ln]->deallocatePatchData(d_onode_dst_data_set);

      }
   
      for (ln = d_finest_level; ln > d_coarsest_level; ln--) {

         tbox::Pointer< hier::PatchLevel<DIM> > level =
            d_hierarchy->getPatchLevel(ln);

         copyNodeToOuternodeOnLevel(level,
                                    d_user_node_data_id,
                                    d_onode_dst_id);

         d_sync_coarsen_schedule[ln]->coarsenData();

         level->deallocatePatchData(d_onode_src_data_set);
         level->deallocatePatchData(d_onode_dst_data_set);
          
      } 

      d_hierarchy->getPatchLevel(d_coarsest_level)->
         deallocatePatchData(d_onode_src_data_set);
      d_hierarchy->getPatchLevel(d_coarsest_level)->
         deallocatePatchData(d_onode_dst_data_set);

   }  // if d_hierarchy_setup_called

}

/*
*************************************************************************
*                                                                       *
* Private member function that performs node sum across single level.   *
*                                                                       *
* 1. Copy node data to local outernode data.                            *
* 2. Transfer and sum outernode data.                                   *
* 3. Copy outernode data back to node data.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::doLevelSum(
   tbox::Pointer<hier::PatchLevel<DIM> > level) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   copyNodeToOuternodeOnLevel(level,
                              d_user_node_data_id,
                              d_onode_src_id);

   int schedule_level_number = 0;
   if (!d_level_setup_called) {
      schedule_level_number = 
         tbox::MathUtilities<int>::Max(0, level->getLevelNumber()); 
   }
   d_single_level_sum_schedule[schedule_level_number]->fillData(0.0, false);

   copyOuternodeToNodeOnLevel(level,
                              d_onode_dst_id,
                              d_user_node_data_id);   
}

/*
*************************************************************************
*                                                                       *
* Private member function to set node node data on a fine level at a    *
* coarse-fine boundary to the sum of the node values and the associated *
* outernode values on a coarsened version of the fine level.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::doLocalCoarseFineBoundarySum(
   tbox::Pointer<hier::PatchLevel<DIM> > fine_level,
   tbox::Pointer<hier::PatchLevel<DIM> > coarsened_fine_level,
   const tbox::Array<int>& node_data_id,
   const tbox::Array<int>& onode_data_id,
   bool fill_hanging_nodes) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fine_level.isNull());
   TBOX_ASSERT(!coarsened_fine_level.isNull());
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());
   for (int i = 0; i < node_data_id.size(); i++) {
      TBOX_ASSERT(fine_level->checkAllocated(node_data_id[i]));
      TBOX_ASSERT(coarsened_fine_level->checkAllocated(onode_data_id[i]));
   }
#endif

   const hier::IntVector<DIM>& ratio = fine_level->getRatioToCoarserLevel();
   const int level_number = fine_level->getLevelNumber();

   for (typename hier::PatchLevel<DIM>::Iterator ip(fine_level); ip; ip++ ) {
      
      const tbox::Array<hier::BoundaryBox<DIM> >& pboundaries =
         d_coarse_fine_boundary[level_number].getBoundaries(ip(), 1);
      const int num_bdry_boxes = pboundaries.getSize();

      if (num_bdry_boxes > 0) {

         tbox::Pointer<hier::Patch<DIM> > fpatch = fine_level->getPatch(ip());
         tbox::Pointer<hier::Patch<DIM> > cfpatch = 
            coarsened_fine_level->getPatch(ip());

         const hier::Index<DIM>& filo = fpatch->getBox().lower();
         const hier::Index<DIM>& fihi = fpatch->getBox().upper();
         const hier::Index<DIM>& cilo = cfpatch->getBox().lower();
         const hier::Index<DIM>& cihi = cfpatch->getBox().upper();

         for (int i = 0; i < node_data_id.size(); i++) {

            tbox::Pointer< pdat::NodeData<DIM,double> > node_data =
               fpatch->getPatchData(node_data_id[i]);
            tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_data =
               cfpatch->getPatchData(onode_data_id[i]);

            const hier::IntVector<DIM>& node_gcw = node_data->getGhostCellWidth();

            pdat::OuternodeData<DIM,double> tmp_onode_data(
               cfpatch->getBox(), onode_data->getDepth());
            tmp_onode_data.fillAll(0.0);

            // Copy "coarse" node values on coarse-fine boundary to
            // temporary outernode data arrays.
            for (int ibb0 = 0; ibb0 < num_bdry_boxes; ibb0++) {

               const hier::BoundaryBox<DIM>& bbox = pboundaries[ibb0];
               const int bbox_loc = bbox.getLocationIndex();

               hier::Box<DIM> node_bbox = 
                  pdat::NodeGeometry<DIM>::toNodeBox(bbox.getBox());
            
               switch(bbox_loc) {
   
                  case 0 : {
                  
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_00 = 
                        tmp_onode_data.getArrayData(0,0);
                     if (tmp_onode_data_side_00.isInitialized()) {
                        tmp_onode_data_side_00.
                           copy(onode_data->getArrayData(0,0), node_bbox);
                     }
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_10 = 
                        tmp_onode_data.getArrayData(1,0);
                     if (tmp_onode_data_side_10.isInitialized()) {
                        tmp_onode_data_side_10.
                           copy(onode_data->getArrayData(1,0), node_bbox);
                     }
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_11 = 
                        tmp_onode_data.getArrayData(1,1);
                     if (tmp_onode_data_side_11.isInitialized()) {
                        tmp_onode_data_side_11.
                           copy(onode_data->getArrayData(1,1), node_bbox);
                     }

                     if (DIM == 3) {
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_20 = 
                           tmp_onode_data.getArrayData(2,0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                              copy(onode_data->getArrayData(2,0), node_bbox);
                        }
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_21 = 
                           tmp_onode_data.getArrayData(2,1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                              copy(onode_data->getArrayData(2,1), node_bbox);
                        }
                     }
                     break;
                  }  // case 0
    
                  case 1 : {
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_01 = 
                        tmp_onode_data.getArrayData(0,1);
                     if (tmp_onode_data_side_01.isInitialized()) {
                        tmp_onode_data_side_01.
                           copy(onode_data->getArrayData(0,1), node_bbox);
                     }
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_10 = 
                        tmp_onode_data.getArrayData(1,0);
                     if (tmp_onode_data_side_10.isInitialized()) {
                        tmp_onode_data_side_10.
                           copy(onode_data->getArrayData(1,0), node_bbox);
                     }
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_11 = 
                        tmp_onode_data.getArrayData(1,1);
                     if (tmp_onode_data_side_11.isInitialized()) {
                        tmp_onode_data_side_11.
                           copy(onode_data->getArrayData(1,1), node_bbox);
                     }

                     if (DIM == 3) {
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_20 = 
                           tmp_onode_data.getArrayData(2,0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                              copy(onode_data->getArrayData(2,0), node_bbox);
                        }
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_21 = 
                           tmp_onode_data.getArrayData(2,1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                              copy(onode_data->getArrayData(2,1), node_bbox);
                        }
                     }
                     break;
                  } // case 1

                  case 2 : {
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_10 = 
                        tmp_onode_data.getArrayData(1,0);
                     if (tmp_onode_data_side_10.isInitialized()) {
                        tmp_onode_data_side_10.
                           copy(onode_data->getArrayData(1,0), node_bbox);
                     }

                     if (DIM == 3) {
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_20 = 
                           tmp_onode_data.getArrayData(2,0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                              copy(onode_data->getArrayData(2,0), node_bbox);
                        }
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_21 = 
                           tmp_onode_data.getArrayData(2,1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                           copy(onode_data->getArrayData(2,1), node_bbox);
                        }
                     }
                     break;
                  }  // case 2

                  case 3 : {
                     pdat::ArrayData<DIM,double>& tmp_onode_data_side_11 = 
                     tmp_onode_data.getArrayData(1,1);
                     if (tmp_onode_data_side_11.isInitialized()) {
                        tmp_onode_data_side_11.
                           copy(onode_data->getArrayData(1,1), node_bbox);
                     }
   
                     if (DIM == 3) {
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_20 = 
                           tmp_onode_data.getArrayData(2,0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                              copy(onode_data->getArrayData(2,0), node_bbox);
                        }
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_21 = 
                           tmp_onode_data.getArrayData(2,1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                                 copy(onode_data->getArrayData(2,1), node_bbox);
                        }
                     }
                     break;
                  }  // case 3

                  case 4 : {
                     if (DIM == 3) {
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_20 = 
                           tmp_onode_data.getArrayData(2,0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                              copy(onode_data->getArrayData(2,0), node_bbox);
                        }
                     }
                     break;
                  }  // case 4

                  case 5 : {
                     if (DIM == 3) {
                        pdat::ArrayData<DIM,double>& tmp_onode_data_side_21 = 
                           tmp_onode_data.getArrayData(2,1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                              copy(onode_data->getArrayData(2,1), node_bbox);
                        }
                     }
                     break;
                  }  // case 5

                  default : {
                     TBOX_ERROR("PatchBoundaryNodeSum::computeSum error...\n"
                        << " object named " << d_object_name
                        << "\n unrecognized coarse-fine boundary box location " 
                        << bbox_loc << std::endl);
                  }

               }  // switch(box_loc)

            } // for (int ibb0 ...  iterate over coarse-fine boundary box regions

            // Sum "coarse" node values on coarse-fine boundary.

            if (DIM == 2) {
         
               if (tmp_onode_data.getArrayData(0,0).isInitialized() &&
                   tmp_onode_data.getArrayData(0,1).isInitialized() &&
                   tmp_onode_data.getArrayData(1,0).isInitialized() &&
                   tmp_onode_data.getArrayData(1,1).isInitialized()) {
               
                  nodeouternodesum2d_(
                     filo(0),filo(1),
                     fihi(0),fihi(1),
                     cilo(0),cilo(1),
                     cihi(0),cihi(1),
                     ratio,
                     node_data->getDepth(),
                     node_gcw(0),node_gcw(1),
                     node_data->getPointer(),     // node data dst
                     tmp_onode_data.getPointer(0,0), // x lower src
                     tmp_onode_data.getPointer(0,1), // x upper src
                     tmp_onode_data.getPointer(1,0), // y lower src
                     tmp_onode_data.getPointer(1,1)); // y upper src
               }

            } // DIM == 2
         
            if (DIM == 3) {
            
               if (tmp_onode_data.getArrayData(0,0).isInitialized() &&
                   tmp_onode_data.getArrayData(0,1).isInitialized() &&
                   tmp_onode_data.getArrayData(1,0).isInitialized() &&
                   tmp_onode_data.getArrayData(1,1).isInitialized() &&
                   tmp_onode_data.getArrayData(2,0).isInitialized() &&
                   tmp_onode_data.getArrayData(2,1).isInitialized()) {
               
                  nodeouternodesum3d_(
                     filo(0),filo(1),filo(2),
                     fihi(0),fihi(1),fihi(2),
                     cilo(0),cilo(1),cilo(2),
                     cihi(0),cihi(1),cihi(2),
                     ratio,
                     node_data->getDepth(),
                     node_gcw(0),node_gcw(1),node_gcw(2),
                     node_data->getPointer(),     // node data dst
                     tmp_onode_data.getPointer(0,0), // x lower src
                     tmp_onode_data.getPointer(0,1), // x upper src
                     tmp_onode_data.getPointer(1,0), // y lower src
                     tmp_onode_data.getPointer(1,1), // y upper src
                     tmp_onode_data.getPointer(2,0), // z lower src
                     tmp_onode_data.getPointer(2,1)); // z upper src
               }

            } // DIM == 3
            
            // If desired, fill "hanging" nodes on fine patch by 
            // linear interpolation between "coarse" nodes on 
            // coarse-fine boundary.
            if (fill_hanging_nodes) {

               for (int ibb1 = 0; ibb1 < num_bdry_boxes; ibb1++) {
   
                  const hier::BoundaryBox<DIM>& bbox = pboundaries[ibb1];
                  const hier::Index<DIM>& bboxilo = bbox.getBox().lower();
                  const hier::Index<DIM>& bboxihi = bbox.getBox().upper();
 
                  const int bbox_loc = bbox.getLocationIndex();
 
                  if (DIM == 2) {
                     nodehangnodeinterp2d_(
                        filo(0),filo(1),
                        fihi(0),fihi(1),
                        cilo(0),cilo(1),
                        cihi(0),cihi(1),
                        bboxilo(0),bboxilo(1),
                        bboxihi(0),bboxihi(1),
                        bbox_loc,
                        ratio,
                        node_data->getDepth(),
                        node_gcw(0),node_gcw(1),
                        node_data->getPointer());
                  }               
 
                  if (DIM == 3) {
                     nodehangnodeinterp3d_(
                        filo(0),filo(1),filo(2),
                        fihi(0),fihi(1),fihi(2),
                        cilo(0),cilo(1),cilo(2),
                        cihi(0),cihi(1),cihi(2),
                        bboxilo(0),bboxilo(1),bboxilo(2),
                        bboxihi(0),bboxihi(1),bboxihi(2),
                        bbox_loc,
                        ratio,
                        node_data->getDepth(),
                        node_gcw(0),node_gcw(1),node_gcw(2),
                        node_data->getPointer());
                  }
                  
               } // iterate over coarse-fine boundary box regions

            } // if fill hanging nodes

         }  //  for (int i ...  iterate over node data 

      }  // if (num_bdry_boxes > 0) .... if patch lies on coarse-fine boundary

   } // iterate over fine level patches

}

/*
*************************************************************************
*                                                                       *
* Private member functions to copy between node and outernode data      *
* over an entire level.                                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryNodeSum<DIM>::copyNodeToOuternodeOnLevel(
   tbox::Pointer<hier::PatchLevel<DIM> > level,
   const tbox::Array<int>& node_data_id,
   const tbox::Array<int>& onode_data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());
#endif

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++ ) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

      for (int i = 0; i < node_data_id.size(); i++) {
         tbox::Pointer< pdat::NodeData<DIM,double> > node_data =
            patch->getPatchData(node_data_id[i]);
         tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_data =
            patch->getPatchData(onode_data_id[i]);

         onode_data->copy(*node_data);
      }
   }

}

template<int DIM> void PatchBoundaryNodeSum<DIM>::copyOuternodeToNodeOnLevel(
   tbox::Pointer<hier::PatchLevel<DIM> > level,
   const tbox::Array<int>& onode_data_id,
   const tbox::Array<int>& node_data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());
#endif

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++ ) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

      for (int i = 0; i < node_data_id.size(); i++) {
         tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_data =
            patch->getPatchData(onode_data_id[i]);
         tbox::Pointer< pdat::NodeData<DIM,double> > node_data =
            patch->getPatchData(node_data_id[i]);

         onode_data->copy2(*node_data);
      }
   }

}


}
}

#endif

