//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/algorithm/femutils/standard/MblkPatchBoundaryNodeSum.C $
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description:	Routines for summing node data at patch boundaries
//

#ifndef included_algs_MblkPatchBoundaryNodeSum_C
#define included_algs_MblkPatchBoundaryNodeSum_C

#include "MblkPatchBoundaryNodeSum.h"

#include "VariableDatabase.h"
#include "NodeData.h"
#include "NodeDataFactory.h"
#include "NodeGeometry.h"
#include "OuternodeData.h"
#include "OuternodeDoubleConstantCoarsen.h"
#include "OuternodeSumTransactionFactory.h"
#include "CoarsenAlgorithm.h"
#include "MultiblockRefineAlgorithm.h"
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

template<int DIM> int MblkPatchBoundaryNodeSum<DIM>::s_instance_counter = 0;

template<int DIM> tbox::Array< tbox::Array<int> > 
   MblkPatchBoundaryNodeSum<DIM>::s_onode_src_id_array =
      tbox::Array< tbox::Array<int> >(0);
template<int DIM> tbox::Array< tbox::Array<int> > 
   MblkPatchBoundaryNodeSum<DIM>::s_onode_dst_id_array =
      tbox::Array< tbox::Array<int> >(0);

/*
*************************************************************************
*                                                                       *
* Constructor patch boundary node sum objects initializes data members  *
* to default (undefined) states.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> MblkPatchBoundaryNodeSum<DIM>::MblkPatchBoundaryNodeSum(
   const std::string& object_name,
   tbox::Pointer< hier::MultiblockPatchHierarchy<DIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
#endif

   d_object_name = object_name;
   d_setup_called = false;

   d_num_reg_sum = 0;

   d_level.setNull();
   d_hierarchy = hierarchy;
   d_coarsest_level = -1;
   d_finest_level = -1;

   d_level_setup_called     = false;
   d_hierarchy_setup_called = false;

   d_sum_transaction_factory = new algs::OuternodeSumTransactionFactory<DIM>();

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

template<int DIM> MblkPatchBoundaryNodeSum<DIM>::~MblkPatchBoundaryNodeSum()
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

template<int DIM> void MblkPatchBoundaryNodeSum<DIM>::registerSum(
   int node_data_id)
{

   if (d_setup_called) {
      TBOX_ERROR("MblkPatchBoundaryNodeSum<DIM>::register error..."
                 << "\nobject named " << d_object_name
                 << "\nCannot call registerSum with this MblkPatchBoundaryNodeSum"
                 << "\nobject since it has already been used to create communication"
                 << "\nschedules; i.e., setupSum() has been called."
                 << std::endl);
   }

   if (node_data_id < 0) {
      TBOX_ERROR("MblkPatchBoundaryNodeSum<DIM> register error..."
                 << "\nobject named " << d_object_name
                 << "\n node_data_id = " << node_data_id
                 << " is an invalid patch data identifier." << std::endl);
   }

   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   tbox::Pointer< pdat::NodeDataFactory<DIM,double> > node_factory =
      var_db->getPatchDescriptor()->getPatchDataFactory(node_data_id);

   if (node_factory.isNull()) {

      TBOX_ERROR("MblkPatchBoundaryNodeSum<DIM> register error..."
                 << "\nobject named " << d_object_name
                 << "\n node_data_id = " << node_data_id
                 << " does not correspond to node data of type double." << std::endl);

   } else {

      static std::string tmp_onode_src_variable_name("MblkPatchBoundaryNodeSum__internal-onode-src");
      static std::string tmp_onode_dst_variable_name("MblkPatchBoundaryNodeSum__internal-onode-dst");

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

      std::string var_suffix = 	 tbox::Utilities::intToString(data_depth_id, 4) + 
	 "__depth=" + tbox::Utilities::intToString(data_depth);

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

template<int DIM> void MblkPatchBoundaryNodeSum<DIM>::setupSum(
   tbox::Pointer<hier::MultiblockPatchLevel<DIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   if (d_hierarchy_setup_called) {

      TBOX_ERROR("MblkPatchBoundaryNodeSum<DIM>::setupSum error...\n"
         << " object named " << d_object_name
         << " already initialized via setupSum() for hierarchy" << std::endl);

   } else {

      d_setup_called = true;
      d_level_setup_called = true;

      d_level = level;

      d_single_level_sum_schedule.resizeArray(1);

      // Communication algorithm for summing outernode values on a level
      tbox::Pointer< xfer::RefineAlgorithm<DIM> > single_level_sum_algorithm =
         new xfer::RefineAlgorithm<DIM>();

      for (int i = 0; i < d_num_reg_sum; i++) {
         single_level_sum_algorithm->registerRefine(
            d_onode_dst_id[i],  // dst data
            d_onode_src_id[i],  // src data
            d_onode_dst_id[i],  // scratch data
            (xfer::RefineOperator<DIM>*)NULL);
      }

      xfer::MultiblockRefineAlgorithm<DIM> mblk_sum_algorithm(
         single_level_sum_algorithm, d_hierarchy);

      d_single_level_sum_schedule[0] =
         mblk_sum_algorithm.createSchedule(
            d_level,
            (xfer::MultiblockRefinePatchStrategy<DIM>*)NULL,
            d_sum_transaction_factory);

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

template<int DIM> void MblkPatchBoundaryNodeSum<DIM>::computeSum(
   const bool fill_hanging_nodes) const
{
   NULL_USE(fill_hanging_nodes);

   if (d_level_setup_called) {

      d_level->allocatePatchData(d_onode_src_data_set);
      d_level->allocatePatchData(d_onode_dst_data_set);

      doLevelSum(d_level);

      d_level->deallocatePatchData(d_onode_src_data_set);
      d_level->deallocatePatchData(d_onode_dst_data_set);

   }

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

template<int DIM> void MblkPatchBoundaryNodeSum<DIM>::doLevelSum(
   tbox::Pointer<hier::MultiblockPatchLevel<DIM> > level) const
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
* Private member functions to copy between node and outernode data      *
* over an entire level.                                                 *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MblkPatchBoundaryNodeSum<DIM>::copyNodeToOuternodeOnLevel(
   tbox::Pointer<hier::MultiblockPatchLevel<DIM> > level,
   const tbox::Array<int>& node_data_id,
   const tbox::Array<int>& onode_data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());
#endif

   for (int bn = 0; bn < level->getNumberOfBlocks(); bn++) {
      tbox::Pointer<hier::PatchLevel<DIM> > patch_level =
         level->getPatchLevelForBlock(bn);

      if (!(patch_level.isNull())) {

         for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
              ip; ip++ ) {
            tbox::Pointer<hier::Patch<DIM> > patch =
               patch_level->getPatch(ip());

            for (int i = 0; i < node_data_id.size(); i++) {
               tbox::Pointer< pdat::NodeData<DIM,double> > node_data =
                  patch->getPatchData(node_data_id[i]);
               tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_data =
                  patch->getPatchData(onode_data_id[i]);

               onode_data->copy(*node_data);
            }
         }
      }
   }
}

template<int DIM> void MblkPatchBoundaryNodeSum<DIM>::copyOuternodeToNodeOnLevel(
   tbox::Pointer<hier::MultiblockPatchLevel<DIM> > level,
   const tbox::Array<int>& onode_data_id,
   const tbox::Array<int>& node_data_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());
#endif

   for (int bn = 0; bn < level->getNumberOfBlocks(); bn++) {
      tbox::Pointer<hier::PatchLevel<DIM> > patch_level =
         level->getPatchLevelForBlock(bn);

      if (!(patch_level.isNull())) {

         for (typename hier::PatchLevel<DIM>::Iterator ip(patch_level);
              ip; ip++ ) {
            tbox::Pointer<hier::Patch<DIM> > patch =
               patch_level->getPatch(ip());

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

}
}

#endif

