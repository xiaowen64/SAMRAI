//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/femutils/standard/PatchBoundaryEdgeSum.C $
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Routines for summing edge data at patch boundaries
//


#ifndef included_algs_PatchBoundaryEdgeSum_C
#define included_algs_PatchBoundaryEdgeSum_C

#include "PatchBoundaryEdgeSum.h"

#include "VariableDatabase.h"
#include "EdgeData.h"
#include "EdgeDataFactory.h"
#include "OuteredgeData.h"
#include "OuteredgeSumTransactionFactory.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"
#include "RefineOperator.h"
#include "tbox/Utilities.h"


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

template<int DIM> int PatchBoundaryEdgeSum<DIM>::s_instance_counter = 0;

template<int DIM> tbox::Array< tbox::Array<int> >
   PatchBoundaryEdgeSum<DIM>::s_oedge_src_id_array = 
      tbox::Array< tbox::Array<int> >(0);
template<int DIM> tbox::Array< tbox::Array<int> >
   PatchBoundaryEdgeSum<DIM>::s_oedge_dst_id_array = 
      tbox::Array< tbox::Array<int> >(0);

/*
*************************************************************************
*                                                                       *
* Static functions to determine number of patch data slots needed       *
* for PatchBoundaryEdgeSum objects.                                     *
*                                                                       *
*************************************************************************
*/
 
template<int DIM>
int
PatchBoundaryEdgeSum<DIM>::getNumSharedPatchDataSlots(
   int max_variables_to_register)
{
   // edge boundary sum requires two internal outeredge variables
   // (source and destination) for each registered variable.
 
   return( 2 * max_variables_to_register );
}
 
template<int DIM>
int
PatchBoundaryEdgeSum<DIM>::getNumUniquePatchDataSlots(
   int max_variables_to_register)
{
   // all patch data slots used by edge boundary sum are static
   // and shared among all objects.
 
   return( 0 );
}

/*
*************************************************************************
*                                                                       *
* Constructor patch boundary edge sum objects initializes data members  *
* to default (undefined) states.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> PatchBoundaryEdgeSum<DIM>::PatchBoundaryEdgeSum(
   const std::string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
#endif

   d_object_name = object_name;
   d_setup_called = false;

   d_num_reg_sum = 0;

   d_level.setNull();

   d_sum_transaction_factory = new OuteredgeSumTransactionFactory<DIM>();

   s_instance_counter++;
}

/*
*************************************************************************
*                                                                       *
* Destructor removes temporary outeredge patch data ids from            *
* variable database, if defined.                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> PatchBoundaryEdgeSum<DIM>::~PatchBoundaryEdgeSum()
{

   s_instance_counter--;
   if (s_instance_counter == 0) {
      const int arr_length_depth = s_oedge_src_id_array.size();
      for (int id = 0; id < arr_length_depth; id++) {
         const int arr_length_nvar = s_oedge_src_id_array[id].size();

         for (int iv = 0; iv < arr_length_nvar; iv++) {

            if (s_oedge_src_id_array[id][iv] >= 0) {
               hier::VariableDatabase<DIM>::getDatabase()->
                  removeInternalSAMRAIVariablePatchDataIndex(
                     s_oedge_src_id_array[id][iv]);
            }
            if (s_oedge_dst_id_array[id][iv] >= 0) {
               hier::VariableDatabase<DIM>::getDatabase()->
                  removeInternalSAMRAIVariablePatchDataIndex(
                     s_oedge_dst_id_array[id][iv]);
            }

            s_oedge_src_id_array[id].resizeArray(0);
            s_oedge_dst_id_array[id].resizeArray(0);

         }
      }

      s_oedge_src_id_array.resizeArray(0);
      s_oedge_dst_id_array.resizeArray(0);
   }

}

/*
*************************************************************************
*                                                                       *
* Register edge patch data index for summation.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryEdgeSum<DIM>::registerSum(
   int edge_data_id)
{

   if (d_setup_called) {
      TBOX_ERROR("PatchBoundaryEdgeSum<DIM>::register error..."
                 << "\nobject named " << d_object_name 
                 << "\nCannot call registerSum with this PatchBoundaryEdgeSum"
                 << "\nobject since it has already been used to create communication"
                 << "\nschedules; i.e., setupSum() has been called."
                 << std::endl);
   }

   if (edge_data_id < 0) {
      TBOX_ERROR("PatchBoundaryEdgeSum<DIM> register error..."
                 << "\nobject named " << d_object_name 
                 << "\n edge_data_id = " << edge_data_id
                 << " is an invalid patch data identifier." << std::endl);
   }

   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   tbox::Pointer< pdat::EdgeDataFactory<DIM,double> > edge_factory =
      var_db->getPatchDescriptor()->getPatchDataFactory(edge_data_id);

   if (edge_factory.isNull()) {

      TBOX_ERROR("PatchBoundaryEdgeSum<DIM> register error..."
                 << "\nobject named " << d_object_name 
                 << "\n edge_data_id = " << edge_data_id
                 << " does not correspond to edge data of type double." << std::endl);

   } else {

      static std::string tmp_oedge_src_variable_name("PatchBoundaryEdgeSum__internal-oedge-src");
      static std::string tmp_oedge_dst_variable_name("PatchBoundaryEdgeSum__internal-oedge-dst");

      const int reg_sum_id = d_num_reg_sum;

      d_num_reg_sum++;

      d_user_edge_data_id.resizeArray(d_num_reg_sum);
         d_user_edge_data_id[reg_sum_id] = ID_UNDEFINED;
      d_user_edge_depth.resizeArray(d_num_reg_sum);
         d_user_edge_depth[reg_sum_id] = ID_UNDEFINED;
      d_tmp_oedge_src_variable.resizeArray(d_num_reg_sum);
      d_tmp_oedge_dst_variable.resizeArray(d_num_reg_sum);
      d_oedge_src_id.resizeArray(d_num_reg_sum);
         d_oedge_src_id[reg_sum_id] = ID_UNDEFINED;
      d_oedge_dst_id.resizeArray(d_num_reg_sum);
         d_oedge_dst_id[reg_sum_id] = ID_UNDEFINED;

      const int data_depth = edge_factory->getDefaultDepth();
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

      if (s_oedge_src_id_array.size() < array_by_depth_size) {
         s_oedge_src_id_array.resizeArray(array_by_depth_size);
         s_oedge_dst_id_array.resizeArray(array_by_depth_size);
      }

      if (s_oedge_src_id_array[data_depth].size() < num_data_at_depth) {
         const int old_size = s_oedge_src_id_array[data_depth].size();
         const int new_size = num_data_at_depth;

         s_oedge_src_id_array[data_depth].resizeArray(new_size);
         s_oedge_dst_id_array[data_depth].resizeArray(new_size);
         for (int i = old_size; i < new_size; i++) {
            s_oedge_src_id_array[data_depth][i] = ID_UNDEFINED;
            s_oedge_dst_id_array[data_depth][i] = ID_UNDEFINED;
         }
      }

      char var_suffix[17];
      sprintf(var_suffix, "%04d__depth=%04d", data_depth_id, data_depth);

      std::string toedge_src_var_name = tmp_oedge_src_variable_name + var_suffix;
      d_tmp_oedge_src_variable[reg_sum_id] = var_db->getVariable(toedge_src_var_name);
      if (d_tmp_oedge_src_variable[reg_sum_id].isNull()) {
         d_tmp_oedge_src_variable[reg_sum_id] =
            new pdat::OuteredgeVariable<DIM,double>(toedge_src_var_name, data_depth);
      }

      std::string toedge_dst_var_name = tmp_oedge_dst_variable_name + var_suffix;
      d_tmp_oedge_dst_variable[reg_sum_id] = var_db->getVariable(toedge_dst_var_name);
      if (d_tmp_oedge_dst_variable[reg_sum_id].isNull()) {
         d_tmp_oedge_dst_variable[reg_sum_id] =
            new pdat::OuteredgeVariable<DIM,double>(toedge_dst_var_name, data_depth);
      }

      if ( s_oedge_src_id_array[data_depth][data_depth_id] < 0 ) {
         s_oedge_src_id_array[data_depth][data_depth_id] =
            var_db->registerInternalSAMRAIVariable(
                    d_tmp_oedge_src_variable[reg_sum_id], 
                    hier::IntVector<DIM>(0));
      }
      if ( s_oedge_dst_id_array[data_depth][data_depth_id] < 0) {
         s_oedge_dst_id_array[data_depth][data_depth_id] =
            var_db->registerInternalSAMRAIVariable(
                    d_tmp_oedge_dst_variable[reg_sum_id], 
                    hier::IntVector<DIM>(0));
      }

      d_user_edge_data_id[reg_sum_id] = edge_data_id;
      d_user_edge_depth[reg_sum_id] = data_depth;

      d_num_registered_data_by_depth[data_depth] = num_data_at_depth;

      d_oedge_src_id[reg_sum_id] = s_oedge_src_id_array[data_depth][data_depth_id];
      d_oedge_dst_id[reg_sum_id] = s_oedge_dst_id_array[data_depth][data_depth_id];

      d_oedge_src_data_set.setFlag(d_oedge_src_id[reg_sum_id]);
      d_oedge_dst_data_set.setFlag(d_oedge_dst_id[reg_sum_id]);

   }

}

/*
*************************************************************************
*                                                                       *
* Set up schedule to sum edge data around patch boundaries              *
* on a single level.                                                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryEdgeSum<DIM>::setupSum(
   tbox::Pointer<hier::PatchLevel<DIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   d_setup_called = true;

   d_level = level;

   // Communication algorithm for summing outeredge values on a level
   xfer::RefineAlgorithm<DIM> single_level_sum_algorithm;

   for (int i = 0; i < d_num_reg_sum; i++) {
      single_level_sum_algorithm.registerRefine(d_oedge_dst_id[i],  // dst data
                                                d_oedge_src_id[i],  // src data
                                                d_oedge_dst_id[i],  // scratch data
                                                (xfer::RefineOperator<DIM>*)NULL);
   }

   d_single_level_sum_schedule =
      single_level_sum_algorithm.createSchedule(
                                 d_level,
                                 (xfer::RefinePatchStrategy<DIM>*)NULL,
                                 d_sum_transaction_factory);

} 

/*
*************************************************************************
*                                                                       *
* Perform patch boundary edge sum across single level or multiple       *
* hierarchy levels depending on how object was initialized.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryEdgeSum<DIM>::computeSum() const
{
   d_level->allocatePatchData(d_oedge_src_data_set);
   d_level->allocatePatchData(d_oedge_dst_data_set);

   doLevelSum(d_level);
   
   d_level->deallocatePatchData(d_oedge_src_data_set);
   d_level->deallocatePatchData(d_oedge_dst_data_set);
   
}

/*
*************************************************************************
*                                                                       *
* Private member function that performs edge sum across single level.   *
*                                                                       *
* 1. Copy edge data to local outeredge data.                            *
* 2. Transfer and sum outeredge data.                                   *
* 3. Copy outeredge data back to edge data.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void PatchBoundaryEdgeSum<DIM>::doLevelSum(
   tbox::Pointer<hier::PatchLevel<DIM> > level) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
#endif

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++ ) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip());

      for (int i = 0; i < d_user_edge_data_id.size(); i++) {
         tbox::Pointer< pdat::EdgeData<DIM,double> > edge_data =
            patch->getPatchData(d_user_edge_data_id[i]);
         tbox::Pointer< pdat::OuteredgeData<DIM,double> > oedge_data =
            patch->getPatchData(d_oedge_src_id[i]);

         oedge_data->copy(*edge_data);

      }
   }

   d_single_level_sum_schedule->fillData(0.0, false);

   for (typename hier::PatchLevel<DIM>::Iterator ip2(level); ip2; ip2++ ) {
      tbox::Pointer<hier::Patch<DIM> > patch = level->getPatch(ip2());

      for (int i = 0; i < d_user_edge_data_id.size(); i++) {
         tbox::Pointer< pdat::EdgeData<DIM,double> > edge_data =
            patch->getPatchData(d_user_edge_data_id[i]);
         tbox::Pointer< pdat::OuteredgeData<DIM,double> > oedge_data =
            patch->getPatchData(d_oedge_dst_id[i]);

         oedge_data->copy2(*edge_data);

      }
   }

}

}
}

#endif


