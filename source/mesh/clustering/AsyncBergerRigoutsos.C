/*
 * File:        $RCSfile$
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 602 $
 * Modified:    $Date: 2005-09-06 11:51:31 -0700 (Tue, 06 Sep 2005) $
 * Description: Asynchronous Berger-Rigoutsos algorithm wrapper
 */

#ifndef included_mesh_AsyncBergerRigoutsos_C
#define included_mesh_AsyncBergerRigoutsos_C

#include "AsyncBergerRigoutsos.h"

#include "BoxComm.h"
#include "AsyncBergerRigoutsosNode.h"
#include "tbox/MPI.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/TimerManager.h"

namespace SAMRAI {
namespace mesh {

using namespace tbox;

template<int DIM>
Pointer<Timer> AsyncBergerRigoutsos<DIM>::t_run_abr;
template<int DIM>
Pointer<Timer> AsyncBergerRigoutsos<DIM>::t_globalize_boxes;


/*
************************************************************************
* Constructor stores parameters of options for ussing                  *
* the asynchronous Berger-Rigoutsos implementation.                    *
************************************************************************
*/
template<int DIM>
AsyncBergerRigoutsos<DIM>::AsyncBergerRigoutsos(
   Pointer<Database> database )
   : d_mpi_communicator(MPI::commNull),
     d_log_node_history(false),
     d_log_cluster_summary(false),
     d_owner_mode("MOST_OVERLAP"),
     d_use_level_boxes(false),
     d_algo_advance_mode("ADVANCE_SOME"),
     d_max_gcw(1)
{

   bool use_private_communicator = true;
   /*
    * Set database-dependent parameters or cache them for use
    * when we construct a dendogram root.
    */
   if ( ! database.isNull() ) {
      use_private_communicator =
         database->getBoolWithDefault( "use_private_communicator",
                                       use_private_communicator );
      d_log_node_history =
         database->getBoolWithDefault( "log_node_history",
                                       d_log_node_history );
      d_log_cluster_summary =
         database->getBoolWithDefault( "log_cluster_summary",
                                       d_log_cluster_summary );
      d_use_level_boxes =
         database->getBoolWithDefault( "use_level_boxes",
                                       d_use_level_boxes );
      d_algo_advance_mode =
         database->getStringWithDefault( "algo_advance_mode",
                                         d_algo_advance_mode );
      if ( database->isInteger("max_gcw") ) {
         database->getIntegerArray( "max_gcw",
                                    (int*)(d_max_gcw),
                                    DIM );
      }
      d_owner_mode =
         database->getStringWithDefault( "owner_mode",
                                         d_owner_mode );
   }

#ifdef HAVE_MPI
   if ( use_private_communicator ) {
      MPI_Comm_dup( MPI::getCommunicator(), &d_mpi_communicator );
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( d_mpi_communicator != MPI::commNull );
#endif
   }
#endif

   if ( t_run_abr.isNull() ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( t_globalize_boxes.isNull() );
#endif
      t_run_abr = TimerManager::getManager()->
         getTimer("mesh::AsyncBergerRigoutsos::run_abr");
      t_globalize_boxes = TimerManager::getManager()->
         getTimer("mesh::AsyncBergerRigoutsos::globalize_boxes");
      ShutdownRegistry::registerShutdownRoutine(freeTimers, 254);
   }

   return;
}


template<int DIM>
AsyncBergerRigoutsos<DIM>::~AsyncBergerRigoutsos(void)
{
#ifdef HAVE_MPI
   if ( d_mpi_communicator != MPI::commNull ) {
      // Free the private communicator (if MPI has not been finalized).
      int flag;
      MPI_Finalized(&flag);
      if ( ! flag ) {
         MPI_Comm_free( &d_mpi_communicator );
      }
   }
#endif

   return;
}


/*
************************************************************************
*                                                                      *
* Implement the mesh::BoxGeneratorStrategy<DIM> interface method using *
* the asynchronous Berger-Rigoutsos implementation.                    *
*                                                                      *
* Create objects for using the ABR recursion tree, set options for     *
* using the ABR implementation, then run it.                           *
*                                                                      *
* The output boxes from the dendogram root is in the form of a         *
* hier::LayerNodeSet<DIM>.  This method postprocess that data to       *
* convert the output to the box list form required by the              *
* box clustering strategy interface.                                   *
*                                                                      *
************************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsos<DIM>::findBoxesContainingTags(
   hier::BoxList<DIM> &boxes,
   const Pointer<hier::PatchLevel<DIM> > level,
   const int tag_data_index,
   const int tag_val,
   const hier::Box<DIM>& bound_box,
   const hier::IntVector<DIM>& min_box,
   const double efficiency_tol,
   const double combine_tol) const
{
   if ( bound_box.empty() ) {
      TBOX_ERROR("AsyncBergerRigoutsos: empty bounding box not allowed.");
   }

   const int nprocs = MPI::getNodes();

   typename AsyncBergerRigoutsosNode<DIM>::CommonParams
      cp( level,
          tag_data_index,
          tag_val,
          min_box,
          efficiency_tol,
          combine_tol,
          d_mpi_communicator == MPI::commNull ?
          MPI::getCommunicator() : d_mpi_communicator );

   AsyncBergerRigoutsosNode<DIM> root( &cp, &bound_box );

   root.setComputeEdges( d_use_level_boxes ? 0 : 2 );
   root.setLogNodeHistory(d_log_node_history);
   root.setOwnerMode( d_owner_mode );
   root.setUseLevelBoxes( d_use_level_boxes );

   cp.job_relauncher.setAlgorithmAdvanceMode( d_algo_advance_mode );

   t_run_abr->start();
   root.runAlgorithm();
   t_run_abr->stop();

#ifdef DEBUG_CHECK_ASSERTIONS
   assertNoMessageForPrivateCommunicator();
#endif

   if ( d_log_cluster_summary ) {
      /*
       * Log summary of clustering and dendogram.
       * Maybe this should have a flag to turn on/off.
       */
      plog << "Async BR on proc " << MPI::getRank() << " owned "
           << root.getMaxOwnership() << " of "
           << root.getMaxNodes() << " nodes ("
           << (double)root.getMaxOwnership()/root.getMaxNodes()
           << ") in " << root.getMaxGeneration() << " generations,"
           << "   " << root.getNumBoxesGenerated() << " boxes generated.\n"
           << "Number of continuations: avg = " << root.getAvgNumberOfCont()
           << "   max = " << root.getMaxNumberOfCont() << "\n";
      ;
   }


   t_globalize_boxes->start();


   hier::LayerNodeSet<DIM> layer_node_set = root.getNewNodes();

   /*
    * In single onwer mode, process 0 owns all the boxes,
    * so we extract the boxes on process 0 and broadcast
    * to other processes.
    *
    * Otherwise, the output boxes are distributed over
    * participating processes.
    * We globalize the output layer node,
    * which duplicates the nodes everywhere.
    * Then each process extract the box list
    * from duplicated set of layer nodes.
    */
   if ( d_owner_mode == "SINGLE_OWNER" ) {

      if ( MPI::getRank() == 0 ) {

         boxes.clearItems();
         int num_boxes = 
            layer_node_set.getNodeContainer().size();
         hier::BoxArray<DIM> new_box_array(num_boxes);
         const typename hier::LayerNodeSet<DIM>::NodeContainer &new_nodes =
            layer_node_set.getNodeContainer();

         typename hier::LayerNodeSet<DIM>::NodeContainer::const_iterator ni;
         for ( ni=new_nodes.begin(); ni!=new_nodes.end(); ++ni ) {
            const typename hier::LayerNodeSet<DIM>::Node &node = *ni;
            boxes.appendItem( node.getBox() );
         }

      }

      hier::BoxComm<DIM>::bcastBoxList(boxes);

   }
   else {

      layer_node_set.setParallelState(hier::LayerNodeSet<DIM>::GLOBALIZED);

      int n;
      int num_boxes = 0;
      for ( n=0; n<nprocs; ++n ) {
         num_boxes += layer_node_set.getNodeContainer(n).size();
      }
      hier::BoxArray<DIM> new_box_array(num_boxes);

      int box_id = 0;
      for ( n=0; n<nprocs; ++n ) {
         const typename hier::LayerNodeSet<DIM>::NodeContainer &new_nodes =
            layer_node_set.getNodeContainer(n);
         typename hier::LayerNodeSet<DIM>::NodeContainer::const_iterator ni;
         for ( ni=new_nodes.begin(); ni!=new_nodes.end(); ++ni ) {
            const typename hier::LayerNodeSet<DIM>::Node &node = *ni;
            new_box_array(box_id++) = node.getBox();
         }
      }
      boxes = hier::BoxList<DIM>(new_box_array);

   }

   t_globalize_boxes->stop();

   if ( d_log_cluster_summary ) {
      int number_cells = 0;
      for ( typename hier::BoxList<DIM>::Iterator bi(boxes); bi; bi++ ) {
         number_cells += bi().size();
      }
      plog << "After globalize: "
           << number_cells << " cells in "
           << boxes.size() << " boxes, averaging "
           << double(number_cells)/boxes.size() << " cells/box\n";
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assertNoMessageForPrivateCommunicator();
#endif

   return;
}



/*
***************************************************************************
*                                                                         *
***************************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsos<DIM>::assertNoMessageForPrivateCommunicator() const
{
#ifdef HAVE_MPI
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if ( d_mpi_communicator != MPI::commNull ) {
      int flag;
      MPI::status mpi_status;
      int mpi_err = MPI_Iprobe( MPI_ANY_SOURCE,
                                MPI_ANY_TAG,
                                d_mpi_communicator,
                                &flag,
                                &mpi_status );
      if ( mpi_err != MPI_SUCCESS ) {
         TBOX_ERROR("Error probing for possible lost messages.");
      }
      if ( flag == true ) {
         int count = -1;
         mpi_err = MPI_Get_count( &mpi_status, MPI_INT, &count );
         TBOX_ERROR("Library error!\n"
                    <<"AsyncBergerRigoutsos detected before or after\n"
                    <<"running AsyncBergerRigoutsosNode that there\n"
                    <<"is a message yet to be received.  This is\n"
                    <<"an error because all messages using the\n"
                    <<"private communicator should have been\n"
                    <<"accounted for.  Message status:\n"
                    <<"source " << mpi_status.MPI_SOURCE << '\n'
                    <<"tag " << mpi_status.MPI_TAG << '\n'
                    <<"count " << count << " (assuming integers)\n");
      }
   }
#endif
   return;
}




/*
***************************************************************************
*                                                                         *
* Release static timers.  To be called by shutdown registry to make sure  *
* memory for timers does not leak.                                        *
*                                                                         *
***************************************************************************
*/
template<int DIM>
void AsyncBergerRigoutsos<DIM>::freeTimers()
{
   t_run_abr.setNull();
   t_globalize_boxes.setNull();
   return;
}


}
}
#endif
