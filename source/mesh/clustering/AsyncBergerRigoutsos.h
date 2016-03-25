/*
 * File:         $RCSfile$
 * Copyright:    (c) 1997-2005 The Regents of the University of California
 * Revision:     $Revision: 279 $
 * Modified:     $Date: 2005-03-31 13:08:56 -0800 (Thu, 31 Mar 2005) $
 * Description:  Asynchronous Berger-Rigoutsos clustering algorithm.
 */

#ifndef included_mesh_AsyncBergerRigoutsos
#define included_mesh_AsyncBergerRigoutsos

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif

#ifndef included_mesh_BoxGeneratorStrategy
#include "BoxGeneratorStrategy.h"
#endif

#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
namespace mesh {


/*!
 * @brief Asynchronous Berger-Rigoutsos implementation.
 * This class is derived from the abstract base class
 * mesh::BoxGeneratorStrategy<DIM>.  Thus, it serves as a concrete
 * implementation of the box generator Strategy pattern interface.
 *
 * This class uses the mesh::AsyncBergerRigoutsosNode class
 * to carry out the Berger-Rigoutsos algorithm in parallel.
 * It handles aspects not related that algorithm.  It:
 * - Implements the box generator Strategy pattern interface.
 * - Provides an interface with the input database for setting
 *   parameters influencing the implementation.
 * - Converts output data into the format required by
 *   mesh::BoxGeneratorStrategy<DIM>.
 * For more details on the parallel implementation,
 * see mesh::AsyncBergerRigoutsosNode.
 *
 * User inputs (default):
 * - bool @b use_private_communicator (true):
 *   Whether to create and cache a private MPI communicator.
 *   A private communicator helps to ensure that the complex
 *   communications used by the asynchronous algorithm is
 *   isolated from other communications.
 *   If a private communicator is used, it is duplicated from
 *   MPI::getCommunicator() at construction time.
 *   If a private communicator is not used, the current value of
 *   MPI::getCommunicator() is used when a communicator is needed.
 *   The private communicator is freed when the object goes out of
 *   scope.  Disabling the use of a private communicator is not
 *   recommended unless you are absolutely sure unintended message
 *   reception will be avoided.
 * - bool @b use_level_boxes (false):
 *   Whether to use and compute
 *   the local copy of all boxes on the patch level.
 * - string @b algo_advance_mode ("ADVANCE_SOME"):
 *   Asynchronous algorithm advance mode.  The default has been
 *   empirically determined to scale best to higher numbers of
 *   processors and work adequately for lower numbers of processors.
 * - string @b owner_mode ("MOST_OVERLAP"):
 *   How to chose the owner from a dendogram node group.
 *   This string is used in AsyncBergerRigoutsosNode::setOwnerMode().
 * - IntVector @b max_gcw (1):
 *   Max ghost width for overlap checks.
 *   It should be set tho the max ghost cell
 *   width used in the problem in order to compute correct edges.
 *   This can be zero for now, as SAMRAI does not use the edge data.
 *
 * Debugging inputs (default):
 * - bool @b log_node_history (false):
 *   Whether to log what certain actions of nodes in the dendogram.
 *   This degrades the performance but is a very useful debugging
 *   tool.
 * - bool @b log_cluster_summary (false):
 *   Whether to log the results of the clustering.
 */
template<int DIM>
class AsyncBergerRigoutsos : public mesh::BoxGeneratorStrategy<DIM>
{

public:

   /*!
    * @brief Constructor.
    */
   AsyncBergerRigoutsos( tbox::Pointer<tbox::Database> database );

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~AsyncBergerRigoutsos(void);


   /*!
    * @brief Implement the mesh::BoxGeneratorStrategy<DIM> interface
    * method of the same name.
    *
    * Create a list of boxes that covers all integer tags on
    * the patch level that match the specified tag value.
    * Each box will be at least as large as the given minimum
    * size and the tolerances will be met.
    *
    * The efficiency tolerance is a threshold value for the percentage of
    * tagged cells in each box.  If this percentage is below the tolerance,
    * the box will continue to be split into smaller boxes.
    *
    * The combine tolerance is a threshold value for the sum of the volumes
    * of two boxes into which a box may be potentially split.  If ratio of
    * that sum and the volume of the original box, the box will not be split.
    *
    * This function is actually a switch for selecting one of several
    * variations of the computational implementation of the algorithm.  
    * See the discussion above for a desription of the different options.
    * All implementations should generate identical results, but the 
    * performance may vary on different systems.  By default, the 
    * ORIGINAL algorithm is used when running on a single processor, and
    * BINARY_TREE is used when running on multiple processors. The user
    * may reset these options through input.  See the discussion above
    * for more information.
    */
   void findBoxesContainingTags(hier::BoxList<DIM>& boxes,
                                const tbox::Pointer<hier::PatchLevel<DIM> > level,
                                const int tag_data_index,
                                const int tag_val,
                                const hier::Box<DIM>& bound_box,
                                const hier::IntVector<DIM>& min_box,
                                const double efficiency_tol,
                                const double combine_tol) const; 
   

private:


   void assertNoMessageForPrivateCommunicator() const;

   /*!
    * Free static timers.
    *
    * To be called by shutdown registry to make sure
    * memory for timers does not leak.
    */
   static void freeTimers();


   //! @brief Communicator, which may be privately created.
   tbox::MPI::comm d_mpi_communicator;

   //! @brief Whether to log execution node allocation and deallocation.
   bool d_log_node_history;

   //! @brief Whether to log cluster summary.
   bool d_log_cluster_summary;

   //! @brief How to select the owner of a node.
   string d_owner_mode;

   //! @brief Whether to take advantage of level boxes data.
   bool d_use_level_boxes;

   //! @brief Asynchronous mode for advancing algorithm.
   string d_algo_advance_mode;

   //! @brief Max ghost cell width for overlap computation.
   hier::IntVector<DIM> d_max_gcw;

   static tbox::Pointer<tbox::Timer> t_run_abr;
   static tbox::Pointer<tbox::Timer> t_globalize_boxes;

};

}
}

#ifndef DEBUG_NO_INLINE
// #include "AsyncBergerRigoutsos.I"
#endif

#endif  // included_mesh_AsyncBergerRigoutsos

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "AsyncBergerRigoutsos.C"
#endif
