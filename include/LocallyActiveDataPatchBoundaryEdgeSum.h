//
// File:	LocallyActiveDataPatchBoundaryEdgeSum.h
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 703 $
// Modified:	$Date: 2005-11-03 14:46:35 -0800 (Thu, 03 Nov 2005) $
// Description:	Routines for summing locally-active edge data at patch boundaries
//
 
#ifndef included_algs_LocallyActiveDataPatchBoundaryEdgeSum
#define included_algs_LocallyActiveDataPatchBoundaryEdgeSum

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_LocallyActiveDataPatchLevelManager
#include "LocallyActiveDataPatchLevelManager.h"
#endif
#ifndef included_pdat_EdgeVariable
#include "EdgeVariable.h"
#endif
#ifndef included_pdat_OuteredgeVariable
#include "OuteredgeVariable.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_xfer_LocallyActiveDataRefineSchedule
#include "LocallyActiveDataRefineSchedule.h"
#endif
#ifndef included_xfer_LocallyActiveDataRefineTransactionFactory
#include "LocallyActiveDataRefineTransactionFactory.h"
#endif

namespace SAMRAI {
    namespace algs {

 /*!
 *  @brief Class LocallyActiveDataPatchBoundaryEdgeSum provides operations summing 
 *  locally-active edge data values at edges that are shared by multiple patches 
 *  on a single level.  Note that this utility only works on a SINGLE patch level, not 
 *  on a multiple levels in an AMR patch hierarchy like the LocallyActiveDataPatchBoundaryNodeSum 
 *  class. Unlike node data, edge data at coarse-fine boundaries are not co-located, so the sum 
 *  operation is not clearly defined.
 *
 *  Usage of a patch boundry edge sum involves the following sequence of steps:
 *
 *  -# Construct a patch boundry edge sum object.  For example,
 *     \verbatim
 *         LocallyActiveDataPatchBoundaryEdgeSum<DIM> my_edge_sum("My Edge Sum");
 *     \endverbatim
 *  -# Register edge data quantities to sum.  For example,
 *     \verbatim
 *         my_edge_sum.registerSum(edge_data_id1);
 *         my_edge_sum.registerSum(edge_data_id2);
 *         etc...
 *     \endverbatim
 *  -# Setup the sum operations for a single level involving patches defined to be active
 *     by a given level manager.  For example,
 *     \verbatim
 *         my_edge_sum.setupSum(level, level_mgr);
 *     \endverbatim
 *  -# Execute the sum operation.  For example,
 *     \verbatim
 *         my_edge_sum.computeSum()
 *     \endverbatim
 *
 *  The result of these operations is that each edge patch data value associated
 *  with the registered ids at patch boundaries on the level is replaced by the 
 *  sum of all data values at the edge.
 */

template<int DIM> class LocallyActiveDataPatchBoundaryEdgeSum
{
public:
      /*!
    *  @brief Constructor initializes object to default (mostly undefined) 
    *  state.
    *
    *  @param object_name const string reference for name of object used 
    *  in error reporting.  When assertion checking is on, the string
    *  cannot be empty.
    */
   LocallyActiveDataPatchBoundaryEdgeSum(const string& object_name);

   /*!
    *  @brief Destructor for the schedule releases all internal storage.
    */
   ~LocallyActiveDataPatchBoundaryEdgeSum<DIM>();

   /*!
    *  @brief Register edge data with given patch data identifier for summing.
    *
    *  @param edge_data_id  integer patch data index for edge data to sum
    *
    *  The edge data id must be a valid patch data id (>=0) and must
    *  correspond to edge-centered double data.  If not, an error will result.
    */
   void registerSum(int edge_data_id);

   /*!
    *  @brief Set up summation operations for edge data across shared edges
    *         on a single level.
    *
    *  @param level         pointer to level on which to perform edge sum
    *  @param level_mgr     pointer to level mgr defining active patches
    *                       for edge data registered with the edge sum object
    *
    *  When assertion checking is active, neither pointer can be null.
    */
   void setupSum(
      tbox::Pointer<hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr);

   /*!
    *  @brief Compute sum of edge values at each shared edge and replace 
    *         each such edge value with the corresponding sum.  
    *
    *  At the end of this method, all values at shared edge locations on 
    *  patch boundaries will have the same value.  
    */
   void computeSum() const;

private:

   /*
    * Private member function to perform edge sum across single level --
    * called from computeSum()
    */
   void doLevelSum(
      tbox::Pointer<hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr) const;

   /*
    * Private member function to set up management of internal work data --
    * called from setupSum()
    */
   void setInternalWorkDataActive(
      tbox::Pointer<hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr);

   /*
    * Static members for managing shared temporary data among multiple
    * LocallyActiveDataPatchBoundaryEdgeSum objects.
    */
   static int s_instance_counter;
   // These arrays are indexed [data depth][number of variables with depth]
   static tbox::Array< tbox::Array<int> > s_oedge_src_id_array;
   static tbox::Array< tbox::Array<int> > s_oedge_dst_id_array;

   enum PATCH_BDRY_EDGE_SUM_DATA_ID { ID_UNDEFINED = -1 };

   string d_object_name;
   bool d_setup_called;

   int d_num_reg_sum;

   // These arrays are indexed [variable registration sequence number]
   tbox::Array<int> d_user_edge_data_id;
   tbox::Array<int> d_user_edge_depth;

   // These arrays are indexed [data depth]
   tbox::Array<int> d_num_registered_data_by_depth;

   /*
    * Edge-centered variables and patch data indices used as internal work quantities.
    */
   // These arrays are indexed [variable registration sequence number]
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_tmp_oedge_src_variable;
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_tmp_oedge_dst_variable;

   // These arrays are indexed [variable registration sequence number]
   tbox::Array<int> d_oedge_src_id;
   tbox::Array<int> d_oedge_dst_id;

   /*
    * Sets of indices for temporary variables to expedite allocation/deallocation.
    */
   hier::ComponentSelector d_oedge_src_data_set;
   hier::ComponentSelector d_oedge_dst_data_set;

   tbox::Pointer< hier::PatchLevel<DIM> > d_level;
   tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > d_level_mgr;

   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > 
      d_sum_transaction_factory;
   
   tbox::Pointer< xfer::LocallyActiveDataRefineSchedule<DIM> > 
      d_single_level_sum_schedule;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataPatchBoundaryEdgeSum.C"
#endif

