//
// File:	PatchBoundaryEdgeSum.h
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 696 $
// Modified:	$Date: 2005-11-03 12:27:01 -0800 (Thu, 03 Nov 2005) $
// Description:	Routines for summing edge data at patch boundaries
//
 
#ifndef included_algs_PatchBoundaryEdgeSum
#define included_algs_PatchBoundaryEdgeSum

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
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
#ifndef included_xfer_RefineSchedule
#include "RefineSchedule.h"
#endif
#ifndef included_xfer_RefineTransactionFactory
#include "RefineTransactionFactory.h"
#endif

namespace SAMRAI {
    namespace algs {

 /*!
 *  @brief Class PatchBoundaryEdgeSum provides operations summing edge data
 *  values at edges that are shared by multiple patches on a single level. 
 *  Note that this utility only works on a SINGLE patch level, not on a multiple 
 *  levels in an AMR patch hierarchy like the PatchBoundaryNodeSum class.   Unlike 
 *  node data, edge data at coarse-fine boundaries are not co-located, so the sum 
 *  operation is not clearly defined.
 *
 *  Usage of a patch boundry edge sum involves the following sequence of steps:
 *
 *  -# Construct a patch boundry edge sum object.  For example,
 *     \verbatim
 *         PatchBoundaryEdgeSum<DIM> my_edge_sum("My Edge Sum");
 *     \endverbatim
 *  -# Register edge data quantities to sum.  For example,
 *     \verbatim
 *         my_edge_sum.registerSum(edge_data_id1);
 *         my_edge_sum.registerSum(edge_data_id2);
 *         etc...
 *     \endverbatim
 *  -# Setup the sum operations for a single level.  For example,
 *     \verbatim
 *         my_edge_sum.setupSum(level);
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

template<int DIM> class PatchBoundaryEdgeSum
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
   PatchBoundaryEdgeSum(const string& object_name);

   /*!
    *  @brief Destructor for the schedule releases all internal storage.
    */
   ~PatchBoundaryEdgeSum<DIM>();

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
    *
    *  When assertion checking is active, the level pointer cannot be null.
    */
   void setupSum(tbox::Pointer<hier::PatchLevel<DIM> > level);

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
   void doLevelSum(tbox::Pointer<hier::PatchLevel<DIM> > level) const;

   /*
    * Static members for managing shared temporary data among multiple
    * PatchBoundaryEdgeSum objects.
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

   tbox::Pointer< xfer::RefineTransactionFactory<DIM> > d_sum_transaction_factory;
   
   tbox::Pointer< xfer::RefineSchedule<DIM> > d_single_level_sum_schedule;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchBoundaryEdgeSum.C"
#endif

