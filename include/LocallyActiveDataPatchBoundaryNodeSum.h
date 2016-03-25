//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/femutils/locally_active/LocallyActiveDataPatchBoundaryNodeSum.h $
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Routines for summing locally-active node data at patch boundaries
//
 
#ifndef included_algs_LocallyActiveDataPatchBoundaryNodeSum
#define included_algs_LocallyActiveDataPatchBoundaryNodeSum

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_pdat_NodeVariable
#include "NodeVariable.h"
#endif
#ifndef included_pdat_OuternodeVariable
#include "OuternodeVariable.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_xfer_LocallyActiveDataCoarsenSchedule
#include "LocallyActiveDataCoarsenSchedule.h"
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
 *  @brief Class LocallyActiveDataPatchBoundaryNodeSum provides operations for summing 
 *  locally-active node data values at nodes shared by multiple patches on a single level 
 *  or across multiple hierarchy levels.  
 *
 *  NOTE: Currently, the functionality to sum data across multiple hierarchy levels is unavailable.
 *
 *  Usage of a patch boundry node sum involves the following sequence of steps:
 *      
 *  -# Construct a patch boundry node sum object.  For example,
 *     \verbatim 
 *         LocallyActiveDataPatchBoundaryNodeSum<DIM> my_node_sum("My Node Sum");
 *     \endverbatim 
 *  -# Register node data quantities to sum.  For example,
 *     \verbatim
 *         my_node_sum.registerSum(node_data_id1);
 *         my_node_sum.registerSum(node_data_id2);
 *         etc...
 *     \endverbatim
 *  -# Setup the sum operations for either single level 
 *     or a range of levels in a patch hierarchy.  
 *     For example,
 *     \verbatim
 *         my_node_sum.setupSum(level, level_mgr);    // single level
 *         -- or --
 *         my_node_sum.setupSum(hierarchy, coarsest_ln, finest_ln);  // multiple levels
 *     \endverbatim 
 *  -# Execute the sum operation.  For example,
 *     \verbatim
 *         my_node_sum.computeSum()
 *     \endverbatim
 *
 *  The result of these operations is that each node patch data value associated 
 *  with the registered ids at patch boundaries, on either the single level or 
 *  range of hierarchy levels, is replaced by the sum of all data values at the node.
 *
 *  Note that only one of the setupSum() functions may be called once a
 *  LocallyActiveDataPatchBoundaryNodeSum<DIM> object is created.
 */
 
template<int DIM> class LocallyActiveDataPatchBoundaryNodeSum
{
public:

   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots ahared among all LocallyActiveDataPatchBoundaryNodeSum
    *         objects (i.e., static members).  To get a correct count,
    *         this routine should only be called once.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with node sum objects.
    */
   static int getNumSharedPatchDataSlots(int max_variables_to_register);
 
   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots unique to each LocallyActiveDataPatchBoundaryNodeSum
    *         object (i.e., non-static members).  To get a correct count,
    *         this routine should only be called exactly once for each object
    *         that will be constructed.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with node sum objects.
    */
   static int getNumUniquePatchDataSlots(int max_variables_to_register);

   /*!
    *  @brief Constructor initializes object to default (mostly undefined) 
    *  state.
    *
    *  @param object_name const std::string reference for name of object used 
    *         in error reporting.  When assertion checking is on, the string 
    *         cannot be empty.
    */
   LocallyActiveDataPatchBoundaryNodeSum(const std::string& object_name);

   /*!
    *  @brief Destructor for the schedule releases all internal storage.
    */
   ~LocallyActiveDataPatchBoundaryNodeSum<DIM>();

   /*!
    *  @brief Register node data with given patch data identifier for summing.
    *
    *  @param node_data_id  integer patch data index for node data to sum
    *
    *  The node data id must be a valid patch data id (>=0) and must 
    *  correspond to node-centered double data.  If not, an error will result.
    */
   void registerSum(int node_data_id);

   /*!
    *  @brief Set up summation operations for node data across shared nodes
    *         on a single level.
    * 
    *  If the other setupSum() function for a range of hierarchy levels has 
    *  been called previously for this object, an error will result.
    * 
    *  @param level         pointer to level on which to perform node sum
    *  @param level_mgr     pointer to level mgr defining active patches
    *                       for node data registered with the node sum object
    *
    *  When assertion checking is active, the level pointer cannot be null. 
    */
   void setupSum(
      tbox::Pointer<hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr);

   /*!
    *  @brief Set up for summation operations for node data at shared nodes
    *         across a range of hierarchy levels.  The active patches for
    *         each node quantity on each level is accumed to be defined by 
    *         the corresponding hier::LocallyActiveDataPatchLevelManager<DIM>
    *         objects held by the hier::LocallyActiveVariableDatabase<DIM> 
    *         Singleton object.
    * 
    *  If the other setupSum() function for a single level has been called 
    *  previously for this object, an error will result.  
    *
    *  @param hierarchy      pointer to hierarchy on which to perform node sum
    *  @param coarsest_level coarsest level number for node sum
    *  @param finest_level   finest level number for node sum
    *
    *  When assertion checking is active, the hierarchy pointer cannot be null,
    *  and the range of levels must be valid.
    */
   void setupSum(tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy,
                 const int coarsest_level,
                 const int finest_level);

   /*!
    *  @brief Compute sum of node values at each shared node and replace 
    *         each such node value with the corresponding sum.  
    *
    *  At the end of this method, all values at shared node locations on
    *  patch boundaries (on levels indicated by the call to one of the 
    *  setupSum() routines) will have the same value.  
    *
    *  When the setupSum() method taking a range of patch levels in a
    *  hierarchy is called, this method will compute the sum of nodal 
    *  quantities at all the specified patch boundaries.  For nodes at a 
    *  coarse-fine boundary, nodal sums will only be performed where the 
    *  coarse and fine nodes overlap.  A node on a fine level that is not
    *  also a node on the next coarser level (a so-called "hanging node") 
    *  will not be summed.
    *
    *  The boolean "fill_hanging_nodes" argument specifies whether the
    *  the hanging nodes should be filled using linearly interpolated values 
    *  from neighboring non-hanging nodes (i.e. those that overlap nodes on 
    *  a coarse level). The correct steps required to deal with hanging
    *  nodes is algorithm dependent so, if left unspecified, values at the
    *  hanging nodes will not be adjusted.  However, because many algorithms 
    *  average hanging nodes we provide the capability to do it here.  Note that
    *  the hanging node interpolation provided does not take into consideration
    *  the spatial location of the nodes.  So the interpolation may not be 
    *  correct for coordinate systems other than standard Cartesian grid geometry.
    *
    *  @param fill_hanging_nodes Optional boolean value specifying whether hanging 
    *         node values should be set to values interpolated from neighboring 
    *         non-hanging node values.  The default is false.
    */
   void computeSum(const bool fill_hanging_nodes = false) const;

private:

   /*
    * Private member function to perform node sum across single level --
    * called from computeSum().
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
    * Private member function to set node node data on a fine level at a 
    * coarse-fine boundary to the sum of the node values and the associated
    * outernode values on a coarsened version of the fine level.
    */
   void doLocalCoarseFineBoundarySum(
      tbox::Pointer<hier::PatchLevel<DIM> > fine_level,
      tbox::Pointer<hier::PatchLevel<DIM> > coarsened_fine_level,
      const tbox::Array<int>& node_data_id,
      const tbox::Array<int>& onode_data_id,
      bool fill_hanging_nodes) const;

   /*
    * Private member function to copy node data to outernode data
    * on all patches on a level.
    */
   void copyNodeToOuternodeOnLevel(
      tbox::Pointer<hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
      const tbox::Array<int>& node_data_id,
      const tbox::Array<int>& onode_data_id) const;

   /*
    * Private member function to copy outernode data to node data
    * on all patches on a level.
    */
   void copyOuternodeToNodeOnLevel(
      tbox::Pointer<hier::PatchLevel<DIM> > level,
      tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > level_mgr,
      const tbox::Array<int>& onode_data_id,
      const tbox::Array<int>& node_data_id) const;

   /*
    * Static members for managing shared temporary data among multiple
    * LocallyActiveDataPatchBoundaryNodeSum objects.
    */
   static int s_instance_counter;
   // These arrays are indexed [data depth][number of variables with depth]
   static tbox::Array< tbox::Array<int> > s_onode_src_id_array;
   static tbox::Array< tbox::Array<int> > s_onode_dst_id_array;

   enum PATCH_BDRY_NODE_SUM_DATA_ID { ID_UNDEFINED = -1 };

   std::string d_object_name;
   bool d_setup_called;

   int d_num_reg_sum;

   // These arrays are indexed [variable registration sequence number]
   tbox::Array<int> d_user_node_data_id;
   tbox::Array<int> d_user_node_depth;

   // These arrays are indexed [data depth]
   tbox::Array<int> d_num_registered_data_by_depth;

   /*
    * Node-centered variables and patch data indices used as internal work quantities.
    */
   // These arrays are indexed [variable registration sequence number]
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_tmp_onode_src_variable;
   tbox::Array< tbox::Pointer< hier::Variable<DIM> > > d_tmp_onode_dst_variable;

   // These arrays are indexed [variable registration sequence number]
   tbox::Array<int> d_onode_src_id;
   tbox::Array<int> d_onode_dst_id;

   /*
    * Sets of indices for temporary variables to expedite allocation/deallocation.
    */
   hier::ComponentSelector d_onode_src_data_set;
   hier::ComponentSelector d_onode_dst_data_set;
   
   tbox::Pointer<hier::PatchLevel<DIM> > d_level;

   tbox::Pointer<hier::PatchHierarchy<DIM> > d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;

   tbox::Array< tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > > d_level_mgr;

   bool d_level_setup_called;
   bool d_hierarchy_setup_called;

   tbox::Pointer< xfer::LocallyActiveDataRefineTransactionFactory<DIM> > 
      d_sum_transaction_factory;

   tbox::Array< tbox::Pointer<xfer::LocallyActiveDataRefineSchedule<DIM> > >  
      d_single_level_sum_schedule;
   tbox::Array< tbox::Pointer<xfer::LocallyActiveDataRefineSchedule<DIM> > >  
      d_cfbdry_copy_schedule;
   tbox::Array< tbox::Pointer<xfer::LocallyActiveDataCoarsenSchedule<DIM> > > 
      d_sync_coarsen_schedule;

   tbox::Array< tbox::Pointer<hier::PatchLevel<DIM> > > d_cfbdry_tmp_level;

   tbox::Array< tbox::Pointer< hier::LocallyActiveDataPatchLevelManager<DIM> > > d_cfbdry_tmp_level_mgr;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataPatchBoundaryNodeSum.C"
#endif

