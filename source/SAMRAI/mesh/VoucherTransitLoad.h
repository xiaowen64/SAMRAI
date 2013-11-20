/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_VoucherTransitLoad
#define included_mesh_VoucherTransitLoad

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <set>

namespace SAMRAI {
namespace mesh {


/*
 * @brief Implementation of TreeLoadBalancer::TransitSet, representing
 * the load with a set of vouchers, sorted from highest load to lowest
 * load.
 *
 * As a container, this class is identical to
 * std::set<Voucher,VoucherMoreLoad>.
 */

class VoucherTransitLoad : public TransitLoad {

public:

   typedef double LoadType;

   VoucherTransitLoad( const PartitioningParams &pparams );

   //! @brief Return the total load contained.
   LoadType getSumLoad() const { return d_sumload; }

   //! @brief Insert all boxes from the given BoxContainer.
   void insertAll( const hier::BoxContainer &other );

   //! @brief Insert all boxes from the given VoucherTransitLoad.
   void insertAll( const TransitLoad &other );

   //! @brief Return number of items in this container.
   size_t getNumberOfItems() const;

   //! @brief Return number of processes contributing to the contents.
   size_t getNumberOfOriginatingProcesses() const;

   //@{
   //! @name Packing/unpacking for communication.

   void putToMessageStream( tbox::MessageStream &msg ) const;

   void getFromMessageStream( tbox::MessageStream &msg );
   //@}


   /*!
    * @brief Adjust the load in this VoucherTransitLoad by moving work
    * between it and another VoucherTransitLoad.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this VoucherTransitLoad.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoad(
      TransitLoad& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load );


   /*!
    * @brief Assign unassigned boxes to local process and generate
    * balanced<==>unbalanced map.
    *
    * This method uses communication to set up the map.
    */
   void
   assignContentToLocalProcessAndGenerateMap(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector &balanced_to_unbalanced,
      hier::MappingConnector &unbalanced_to_balanced ) const;


   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named "mesh::VoucherTransitLoad::*",
    * where the third field is the specific steps performed
    * by the Schedule.  You can override the first two
    * fields with this method.  Conforming to the timer
    * naming convention, timer_prefix should have the form
    * "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

   void recursivePrint(
      std::ostream &co=tbox::plog,
      const std::string &border=std::string(),
      int detail_depth=1 ) const;


private:

   static const int VoucherTransitLoad_EDGETAG0 = 3;
   static const int VoucherTransitLoad_EDGETAG1 = 4;


   /*!
    * @brief Voucher's load value.
    */
   struct VoucherValue {

      /*!
       * @brief Constructor
       */
      VoucherValue() : d_load(0) {}

      /*!
       * @brief Constructor
       */
      VoucherValue( double load ) : d_load(load) {}

      //! @brief Put self into a MessageStream.
      void putToMessageStream(tbox::MessageStream &mstream) const {
         mstream << d_load;
         return;
      }

      //! @brief Set attributes according to data in a MessageStream.
      void getFromMessageStream(tbox::MessageStream &mstream) {
         mstream >> d_load;
         return;
      }

      //! @brief Work load in this d_box.
      LoadType d_load;
   };


   //@{
   //! @name Interfaces like the C++ standard stl::map, to help readability.
   typedef std::map<int,LoadType>::iterator iterator;
   typedef std::map<int,LoadType>::const_iterator const_iterator;
   typedef std::map<int,LoadType>::reverse_iterator reverse_iterator;
   typedef std::map<int,LoadType>::key_type key_type;
   typedef std::map<int,LoadType>::value_type value_type;
   iterator begin() { return d_vouchers.begin(); }
   iterator end() { return d_vouchers.end(); }
   const_iterator begin() const { return d_vouchers.begin(); }
   const_iterator end() const { return d_vouchers.end(); }
   size_t size() const { return d_vouchers.size(); }
   bool empty() const { return d_vouchers.empty(); }
   void clear() { d_sumload = 0; d_vouchers.clear(); }
   void swap( VoucherTransitLoad &other ) {
      const LoadType tl = d_sumload;
      d_sumload = other.d_sumload;
      other.d_sumload = tl;
      d_vouchers.swap(other.d_vouchers);
   }
   //@}


   /*!
    * @brief Insert Voucher into an output stream.
    */
   friend std::ostream& operator << ( std::ostream& co,
                                      const VoucherValue& r ) {
      co << r.d_load;
      return co;
   }


   /*!
    * @brief Construct semilocal relationships in
    * unbalanced--->balanced Connector.
    *
    * Constructing semilocal unbalanced--->balanced relationships
    * require communication to determine where exported work ended up.
    *
    * @param [out] unbalanced_to_balanced Connector to store
    * relationships in.
    *
    * @param [in] kept_imports Work that was imported and locally kept.
    */
   void constructSemilocalUnbalancedToBalanced(
      hier::MappingConnector &unbalanced_to_balanced,
      const VoucherTransitLoad &kept_imports ) const;


   /*!
    * @brief Re-cast a TransitLoad object to a BoxTransitSet.
    */
   const VoucherTransitLoad &recastTransitLoad( const TransitLoad &transit_load ) {
      const VoucherTransitLoad *ptr = static_cast<const VoucherTransitLoad*>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }


   /*!
    * @brief Re-cast a TransitLoad object to a BoxTransitSet.
    */
   VoucherTransitLoad &recastTransitLoad( TransitLoad &transit_load ) {
      VoucherTransitLoad *ptr = static_cast<VoucherTransitLoad*>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }


   //! @brief Compute the load for a Box.
   LoadType computeLoad( const hier::Box& box) const {
      return LoadType(box.size());
   }


   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void initializeCallback() {
      TimerStruct& timers(s_static_timers[s_default_timer_prefix]);
      getAllTimers(s_default_timer_prefix, timers);
   }


   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void finalizeCallback() {
      s_static_timers.clear();
   }


   //! @brief Balance penalty is proportional to imbalance.
   double computeBalancePenalty( double imbalance) const {
      return tbox::MathUtilities<double>::Abs(imbalance);
   }


   //! @brief Work load, sorted by originating rank.
   std::map<int,LoadType> d_vouchers;

   LoadType d_sumload;

   const PartitioningParams *d_pparams;


   //@{
   //! @name Debugging stuff.
   bool d_print_steps;
   bool d_print_edge_steps;
   //@}

   //@{
   //! @name Timer data for Schedule class.

   /*
    * @brief Structure of timers used by this class.
    *
    * Each Schedule object can set its own timer names through
    * setTimerPrefix().  This leads to many timer look-ups.  Because
    * it is expensive to look up timers, this class caches the timers
    * that has been looked up.  Each TimerStruct stores the timers
    * corresponding to a prefix.
    */
   struct TimerStruct {
      boost::shared_ptr<tbox::Timer> t_adjust_load;
      boost::shared_ptr<tbox::Timer> t_adjust_load_by_popping;
      boost::shared_ptr<tbox::Timer> t_adjust_load_by_swapping;
      boost::shared_ptr<tbox::Timer> t_shift_loads_by_breaking;
      boost::shared_ptr<tbox::Timer> t_find_swap_pair;
      boost::shared_ptr<tbox::Timer> t_assign_content_to_local_process_and_generate_map;
      boost::shared_ptr<tbox::Timer> t_construct_semilocal;
      boost::shared_ptr<tbox::Timer> t_construct_semilocal_comm_wait;
      boost::shared_ptr<tbox::Timer> t_construct_semilocal_send_edges;
      boost::shared_ptr<tbox::Timer> t_pack_edge;
      boost::shared_ptr<tbox::Timer> t_unpack_edge;
      boost::shared_ptr<tbox::Timer> t_post_load_distribution_barrier;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

   /*!
    * @brief Structure of timers in s_static_timers, matching this
    * object's timer prefix.
    */
   TimerStruct* d_object_timers;

   /*!
    * @brief Get all the timers defined in TimerStruct.  The timers
    * are named with the given prefix.
    */
   static void
   getAllTimers(
      const std::string& timer_prefix,
      TimerStruct& timers);

   //@}

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;
};



}
}

#endif
