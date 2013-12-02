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
#include "SAMRAI/hier/SequentialLocalIdGenerator.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <set>

namespace SAMRAI {
namespace mesh {


/*
 * @brief Implementation of TransitLoad, representing the load with a
 * set of vouchers, each of which has a work value and the rank of the
 * process that issued the voucher.
 *
 * As a container, this class is nearly identical to
 * std::set<Voucher,VoucherMoreLoad>.
 *
 * Terminology: To avoid confusing the sending and receiving of
 * messages, the sending and receiving of work use the terms supply
 * and demand.  During redemption of vouchers, the holder of a voucher
 * demands the work.  The issuer of that voucher supplies the work.
 * Both send and receive messages to accomplish this.
 */

class VoucherTransitLoad : public TransitLoad {

public:

   typedef double LoadType;

   VoucherTransitLoad( const PartitioningParams &pparams );

   VoucherTransitLoad( const VoucherTransitLoad &other );

   //! @brief Clone object according to design pattern.
   boost::shared_ptr<TransitLoad> clone() const;

   //! @brief Initialize object.
   void initialize();

   //! @brief Return the total load contained.
   LoadType getSumLoad() const { return d_sumload; }

   //! @brief Insert all boxes from the given BoxContainer.
   void insertAll( const hier::BoxContainer &box_container );

   //! @brief Insert all boxes from the given TransitLoad.
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
    * @brief Assign unassigned boxes to local process and populate
    * balanced<==>unbalanced.
    *
    * This method uses communication to set up the map.
    */
   void
   assignContentToLocalProcessAndPopulateMaps(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector &balanced_to_unbalanced,
      hier::MappingConnector &unbalanced_to_balanced );


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


   //! @brief Voucher.
   struct Voucher {
      //! @brief Default constructor constructs voucher of zero value.
      Voucher() :
         d_issuer_rank(tbox::SAMRAI_MPI::getInvalidRank()),
         d_load(0.0) {
      }
      //! @brief Initializing constructor sets voucher value and issuer.
      Voucher( const LoadType &load, int issuer_rank ) :
         d_issuer_rank(issuer_rank),
         d_load(load) {
      }
      //! @brief Construct a voucher by combining two vouchers from the same issuer.
      Voucher( const Voucher &a, const Voucher &b ) :
         d_issuer_rank(a.d_issuer_rank),
         d_load( a.d_load + b.d_load ) {
         if ( a.d_issuer_rank != b.d_issuer_rank ) {
            TBOX_ERROR("VoucherTransitLoad: Cannot combine vouchers from different issuers.");
         }
      }
      //@ @brief Construct a Voucher by taking value from an existing Voucher.
      Voucher ( LoadType load, Voucher &src ) :
         d_issuer_rank(src.d_issuer_rank),
         d_load(load <= src.d_load ? load : src.d_load) {
         src.d_load -= d_load;
      }
      friend std::ostream &operator<<( std::ostream &co, const Voucher &v ) {
         co << v.d_issuer_rank << ':' << v.d_load;
         return co;
      }
      /*!
       * @brief Adjust load by taking work from or giving work to
       * another Voucher.
       *
       * Similar to the interface defined in TransitLoad but working
       * with individual vouchers instead of containers.
       */
      LoadType adjustLoad( Voucher &other, LoadType ideal_load );
      int d_issuer_rank;
      LoadType d_load;
   };


   //! @brief Comparison functor for sorting Vouchers by issuer rank and load.
   struct VoucherRankCompare {
      bool operator () ( const Voucher& a, const Voucher& b) const {
         if (a.d_issuer_rank != b.d_issuer_rank) {
            return a.d_issuer_rank < b.d_issuer_rank;
         }
         return a.d_load < b.d_load;
      }
   };


   //! @brief Encapsulates voucher redemption for both demander and supplier.
   struct VoucherRedemption {
      VoucherRedemption() :
         d_demander_rank(tbox::SAMRAI_MPI::getInvalidRank()),
         d_mpi(MPI_COMM_NULL),
         d_mpi_request(MPI_REQUEST_NULL) {};
      ~VoucherRedemption() {
         if ( d_mpi_request != MPI_REQUEST_NULL ) {
            TBOX_ERROR("Need to finish communication before destructing d_msg.");
            tbox::SAMRAI_MPI::Request_free(&d_mpi_request);
         }
      }

      //@{
      //! @name Demanding and supplying work for based on a voucher.
      void sendWorkDemand(
         const Voucher &v,
         const hier::SequentialLocalIdGenerator &id_gen,
         const tbox::SAMRAI_MPI &mpi );

      void recvWorkDemand(
         int demander_rank,
         int message_length,
         const tbox::SAMRAI_MPI &mpi );

      void sendWorkSupply(
         BoxTransitSet &reserve,
         hier::MappingConnector &unbalanced_to_balanced,
         hier::MappingConnector &balanced_to_unbalanced /* unused, can be removed */,
         const PartitioningParams &pparams );

      void recvWorkSupply(
         int message_length,
         hier::BoxLevel &balanced_box_level,
         hier::MappingConnector &unbalanced_to_balanced /* unused, can be removed */,
         hier::MappingConnector &balanced_to_unbalanced,
         const PartitioningParams &pparams );
      //@}

      Voucher d_voucher;
      int d_demander_rank;
      //! @brief Demander-specified LocalId generator to avoid ID clashes.
      hier::SequentialLocalIdGenerator d_id_gen;
      boost::shared_ptr<tbox::MessageStream> d_msg;
      tbox::SAMRAI_MPI d_mpi;
      MPI_Request d_mpi_request;
   };


   //@{
   //! @name Interfaces like the C++ standard stl::set, to help readability.
   typedef std::set<Voucher,VoucherRankCompare>::iterator iterator;
   typedef std::set<Voucher,VoucherRankCompare>::const_iterator const_iterator;
   typedef std::set<Voucher,VoucherRankCompare>::reverse_iterator reverse_iterator;
   typedef std::set<Voucher,VoucherRankCompare>::key_type key_type;
   typedef std::set<Voucher,VoucherRankCompare>::value_type value_type;
   iterator begin() { return d_voucher_set.begin(); }
   iterator end() { return d_voucher_set.end(); }
   const_iterator begin() const { return d_voucher_set.begin(); }
   const_iterator end() const { return d_voucher_set.end(); }
   reverse_iterator rbegin() { return d_voucher_set.rbegin(); }
   reverse_iterator rend() { return d_voucher_set.rend(); }
   size_t size() const { return d_voucher_set.size(); }
   bool empty() const { return d_voucher_set.empty(); }
   void clear() { d_sumload = 0; d_voucher_set.clear(); }
   std::pair<iterator, bool> insert( const Voucher &v ) {
      iterator itr = d_voucher_set.lower_bound( Voucher(0,v.d_issuer_rank) );
      if ( itr != d_voucher_set.end() &&
           itr->d_issuer_rank == v.d_issuer_rank ) {
         TBOX_ERROR("Cannot insert Voucher issued by rank " << v.d_issuer_rank
                    << " because there is already one from that issuer.\n"
                    << "To combine the vouchers, use insertCombine().");
      }
      itr = d_voucher_set.insert( itr, v );
      return std::pair<iterator, bool>(itr,true);
   }
   size_t erase( const Voucher &v ) {
      iterator vi = d_voucher_set.lower_bound(v);
      if ( vi != d_voucher_set.end() ) {
         d_sumload -= vi->d_load;
         erase(vi);
         return 1;
      }
      return 0;
   }
   void erase( iterator pos) {
      d_sumload -= pos->d_load;
      return d_voucher_set.erase(pos);
   }
   void swap( VoucherTransitLoad &other ) {
      const LoadType tmpload = d_sumload;
      d_sumload = other.d_sumload;
      other.d_sumload = tmpload;
      d_voucher_set.swap(other.d_voucher_set);
   }
   //@}


   /*!
    * @brief Insert voucher, combining with existing voucher from same issuer.
    */
   iterator insertCombine( const Voucher &v ) {
      iterator itr = d_voucher_set.lower_bound(v);
      if ( itr == d_voucher_set.end() ||
           v.d_issuer_rank != itr->d_issuer_rank ) {
         itr = d_voucher_set.insert( itr, v );
      }
      else {
         iterator pre_combine = itr;
         itr = d_voucher_set.insert( itr, Voucher(*itr,v) );
         d_voucher_set.erase(pre_combine);
      }
      d_sumload += v.d_load;
      return itr;
   }


   /*!
    * @brief Erase voucher issued by the given process.
    *
    * @return Whether there was a Voucher to be erased.
    */
   bool eraseIssuer( int issuer_rank );


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
    * @brief Raise load of dst container by shifing load from src.
    */
   static LoadType
   raiseDstLoad(
      VoucherTransitLoad& src,
      VoucherTransitLoad& dst,
      LoadType ideal_dst_load );


   /*!
    * @brief Re-cast a TransitLoad object to a VoucherTransitLoad.
    */
   const VoucherTransitLoad &recastTransitLoad( const TransitLoad &transit_load ) {
      const VoucherTransitLoad *ptr = static_cast<const VoucherTransitLoad*>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }


   /*!
    * @brief Re-cast a TransitLoad object to a VoucherTransitLoad.
    */
   VoucherTransitLoad &recastTransitLoad( TransitLoad &transit_load ) {
      VoucherTransitLoad *ptr = static_cast<VoucherTransitLoad*>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }


   /*!
    * @brief Find the voucher value issued by the given process.
    */
   LoadType findIssuerValue( int issuer_rank ) const;


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
   std::set<Voucher,VoucherRankCompare> d_voucher_set;

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
      boost::shared_ptr<tbox::Timer> t_assign_content_to_local_process_and_generate_map;
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
