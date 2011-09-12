/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A class describing the adjacency of Boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxNeighborhoodCollection_C
#define included_hier_BoxNeighborhoodCollection_C

#include "SAMRAI/hier/BoxNeighborhoodCollection.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxSetSingleBlockIterator.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxNeighborhoodCollection.I"
#endif

namespace SAMRAI {
namespace hier {

BoxNeighborhoodCollection::BoxNeighborhoodCollection(
   int rank,
   const BoxContainer& roots) :
   d_rank(rank)
{
   // For each root in roots create an empty neighborhood.
   for (BoxContainer::ConstIterator itr(roots.begin()); itr; ++itr) {
      d_roots.insert(&(d_adj_list.insert(std::make_pair(*itr, Neighbors())).first->first));
   }
}

BoxNeighborhoodCollection::BoxNeighborhoodCollection(
   const BoxNeighborhoodCollection& other)
{
   // Iterate through the other collection and create in this the same
   // neighborhoods that the other contains.
   for (ConstIterator root_itr(other.begin());
        root_itr != other.end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(other.begin(root_itr));
           nbr_itr != other.end(root_itr); ++nbr_itr) {
         insert(*root_itr, *nbr_itr);
      }
   }
}

BoxNeighborhoodCollection& BoxNeighborhoodCollection::operator = (
   const BoxNeighborhoodCollection& rhs)
{
   // Empty this container then iterate through the other collection and
   // create in this the same neighborhoods that the other contains.
   clear();
   for (ConstIterator root_itr(rhs.begin());
        root_itr != rhs.end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(rhs.begin(root_itr));
           nbr_itr != rhs.end(root_itr); ++nbr_itr) {
         insert(*root_itr, *nbr_itr);
      }
   }
   return *this;
}

bool BoxNeighborhoodCollection::operator == (
   const BoxNeighborhoodCollection& rhs) const
{
   // Both collections must have the same number of box neighborhoods.
   bool result = numBoxNeighborhoods() == rhs.numBoxNeighborhoods();

   if (result) {
      // Both collections have the same number of box neighborhoods so now
      // check that the neighborhoods are the same.
      ConstIterator root_itr(begin());
      ConstIterator rhs_root_itr(rhs.begin());
      for (; root_itr != end(); ++root_itr, ++rhs_root_itr) {
         // Check that this root is the same in each collection and that the
         // number of neighbors of this root is the same in each collection.
         if (!(*root_itr).isSpatiallyEqual(*rhs_root_itr) ||
             !(*root_itr).isIdEqual(*rhs_root_itr)) {
            result = false;
            break;
         }
         if (numNeighbors(root_itr) != rhs.numNeighbors(rhs_root_itr)) {
            result = false;
            break;
         }

         // This root and its number of neighbors match so check that each
         // neighbor from each collection matches.
         ConstNeighborIterator nbr_itr(begin(root_itr));
         ConstNeighborIterator rhs_nbr_itr(rhs.begin(rhs_root_itr));
         for (; nbr_itr != end(root_itr); ++nbr_itr, ++rhs_nbr_itr) {
            if (!(*nbr_itr).isSpatiallyEqual(*rhs_nbr_itr) ||
                !(*nbr_itr).isIdEqual(*rhs_nbr_itr)) {
               result = false;
               break;
            }
         }

         // If one of the neighbors of this root didn't match bail out.
         if (!result) {
            break;
         }
      }
   }

   return result;
}

BoxNeighborhoodCollection::ConstIterator BoxNeighborhoodCollection::find(
   const Box& root) const
{
   // Look for "root" in the set of roots.  If it's not there then it's not a
   // root and return an iterator pointing to the end.  Otherwise return an
   // iterator pointing to it.
   RootsItr root_itr(d_roots.find(&root));
   if (root_itr == d_roots.end()) {
      ConstIterator result(end());
      return result;
   }
   else {
      ConstIterator result(*this, d_adj_list.find(root));
      return result;
   }
}

BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::find(
   const Box& root)
{
   // Look for "root" in the set of roots.  If it's not there then it's not a
   // root and return an iterator pointing to the end.  Otherwise return an
   // iterator pointing to it.
   RootsItr root_itr(d_roots.find(&root));
   if (root_itr == d_roots.end()) {
      Iterator result(end());
      return result;
   }
   else {
      Iterator result(*this, d_adj_list.find(root));
      return result;
   }
}

int BoxNeighborhoodCollection::numBoxNeighborhoods() const
{
   // Count the number of roots.
   int ct = 0;
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      ++ct;
   }
   return ct;
}

int BoxNeighborhoodCollection::totalNumNeighbors() const
{
   // Count the neighbors in each root.
   int ct = 0;
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      ct += numNeighbors(root_itr);
   }
   return ct;
}

bool BoxNeighborhoodCollection::hasPeriodicNeighborhood(
   const ConstIterator& root_itr) const
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::hasPeriodicNeighborhood root_itr points to end");
   }

   // See if any neighbor of the supplied root is periodic.
   bool result = false;
   for (ConstNeighborIterator nbr_itr(begin(root_itr));
        nbr_itr != end(root_itr); ++nbr_itr) {
      if ((*nbr_itr).isPeriodicImage()) {
         result = true;
         break;
      }
   }
   return result;
}

bool BoxNeighborhoodCollection::isLocal() const
{
   // See if all the neighbors of all roots belong to the same processor as
   // this object.
   bool is_local = true;
   for (ConstIterator root_itr(begin()); is_local && root_itr != end();
        ++root_itr) {
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++nbr_itr) {
         if ((*nbr_itr).getOwnerRank() != d_rank) {
            is_local = false;
            break;
         }
      }
   }
   return is_local;
}

void BoxNeighborhoodCollection::getOwners(
   std::set<int>& owners) const
{
   // Put the processor id of all neigbhors of all roots into owners.
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++nbr_itr) {
         owners.insert((*nbr_itr).getOwnerRank());
      }
   }
   return;
}

BoxNeighborhoodCollection::InsertRetType  BoxNeighborhoodCollection::insert(
   const Box& root,
   const Box& new_nbr)
{
   // First, add the root to the collection.  If it's already a root this is a
   // no_op.
   InsertRetType result(insert(root));

   // Now call insert with the iterator pointing to the root.
   insert(result.first, new_nbr);

   return result;
}

void BoxNeighborhoodCollection::insert(
   Iterator& root_itr,
   const Box& new_nbr)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::insert root_itr points to end");
   }

   // First add the new_nbr to the collection if it is not there.
   AdjListItr nbr_itr(d_adj_list.find(new_nbr));
   if (nbr_itr == d_adj_list.end()) {
      nbr_itr = d_adj_list.insert(
         nbr_itr,
         std::make_pair(new_nbr, Neighbors()));
   }

   // Now add the neighbor to the root's neighborhood.
   root_itr.d_itr->second.insert(&(nbr_itr->first));

   return;
}

BoxNeighborhoodCollection::InsertRetType BoxNeighborhoodCollection::insert(
   const Box& root,
   const BoxContainer& new_nbrs)
{
   // First, add the root to the collection.  If it's already a root this is a
   // no_op.
   InsertRetType result(insert(root));

   // Now call insert with the iterator pointing to the root.
   insert(result.first, new_nbrs);

   return result;
}

void BoxNeighborhoodCollection::insert(
   Iterator& root_itr,
   const BoxContainer& new_nbrs)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::insert root_itr points to end");
   }

   // Add each neighbor in the container to the root.
   for (BoxContainer::ConstIterator nbr_itr(new_nbrs.begin());
        nbr_itr; ++nbr_itr) {

      // First add this new neighbor to the collection if it is not there.
      AdjListItr itr(d_adj_list.find(*nbr_itr));
      if (itr == d_adj_list.end()) {
         itr = d_adj_list.insert(itr, std::make_pair(*nbr_itr, Neighbors()));
      }

      // Now add this new neighbor to the root's neighborhood.
      root_itr.d_itr->second.insert(&itr->first);
   }

   return;
}

void BoxNeighborhoodCollection::erase(
   Iterator& root_itr,
   const BoxContainer& nbrs)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::erase root_itr points to end");
   }

   // Remove each neighbor in the container from the root.
   for (BoxContainer::ConstIterator nbr_itr(nbrs.begin());
        nbr_itr; ++nbr_itr) {

      // Remove this neighbor from the root's neighborhood.
      root_itr.d_itr->second.erase(&(*nbr_itr));
   }
   return;
}

BoxNeighborhoodCollection::InsertRetType BoxNeighborhoodCollection::insert(
   const Box& new_root)
{
   bool added_root;

   // First, add the root to the collection if it is not there.
   AdjListItr root_itr(d_adj_list.find(new_root));
   if (root_itr == d_adj_list.end()) {
      root_itr = d_adj_list.insert(
         root_itr,
         std::make_pair(new_root, Neighbors()));
      added_root = true;
   }
   else {
      added_root = false;
   }

   // Second, add the root to the set of roots.  If it's already there this is
   // a no-op.
   d_roots.insert(&root_itr->first);

   return std::make_pair(Iterator(*this, root_itr), added_root);
}

void BoxNeighborhoodCollection::erase(
   Iterator& root_itr,
   bool erase_root)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::erase root_itr points to end");
   }

   if (erase_root) {
      // Erasing root so first remove it from neighborhoods of other roots.
      for (Iterator other_roots(begin());
           other_roots != end(); ++other_roots) {
         for (NeighborIterator other_nbrs(begin(other_roots));
              other_nbrs != end(other_roots); ++other_nbrs) {
            if ((*other_nbrs).isIdEqual(*root_itr)) {
               other_roots.d_itr->second.erase(other_nbrs.d_itr);
               break;
            }
         }
      }

      // Now remove the root.
      d_roots.erase(&(*root_itr));
      d_adj_list.erase(root_itr.d_itr);
   }
   else {
      // Not erasing root so just make its neighborhood empty.
      root_itr.d_itr->second.clear();
   }
   return;
}

void BoxNeighborhoodCollection::erase(
   Iterator& first_root_itr,
   Iterator& last_root_itr,
   bool erase_roots)
{
   // For each root in the range erase it.
   for (Iterator root_itr(first_root_itr); root_itr != last_root_itr;) {
      Iterator current_root(root_itr);
      ++root_itr;
      erase(current_root, erase_roots);
   }
   return;
}

void BoxNeighborhoodCollection::localize()
{
   // Find all roots which do not belong to the same processor as this object
   // and remove them entirely.
   for (Iterator root_itr(begin()); root_itr != end();) {
      if ((*root_itr).getOwnerRank() != d_rank) {
         Iterator current_root(root_itr);
         ++root_itr;
         erase(current_root, true);
      }
      else {
         ++root_itr;
      }
   }
   return;
}

void BoxNeighborhoodCollection::eraseEmptyNeighborhoods()
{
   // Find all roots which have no neighbors and remove them entirely.
   for (Iterator root_itr(begin()); root_itr != end();) {
      if (root_itr.d_itr->second.size() == 0) {
         Iterator current_root(root_itr);
         ++root_itr;
         erase(current_root, true);
      }
      else {
         ++root_itr;
      }
   }
   return;
}

void BoxNeighborhoodCollection::clear()
{
   // Entirely remove all roots.
   for (Iterator root_itr(begin()); root_itr != end();) {
      Iterator current_root(root_itr);
      ++root_itr;
      erase(current_root, true);
   }
   return;
}

void BoxNeighborhoodCollection::coarsenNeighbors(
   BoxNeighborhoodCollection& result,
   const IntVector& ratio) const
{
   if (!result.empty()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::coarsenNeighbors result not empty");
   }

   // First, coarsen each neighbor and put the coarsened Boxes into
   // coarsened_nbrs.
   Box nbr(ratio.getDim());
   std::set<Box, box_less> coarsened_nbrs;
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++root_itr) {
         if (coarsened_nbrs.find(*nbr_itr) != coarsened_nbrs.end()) {
            nbr = *nbr_itr;
            nbr.coarsen(ratio);
            coarsened_nbrs.insert(nbr);
         }
      }
   }

   // Next, replicate this collection of neighborhoods in the result using the
   // Boxes in coarsened_nbrs.  Look for roots in coarsened_nbrs as a root may
   // be a neighbor of another root and we want the coarsened entity in the
   // result.
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      std::set<Box, box_less>::iterator root_as_nbr(coarsened_nbrs.find(*root_itr));
      bool root_is_nbr = root_as_nbr != coarsened_nbrs.end();
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++root_itr) {
         if (root_is_nbr) {
            result.insert(*root_as_nbr, *(coarsened_nbrs.find(*nbr_itr)));
         }
         else {
            result.insert(*root_itr, *(coarsened_nbrs.find(*nbr_itr)));
         }
      }
   }
   return;
}

void BoxNeighborhoodCollection::refineNeighbors(
   BoxNeighborhoodCollection& result,
   const IntVector& ratio) const
{
   if (!result.empty()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::refineNeighbors result not empty");
   }

   // First, refine each neighbor and put the refined Boxes into refined_nbrs.
   Box nbr(ratio.getDim());
   std::set<Box, box_less> refined_nbrs;
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++root_itr) {
         if (refined_nbrs.find(*nbr_itr) != refined_nbrs.end()) {
            nbr = *nbr_itr;
            nbr.refine(ratio);
            refined_nbrs.insert(nbr);
         }
      }
   }

   // Next, replicate this collection of neighborhoods in the result using the
   // Boxes in refined_nbrs.  Look for roots in refined_nbrs as a root may
   // be a neighbor of another root and we want the refined entity in the
   // result.
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      std::set<Box, box_less>::iterator root_as_nbr(refined_nbrs.find(*root_itr));
      bool root_is_nbr = root_as_nbr != refined_nbrs.end();
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++root_itr) {
         if (root_is_nbr) {
            result.insert(*root_as_nbr, *(refined_nbrs.find(*nbr_itr)));
         }
         else {
            result.insert(*root_itr, *(refined_nbrs.find(*nbr_itr)));
         }
      }
   }
   return;
}

void BoxNeighborhoodCollection::growNeighbors(
   BoxNeighborhoodCollection& result,
   const IntVector& growth) const
{
   if (!result.empty()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::growNeighbors result not empty");
   }

   // First, grow each neighbor and put the grown Boxes into grown_nbrs.
   Box nbr(growth.getDim());
   std::set<Box, box_less> grown_nbrs;
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++root_itr) {
         if (grown_nbrs.find(*nbr_itr) != grown_nbrs.end()) {
            nbr = *nbr_itr;
            nbr.grow(growth);
            grown_nbrs.insert(nbr);
         }
      }
   }

   // Next, replicate this collection of neighborhoods in the result using the
   // Boxes in grown_nbrs.  Look for roots in grown_nbrs as a root may be a
   // neighbor of another root and we want the grown entity in the result.
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      std::set<Box, box_less>::iterator root_as_nbr(grown_nbrs.find(*root_itr));
      bool root_is_nbr = root_as_nbr != grown_nbrs.end();
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++root_itr) {
         if (root_is_nbr) {
            result.insert(*root_as_nbr, *(grown_nbrs.find(*nbr_itr)));
         }
         else {
            result.insert(*root_itr, *(grown_nbrs.find(*nbr_itr)));
         }
      }
   }
   return;
}

void BoxNeighborhoodCollection::getNeighbors(
   BoxSet& neighbors) const
{
   // Iterate through the neighbors of each neighborhood and dump them into the
   // BoxSet.
   for (AdjListConstItr root_itr(d_adj_list.begin());
        root_itr != d_adj_list.end(); ++root_itr) {
      for (NeighborsConstItr nbr_itr(root_itr->second.begin());
           nbr_itr != root_itr->second.end(); ++nbr_itr) {
         neighbors.insert(*(*nbr_itr));
      }
   }
   return;
}

void BoxNeighborhoodCollection::getNeighbors(
   BoxList& neighbors) const
{
   // This is a 2 step process so it's more expensive than the version using a
   // BoxSet for the result.  This should go away when the unfied BoxContainer
   // is in place.

   // First call the version taking a BoxSet so that there are no duplicates.
   // Then dump the neighbors into the supplied BoxList.
   BoxSet tmp_nbrs;
   getNeighbors(tmp_nbrs);
   for (BoxSet::const_iterator itr(tmp_nbrs.begin());
        itr != tmp_nbrs.end(); ++itr) {
      neighbors.appendItem(*itr);
   }
   return;
}

void BoxNeighborhoodCollection::getNeighbors(
   BoxList& neighbors,
   const BlockId& block_id) const
{
   // This is a 2 step process so it's more expensive than the version using a
   // BoxSet for the result.  This should be changed to take a BoxSet when the
   // unfied BoxContainer is in place.

   // First call the version taking a BoxSet so that there are no duplicates.
   // Then dump the neighbors with the desired BlockId into the supplied
   // BoxList.
   BoxSet tmp_nbrs;
   getNeighbors(tmp_nbrs);
   for (BoxSetSingleBlockIterator itr(tmp_nbrs, block_id);
        itr.isValid(); ++itr) {
      neighbors.appendItem(*itr);
   }
   return;
}

void BoxNeighborhoodCollection::getNeighbors(
   std::map<BlockId, BoxList>& neighbors) const
{
   // This is a 2 step process so it's more expensive than the version using a
   // BoxSet for the result.  This should be changed to take a BoxSet when the
   // unfied BoxContainer is in place.

   // First call the version taking a BoxSet so that there are no duplicates.
   // Then dump each neighbor into the BoxList associated with its BlockId.
   BoxSet tmp_nbrs;
   getNeighbors(tmp_nbrs);
   for (BoxSet::const_iterator itr(tmp_nbrs.begin());
        itr != tmp_nbrs.end(); ++itr) {
      neighbors[itr->getBlockId()].appendItem(*itr);
   }
   return;
}

void BoxNeighborhoodCollection::getPeriodicNeighbors(
   BoxSet& result) const
{
   // Iterate through each neighbor of each root and place each neighbor which
   // is a periodic image into the supplied BoxSet.
   for (ConstIterator root_itr(begin()); root_itr != end(); ++root_itr) {
      for (ConstNeighborIterator nbr_itr(begin(root_itr));
           nbr_itr != end(root_itr); ++nbr_itr) {
         if ((*nbr_itr).isPeriodicImage()) {
            result.insert(*nbr_itr);
         }
      }
   }
   return;
}

void BoxNeighborhoodCollection::putToIntBuffer(
   std::vector<int>& send_mesg,
   const tbox::Dimension& dim,
   int offset,
   int buff_init) const
{
   (void)offset;

   /*
    * Pack neighborhood info into send_mesg starting at the offset location.
    */
   /*
    * Information to be packed:
    *   - Number of local boxes with neighbors (no info to send
    *     for those boxes without neighbors)
    *   - For each local box,
    *     - box's local index
    *     - number of neighbors
    *     - neighbors
    */
   int box_com_buf_size = Box::commBufferSize(dim);
   const int num_box_neighborhoods = numBoxNeighborhoods();
   int num_nabrs = totalNumNeighbors();

   /* Message size is 1 for number of local roots with neighbor lists plus
    * the number of neighbors for each root plus the amount of box data
    * written for each root and each root's neighbors. */
   const int mesg_size = 1 + num_box_neighborhoods +
      (num_box_neighborhoods + num_nabrs) * box_com_buf_size;

   send_mesg.resize(mesg_size, buff_init);
   send_mesg[0] = num_box_neighborhoods;
   int imesg = 1;

   for (ConstIterator ci(begin()); ci != end(); ++ci) {

      /* Info for each root. */
      (*ci).putToIntBuffer(&send_mesg[imesg]);
      imesg += box_com_buf_size;

      /* Number of neighbors for this root. */
      send_mesg[imesg++] = numNeighbors(ci);

      for (ConstNeighborIterator ni(begin(ci));
           ni != end(ci); ++ni) {
         /* Info for each neighbor. */
         (*ni).putToIntBuffer(&send_mesg[imesg]);
         imesg += box_com_buf_size;
      }

   }

   TBOX_ASSERT(imesg == mesg_size);

   return;
}

void BoxNeighborhoodCollection::getFromIntBuffer(
   const std::vector<int>& recv_mesg,
   const std::vector<int>& proc_offset,
   const tbox::Dimension& dim,
   int num_proc)
{
   /*
    * Unpack neighborhood info from recv_mesg.
    */
   int box_com_buf_size = Box::commBufferSize(dim);

   Box root(dim);
   Box nabr(dim);
   for (int n = 0; n < num_proc; ++n) {

      if (n != d_rank) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int num_nabr_lists = *(ptr++);

         for (int i = 0; i < num_nabr_lists; ++i) {

            root.getFromIntBuffer(ptr);
            ptr += box_com_buf_size;

            const int num_nabrs = (*ptr++);

            for (int nn = 0; nn < num_nabrs; ++nn) {
               nabr.getFromIntBuffer(ptr);
               ptr += box_com_buf_size;
               insert(root, nabr);
            }

         }

      }
   }
   return;
}

// These are defined in NeighborhoodSet but it's not clear that they
// are ever called or even needed.
void BoxNeighborhoodCollection::putToDatabase(
   tbox::Database& database) const
{
   return;
}

void BoxNeighborhoodCollection::getFromDatabase(
   const tbox::Database& database)
{
   return;
}

BoxNeighborhoodCollection::Iterator::Iterator(
   BoxNeighborhoodCollection& nbrhds,
   bool from_start) :
   d_collection(&nbrhds),
   d_itr(from_start ? nbrhds.d_adj_list.begin() :
                      nbrhds.d_adj_list.end())
{
   // If the start of the iteration is requested then find the first member in
   // the map of roots and neighbors which corresponds to the first root.
   if (from_start) {
      const Box* first_root = *(nbrhds.d_roots.begin());
      for (; d_itr != nbrhds.d_adj_list.end() && !d_itr->first.isIdEqual(*first_root);
           ++d_itr) {
      }
   }
}

BoxNeighborhoodCollection::Iterator&
BoxNeighborhoodCollection::Iterator::operator ++ ()
{
   // Find which box is the next root.
   RootsItr roots(d_collection->d_roots.find(&(d_itr->first)));
   ++roots;

   // If there are no more roots then move the iterator to the end of the map
   // of roots and neighbors.  Otherwise, move the iterator to the next root in
   // the map of roots and neighbors.
   if (roots == d_collection->d_roots.end()) {
      for (; d_itr != d_collection->d_adj_list.end(); ++d_itr) {
      }
   }
   else {
      const Box* next_root = *roots;
      ++d_itr;
      for (; d_itr != d_collection->d_adj_list.end() && !d_itr->first.isIdEqual(*next_root);
           ++d_itr) {
      }
   }
   return *this;
}

BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const BoxNeighborhoodCollection& nbrhds,
   bool from_start) :
   d_collection(&nbrhds),
   d_itr(from_start ? nbrhds.d_adj_list.begin() :
                      nbrhds.d_adj_list.end())
{
   // If the start of the iteration is requested then find the first member in
   // the map of roots and neighbors which corresponds to the first root.
   if (from_start) {
      const Box* first_root = *(nbrhds.d_roots.begin());
      for (; d_itr != nbrhds.d_adj_list.end() && !d_itr->first.isIdEqual(*first_root);
           ++d_itr) {
      }
   }
}

BoxNeighborhoodCollection::ConstIterator&
BoxNeighborhoodCollection::ConstIterator::operator ++ ()
{
   // Find which box is the next root.
   RootsItr roots(d_collection->d_roots.find(&(d_itr->first)));
   ++roots;

   // If there are no more roots then move the iterator to the end of the map
   // of roots and neighbors.  Otherwise, move the iterator to the next root in
   // the map of roots and neighbors.
   if (roots == d_collection->d_roots.end()) {
      for (; d_itr != d_collection->d_adj_list.end(); ++d_itr) {
      }
   }
   else {
      const Box* next_root = *roots;
      ++d_itr;
      for (; d_itr != d_collection->d_adj_list.end() && d_itr->second.size() == 0;
           ++d_itr) {
      }
   }
   return *this;
}

}
}

#endif
