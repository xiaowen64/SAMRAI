//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/multiblock/MultiblockPatchHierarchy.h $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: class to manage multiblocks
//

#ifndef included_hier_MultiblockPatchHierarchy
#define included_hier_MultiblockPatchHierarchy

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "BoxArray.h"
#include "BasePatchHierarchy.h"
#include "PatchHierarchy.h"
#include "MultiblockGridGeometry.h"
#include "MultiblockPatchLevel.h"
#include "tbox/Array.h"
#include "tbox/List.h"


namespace SAMRAI {
    namespace hier {

/*!
 * @brief Class MultiblockPatchHierarchy<DIM> manages an array of patch
 * hierarchies that represent a multiblock domain, and describes the
 * relationship between these hierarchies
 *
 * Each patch hierarchy contained in the array of patch hierarchies
 * represents a logically rectangular block of a multiblock domain.  This
 * class contains this array, and also contains information that describes
 * the relationship between neighboring blocks.  It also contains
 * information about the point or points of singularity that exist in the
 * multiblock domain.
 *
 * @see hier::BasePatchHierarchy
 * @see hier::PatchHierarchy
 * @see hier::MultiblockPatchLevel
 * @see hier::PatchLevel
 */

template<int DIM>
class MultiblockPatchHierarchy
: public hier::BasePatchHierarchy<DIM>
{

public:

   /*!
    * @brief Constructor for MultiblockPatchHierarchy.
    *
    * @param object_name String identifier for database operations
    * @param input_db    Input Database
    * @param geometry    Geometry
    * @param register_for_restart  Boolean switch for restart registration
    */
   MultiblockPatchHierarchy(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db,
      tbox::Pointer< hier::MultiblockGridGeometry<DIM> >& geometry,
      const bool register_for_restart = true);

   /*!
    * @brief Destructor for multiblock
    */
   ~MultiblockPatchHierarchy<DIM>();

   /*!
    * @brief Get a pointer to a single patch hierarchy represented by
    * the specified block number
    *
    * @param block_num Integer block number, corresponding to array index
    *                  of the hiearchy array passed into the constructor
    */
   tbox::Pointer< hier::PatchHierarchy<DIM> >& getHierarchy(
      const int block_num);

   /*!
    * @brief Get a pointer to a MultiblockPatchLevel<DIM> associated with the
    * given level number.
    *
    * The MultiblockPatchLevel<DIM> will give access to an array of
    * pointers to hier::PatchLevel<DIM>, each member of which represents a
    * level from one block of the hierarchy corresponding to the specified
    * level number.
    *
    * @param level_num Level number in AMR hierachy of desired multiblock
    *                  patch level.
    */ 
   tbox::Pointer< MultiblockPatchLevel<DIM> >
   getMultiblockPatchLevel(const int level_num) const;

   /*!
    * @brief Get a pointer to a hier::BasePatchLevel<DIM> association with the
    * given level number.
    *
    * Returns a pointer to an object of the virtual base class
    * hier::BasePatchLevel<DIM>.  The object being pointed to can be cast to
    * MultiblockPatchLevel<DIM>.
    *
    * @param level_num Level number in AMR hierachy of desired multiblock
    *                  patch level.
    */
   tbox::Pointer< hier::BasePatchLevel<DIM> > getPatchLevel(
      const int level_num) const;

   /*!
    * @brief Adjust boundary data of a level to be consistent with the
    * multiblock nature of the domain.
    *
    * On a MultiblockPatchHierarchy, each hier::PatchHierarchy<DIM> will
    * contain PatchLevels with patches that were constructed independent of
    * any knowledge of the multiblock nature of the complete domain.  Thus
    * the patches will contain boundary data that recognizes no difference
    * between a physical domain boundary and a block boundary.  Calling this
    * routine will adjust the boundary data on all patches in the given level
    * such that the true boundaries of the domain are represented.
    *
    * @param level Multiblock level where boundaries need to be adjusted.
    */ 
   void adjustMultiblockPatchLevelBoundaries(
      tbox::Pointer< MultiblockPatchLevel<DIM> > level);

   /*!
    * @brief Returns the number of blocks in the multiblock domain.
    */
   int getNumberOfBlocks() const;

   /*!
    * @brief Return the number of neighbors a specific block of the Multiblock
    * domain has.
    * 
    * The block_number argument identifies a specific block
    * by the array index of the array of pointers to hierarchies given
    * the constructor of this class.  A block is the neighbor of
    * another block if the two blocks abut at a point, a 1D line, or
    * a 2D plane
    *
    * @param block_number The block of which the number of neighbors is sought
    */
   int getNumberOfNeighbors(const int block_number);

   enum RotationIdentifier
   {
      NO_ROTATE = 0,
      IUP_JUP = 0,
      JUP_IDOWN = 1,
      IDOWN_JDOWN = 2,
      JDOWN_IUP = 3,
      IUP_JUP_KUP = 0,
      KUP_IUP_JUP = 1,
      JUP_KUP_IUP = 2,
      IDOWN_KUP_JUP = 3,
      KUP_JUP_IDOWN = 4,
      JUP_IDOWN_KUP = 5,
      KDOWN_JUP_IUP = 6,
      IUP_KDOWN_JUP = 7,
      JUP_IUP_KDOWN = 8,
      KDOWN_IDOWN_JUP = 9,
      IDOWN_JUP_KDOWN = 10,
      JUP_KDOWN_IDOWN = 11,
      JDOWN_IUP_KUP = 12,
      IUP_KUP_JDOWN = 13,
      KUP_JDOWN_IUP = 14,
      JDOWN_KUP_IDOWN = 15,
      IDOWN_JDOWN_KUP = 16,
      KUP_IDOWN_JDOWN = 17,
      JDOWN_KDOWN_IUP = 18,
      KDOWN_IUP_JDOWN = 19,
      IUP_JDOWN_KDOWN = 20,
      JDOWN_IDOWN_KDOWN = 21,
      KDOWN_JDOWN_IDOWN = 22,
      IDOWN_KDOWN_JDOWN = 23
   };

   /*!
    * @brief Structure to represent the neighbor of a given block.
    *
    * @param d_id Integer identifier of the neighboring block, based on the
    *             ordering of the array given to this class's constructor
    * @param  d_translated_domain BoxArray describing coarse-level domain of
    *                             the neighboring block in terms of the index
    *                             space of the given block
    * @param d_rotation The rotation needed to rotate the index space of the
    *                   neighboring block so that it is aligned with the index
    *                   space of the given block
    * @param d_shift The shift needed to move the domain of the
    *                neighboring block from its post-rotation location
    *                to its actual location, from the point of view of the
    *                index space of the given block.
    * @param d_is_singularity Whether or not the two blocks touch at a
    *                         singularity
    */
   struct Neighbor {
      int d_id;
      hier::BoxArray<DIM> d_translated_domain;
      RotationIdentifier d_rotation;
      hier::IntVector<DIM> d_shift;
      bool d_is_singularity;
   };

   /*!
    * @brief Return a list of Neighbor objects describing all of the neighbors
    * of the block indicated by the block_number.
    *
    * @param block_number A list of neighbors of the given block is returned
    */
   tbox::List<Neighbor>& getNeighbors(const int block_number);

   /*!
    * @brief Return a BoxList that describes all of the singularities
    * touched by the block indicated by block_number.
    *
    * For every singularity point the block touches, the BoxList will contain
    * a single-cell box that lies just outside the block domain, touching
    * the block only at the singularity point.  For line singularities,
    * the BoxList will contain boxes of width 1 in all dimensions except 1,
    * lying outside the block's coarse-level domain and touching the domain
    * only along the line of singularity.
    *
    * @param block_number Identifies block for which singularityes will be
    *                     retrieved
    */
   hier::BoxList<DIM>& getSingularityBoxList(const int block_number);

   /*!
    * @brief Get a box array that describes the coarse-level domain of the
    * translated_block in terms of the index space of base_block.
    *
    * @param block_boxes Output consisting of the coarse-leve domain of
    *                    the block identified by translated_block, represented
    *                    in the index space of the block identified by
    *                    base_block
    * @param base_block  Integer identifier of the block whose index space
    *                    will be used for the output boxes
    * @param translated_block Integer identifier of another block whose
    *                         domain will be represented in the index space of
    *                         the base block
    */
   void getTranslatedBlock(hier::BoxArray<DIM>& block_boxes,
                           const int base_block,
                           const int translated_block);


   /*!
    * @brief Modify boxes by rotating and shifting from the index space of
    * the translated_block to the index space of the base_block at the
    * resolution level defined by ratio.
    *
    * @param boxes Input and output.  The boxes will be translated from
    *              the translated_block index space to the base_block
    *              index space.
    * @param ratio The boxes will be refined to the ratio given.  This
    *              should represent the refinement ratio of a level in
    *              the AMR hierarchy.
    * @param base_block Integer identifier of the block whose index space
    *                   will be represented in the boxes output
    * @param translated_block Integer identifier of the block whose index
    *                         space is represented in the boxes input
    */
   void translateBoxArray(hier::BoxArray<DIM>& boxes,
                          const hier::IntVector<DIM>& ratio,
                          const int base_block,
                          const int translated_block);

   /*!
    * @brief static routine to get a reverse rotation identifier
    *
    * A rotation identifier signifies a specific rotation of an index space.
    * For each rotation there is another rotation that rotates in the exact
    * opposite manner.  This routine returns the identifier of the reverse
    * rotation corresponding to the given rotation.
    *
    * @param rotation Rotation for which the reverse rotation is sought
    */ 
   static RotationIdentifier getReverseRotationIdentifier(
      const RotationIdentifier rotation);

   /*!
    * @brief static routine to get a reverse shift
    */
   static void calculateReverseShift(hier::IntVector<DIM>& back_shift,
                                     const hier::IntVector<DIM>& shift,
                                     const RotationIdentifier back_rotation);

   /*!
    * @brief Get finest level number existing in multiblock patch hierarchy
    */
   int getFinestLevelNumber() const;

   /*!
    * Return the number of levels that currently exist in the hierarchy.
    */
   int getNumberOfLevels() const;

   /*!
    * @brief Returns true if the array of patch levels contains a patch level
    * finer than the specified patch level. Otherwise, false is returned.
    *
    * @param ln  A given level number
    */
   bool finerLevelExists(const int ln) const;

   /*!
    * Get a BoxList that contains all of the index space of all other blocks
    * in the multiblock domain.
    *
    * A BoxList will be constructed that contains the full set of the
    * coarse level domains of all blocks except the one identified by
    * block_number.  The domains will all be translated into the index space
    * represented by block_number.
    *
    * @param domain_outside_block Output box list
    * @param block_number         domain_outside_block will represent the
    *                             domains of all blocks except this one.
    * 
    */
   void getDomainOutsideBlock(hier::BoxList<DIM>& domain_outside_block,
                              const int block_number);

   /*!
    * @brief Return true if block represented by block_number touches
    * a reduced-connectivity singularity
    *
    * @param block_number Number of block to be tested
    */
   bool reducedConnectivityExists(const int block_number) const;

   /*!
    * @brief Writes the state of the PatchHierarchy object and the PatchLevels
    * it contains to the database.
    * It should be noted that only those patch data which have been
    * registered for restart with the hier::VariableDatabase<DIM> will be
    * written to the database.  This method implements the pure virtual method
    * in tbox::Serializable class which is used by the tbox::RestartManager
    * for writing the MultiblockPatchHierarchy to a restart file.
    *
    * When assertion checking is active, the database pointer must be
    * non-null.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database);

   /*!
    * !brief Writes the state of the MultiblockPatchHierarchy object and the
    * PatchHierarchies it contains to the database.  Only those patchdata
    * corresponding to the set bits in the hier::ComponentSelector are written
    * to the specified database.
    *
    * When assertion checking is active, the database pointer must be
    * non-null.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> database,
                      const hier::ComponentSelector& patchdata_write_table);

   /*!
    * @brief Read in the entire hierarchy from the restart file.
    * The database from which the restart data is read is determined by the
    * object_name specified in the constructor.
    *
    * Notes:
    *
    *
    *    -
    *        This method handles the memory allocation for each PatchLevel
    *        it reads in.
    *
    *    -
    *        The number of levels read in is the minimum of the max_levels
    *        argument and the number of levels stored in the database.
    *
    *
    *
    * When assertion checking is active, the max_levels argument must be
    * greater than zero.  An unrecoverable exception will result if the
    * database cannot be found in the restart file or the data in the
    * restart file is erroneous.
    *
    * @param max_levels maximum number of levels to read in.
    */
   void getFromRestart(const int max_levels);

private:

   /*!
    * @brief Returns true if block a and block b are neighbors
    */
   bool areNeighbors(const int a, const int b);


   /*!
    * @brief Create a neighbor.
    *
    * @param id         The block number of the neighboring block
    * @param domain     The neighboring block's domain in the current block's
    *                   index space
    * @param rotation   The rotation to represent the neighboring block
    *                   in the current index space.
    * @param shift      The shift needed to shift the neighboring block to
    *                   its correct position in the current index space.
    * @param is_singularity Boolean telling whether the current block and
    *                       the neighboring block abut at a singularity
    */
   Neighbor createNeighbor(int id, hier::BoxArray<DIM>& domain,
                           const RotationIdentifier rotation,
                           const hier::IntVector<DIM>& shift, 
                           const bool is_singularity);

   /*!
    * @brief Register a relationship between two neighboring blocks of a
    * multiblock domain.
    *
    * @param block_a         One block in the relationship
    * @param block_b         The other block
    * @param rotation_b_to_a The rotation that align's block b's index space
    *                        with block a's
    * @param shift_b_to_a    The post-rotation shift to move b into its
    *                        correct location within a's index space
    * @param neighbor_type   The type (codimension) of the neighbor
    *                        relationship
    */
   void registerNeighbors(int block_a,
                          int block_b,
                          RotationIdentifier rotation_b_to_a,
                          hier::IntVector<DIM>& shift_b_to_a,
                          const int neighbor_type);

   /*!
    * @brief Map a string identifier of a rotation operation to an integer id.
    *
    * @param rotation_string
    */
   RotationIdentifier
   getRotationIdentifier(const tbox::Array<std::string>& rotation_string) const;

   void adjustBoundaryBoxesOnPatch(
      const hier::Patch<DIM>& patch,
      const tbox::Pointer<hier::GridGeometry<DIM> > geometry,
      const hier::BoxList<DIM>& pseudo_domain,
      const hier::IntVector<DIM>& gcw,
      const hier::BoxList<DIM>& singularity);

   std::string d_object_name;
   bool   d_registered_for_restart;
   

   /*
    * Number of blocks in the hierarchy.
    */
   int d_number_blocks;

   /*
    * Array of hierarchies that cover the multiblock domain
    */
   tbox::Array< tbox::Pointer<hier::PatchHierarchy<DIM> > > d_hierarchies;

   /*
    * Associated with each hierarchy is a list of Neighbors that
    * it shares a block boundary with.   
    */
   tbox::Array< tbox::List<Neighbor> > d_block_neighbors;

   /*
    * An array of BoxLists defining the singularities of a multiblock
    * domain.  Each BoxList element defines the singularities that a single
    * block touches.
    */
   tbox::Array< hier::BoxList<DIM> > d_singularity;

   tbox::Pointer< hier::MultiblockGridGeometry<DIM> > d_geometry;

   tbox::Array<bool> d_reduced_connect;
};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockPatchHierarchy.C"
#endif

#endif
