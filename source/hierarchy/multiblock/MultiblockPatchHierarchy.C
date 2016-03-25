//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/multiblock/MultiblockPatchHierarchy.C $
// Package:     SAMRAI multiblock package
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1887 $
// Modified:    $LastChangedDate: 2008-01-22 15:58:40 -0800 (Tue, 22 Jan 2008) $
// Description: Base class for geometry management on patches
//

#ifndef included_hier_MultiblockPatchHierarchy_C
#define included_hier_MultiblockPatchHierarchy_C

#include "MultiblockPatchHierarchy.h"

#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define MBLK_PATCH_HIERARCHY_VERSION (2)

#define BLOCK_NAME_BUF_SIZE (16)

namespace SAMRAI {
    namespace hier {

/*
 * ************************************************************************
 *                                                                        *
 * The constructor initializes the arry of hierarchies, and reads from    *
 * input all of the neighbor relationships between the various blocks of  *
 * the domain.                                                            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
MultiblockPatchHierarchy<DIM>::MultiblockPatchHierarchy(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   tbox::Pointer< hier::MultiblockGridGeometry<DIM> >& geometry,
   bool register_for_restart)
{
   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;
   
   d_number_blocks = input_db->getInteger("num_blocks");
   d_hierarchies.resizeArray(d_number_blocks);
   d_block_neighbors.resizeArray(d_number_blocks);
   d_geometry = geometry;

   for (int g = 0; g < d_number_blocks; g++) {

      char hier_name[32];
      sprintf(hier_name, "PatchHierarchy%d", g);

      d_hierarchies[g] = new hier::PatchHierarchy<DIM>(
                                hier_name,
                                geometry->getBlockGeometry(g)); 

   }

   d_singularity.resizeArray(d_number_blocks);
   d_reduced_connect.resizeArray(d_number_blocks);

   char sing_name[32];
   char neighbor_name[32];

   for (int i = 0; i < d_number_blocks; i++) {

      d_reduced_connect[i] = false;

   }

   for (int sn = 0; true; sn++) {

      sprintf(sing_name, "Singularity%d", sn);

      if (!input_db->keyExists(sing_name)) {
         break;
      }

      tbox::Pointer<tbox::Database> sing_db =
         input_db->getDatabase(sing_name);

      tbox::Array<int> blocks = sing_db->getIntegerArray("blocks");

      char block_box_name[32];

      for (int nb = 0; nb < blocks.size(); nb++) {
         sprintf(block_box_name, "sing_box_%d", blocks[nb]);

         hier::Box<DIM> sing_box = sing_db->getDatabaseBox(block_box_name);

         d_singularity[blocks[nb]].unionBoxes(sing_box);
      }
   }
 
   for (int bn = 0; true; bn++) {
      sprintf(neighbor_name, "BlockNeighbors%d", bn);

      if (!input_db->keyExists(neighbor_name)) {
         break;
      }
      tbox::Pointer<tbox::Database> pair_db =
         input_db->getDatabase(neighbor_name);

      int block_a = pair_db->getInteger("block_a");
      int block_b = pair_db->getInteger("block_b");
      RotationIdentifier rotation_b_to_a;

      hier::IntVector<DIM> shift(0);
      if (DIM == 1) {
         rotation_b_to_a = NO_ROTATE;
      } else {
         tbox::Array<std::string> rstr =
            pair_db->getStringArray("rotation_b_to_a");
         rotation_b_to_a = getRotationIdentifier(rstr);

         tbox::Array<int> b_array = 
            pair_db->getIntegerArray("point_in_b_space");
         tbox::Array<int> a_array =
            pair_db->getIntegerArray("point_in_a_space");

         hier::Index<DIM> b_index;
         hier::Index<DIM> a_index;

         for (int p = 0; p < DIM; p++) {
            b_index(p) = b_array[p];
            a_index(p) = a_array[p];
         }

         hier::Box<DIM> b_box(b_index, b_index);
         hier::Box<DIM> a_box(a_index, a_index);

         b_box.rotate(rotation_b_to_a);
         hier::Index<DIM> b_rotated_point(b_box.lower());
         hier::Index<DIM> a_point = (a_box.lower());
   
         shift = a_point - b_rotated_point;
      }

      bool is_singularity =
         pair_db->getBoolWithDefault("is_singularity", false);

      registerNeighbors(block_a, block_b,
                        rotation_b_to_a, shift, is_singularity);

   }

   for (int b = 0; b < d_number_blocks; b++) {
      hier::BoxList<DIM> pseudo_domain;
      getDomainOutsideBlock(pseudo_domain, b);

      pseudo_domain.unionBoxes(d_hierarchies[b]->getGridGeometry()->
                                  getPhysicalDomain());

      for (typename hier::BoxList<DIM>::Iterator
           si(d_singularity[b]); si; si++) {
         hier::BoxList<DIM> test_domain(pseudo_domain);
         test_domain.intersectBoxes(si());
         if (test_domain.size() == 0) {
            d_reduced_connect[b] = true;
            break;
         }
      }  
   }

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * The destructor for the multiblock class implicitly deallocates all of  *
 * the data associated with the multiblock.                               *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
MultiblockPatchHierarchy<DIM>::~MultiblockPatchHierarchy()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Register a relationship between two blocks signified by the integer    *
 * identifiers, with the rotation, shift and neighbor type arguments      *
 * describing the exact nature of the neighbor relation.                  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void MultiblockPatchHierarchy<DIM>::registerNeighbors(
   int a,
   int b,
   RotationIdentifier rotation,
   hier::IntVector<DIM>& shift,
   const int is_singularity)
{

   hier::BoxArray<DIM> a_domain = d_hierarchies[a]->getGridGeometry()->
                                getPhysicalDomain(); 
   hier::BoxArray<DIM> b_domain = d_hierarchies[b]->getGridGeometry()->
                                getPhysicalDomain(); 

   hier::BoxArray<DIM> b_domain_in_a_space(b_domain);
   hier::BoxArray<DIM> a_domain_in_b_space(a_domain);

   RotationIdentifier back_rotation = NO_ROTATE;
   hier::IntVector<DIM> back_shift;

   if (DIM == 2) {

      if (rotation == IUP_JUP) {
         back_rotation = IUP_JUP;
         back_shift = -shift;
      } else if (rotation == JUP_IDOWN) {
         back_rotation = JDOWN_IUP;
         back_shift(0) = shift(1);
         back_shift(1) = -shift(0);
      } else if (rotation == IDOWN_JDOWN) {
         back_rotation = IDOWN_JDOWN;
         back_shift(0) = shift(0);
         back_shift(1) = shift(1);
      }  else if (rotation == JDOWN_IUP) {
         back_rotation = JUP_IDOWN;
         back_shift(0) = -shift(1);
         back_shift(1) = shift(0);
      } else {
          TBOX_ERROR("MultiblockPatchHierarchy<DIM>::registerNeighbors error...\n"
          << "  object name = " << d_object_name
          << " Invalid RotationIdentifier value given" << std::endl);
      }
   } else if (DIM == 3) {
      back_rotation = getReverseRotationIdentifier(rotation);
      calculateReverseShift(back_shift, shift, back_rotation);
   } else {
      TBOX_ERROR("MultiblockPatchHierarchy<DIM>::registerNeighbors error...\n"
          << "  object name = " << d_object_name
          << " Multiblock only works for 2D and 3D" << std::endl);
     
   }


   bool rotation_needed;
   if (rotation != 0) {
      rotation_needed = true;
   } else {
      rotation_needed = false;
   }
 
   if (rotation_needed) {
      b_domain_in_a_space.rotate(rotation); 
      a_domain_in_b_space.rotate(back_rotation);
   }
   b_domain_in_a_space.shift(shift);
   a_domain_in_b_space.shift(back_shift);

   Neighbor neighbor_of_b = createNeighbor(a, a_domain_in_b_space,
                                           back_rotation, back_shift,
                                           is_singularity);
   Neighbor neighbor_of_a = createNeighbor(b, b_domain_in_a_space,
                                           rotation, shift,
                                           is_singularity);

   d_block_neighbors[a].addItem(neighbor_of_a);
   d_block_neighbors[b].addItem(neighbor_of_b);

}

template<int DIM> void
MultiblockPatchHierarchy<DIM>::calculateReverseShift(
                                        hier::IntVector<DIM>& back_shift,
                                        const hier::IntVector<DIM>& shift,
                                        const RotationIdentifier back_rotation)
{
   if (DIM == 3) {
      if (back_rotation == IUP_JUP_KUP) {
         back_shift = -shift;
      } else if (back_rotation == KUP_IUP_JUP) {
         back_shift(0) = -shift(2);
         back_shift(1) = -shift(0);
         back_shift(2) = -shift(1);
      } else if (back_rotation == JUP_KUP_IUP) {
         back_shift(0) = -shift(1);
         back_shift(1) = -shift(2);
         back_shift(2) = -shift(0);
      } else if (back_rotation == IDOWN_KUP_JUP) {
         back_shift(0) = shift(0);
         back_shift(1) = -shift(2);
         back_shift(2) = -shift(1);
      } else if (back_rotation == KUP_JUP_IDOWN) {
         back_shift(0) = -shift(2);
         back_shift(1) = -shift(1);
         back_shift(2) = shift(0);
      } else if (back_rotation == JUP_IDOWN_KUP) {
         back_shift(0) = -shift(1);
         back_shift(1) = shift(0);
         back_shift(2) = -shift(2);
      } else if (back_rotation == KDOWN_JUP_IUP) {
         back_shift(0) = shift(2);
         back_shift(1) = -shift(1);
         back_shift(2) = -shift(0);
      } else if (back_rotation == IUP_KDOWN_JUP) {
         back_shift(0) = -shift(0);
         back_shift(1) = shift(2);
         back_shift(2) = -shift(1);
      } else if (back_rotation == JUP_IUP_KDOWN) {
         back_shift(0) = -shift(1);
         back_shift(1) = -shift(0);
         back_shift(2) = shift(2);
      } else if (back_rotation == KDOWN_IDOWN_JUP) {
         back_shift(0) = shift(2);
         back_shift(1) = shift(0);
         back_shift(2) = -shift(1);
      } else if (back_rotation == IDOWN_JUP_KDOWN) {
         back_shift(0) = shift(0);
         back_shift(1) = -shift(1);
         back_shift(2) = shift(2);
      } else if (back_rotation == JUP_KDOWN_IDOWN) {
         back_shift(0) = -shift(1);
         back_shift(1) = shift(2);
         back_shift(2) = shift(0);
      } else if (back_rotation == JDOWN_IUP_KUP) {
         back_shift(0) = shift(1);
         back_shift(1) = -shift(0);
         back_shift(2) = -shift(2);
      } else if (back_rotation == IUP_KUP_JDOWN) {
         back_shift(0) = -shift(0);
         back_shift(1) = -shift(2);
         back_shift(2) = shift(1);
      } else if (back_rotation == KUP_JDOWN_IUP) {
         back_shift(0) = -shift(2);
         back_shift(1) = shift(1);
         back_shift(2) = -shift(0);
      } else if (back_rotation == JDOWN_KUP_IDOWN) {
         back_shift(0) = shift(1);
         back_shift(1) = -shift(2);
         back_shift(2) = shift(0);
      } else if (back_rotation == IDOWN_JDOWN_KUP) {
         back_shift(0) = shift(0);
         back_shift(1) = shift(1);
         back_shift(2) = -shift(2);
      } else if (back_rotation == KUP_IDOWN_JDOWN) {
         back_shift(0) = -shift(2);
         back_shift(1) = shift(0);
         back_shift(2) = shift(1);
      } else if (back_rotation == JDOWN_KDOWN_IUP) {
         back_shift(0) = shift(1);
         back_shift(1) = shift(2);
         back_shift(2) = -shift(0);
      } else if (back_rotation == KDOWN_IUP_JDOWN) {
         back_shift(0) = shift(2);
         back_shift(1) = -shift(0);
         back_shift(2) = shift(1);
      } else if (back_rotation == IUP_JDOWN_KDOWN) {
         back_shift(0) = -shift(0);
         back_shift(1) = shift(1);
         back_shift(2) = shift(2);
      } else if (back_rotation == JDOWN_IDOWN_KDOWN) {
         back_shift(0) = shift(1);
         back_shift(1) = shift(0);
         back_shift(2) = shift(2);
      } else if (back_rotation == KDOWN_JDOWN_IDOWN) {
         back_shift(0) = shift(2);
         back_shift(1) = shift(1);
         back_shift(2) = shift(0);
      } else if (back_rotation == IDOWN_KDOWN_JDOWN) {
         back_shift(0) = shift(0);
         back_shift(1) = shift(2);
         back_shift(2) = shift(1);
      }
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Create a neighbor structure to store the nature of a relationship      *
 * between neighboring blocks.                                            *
 *                                                                        *
 * ************************************************************************
 */
template<int DIM>
typename MultiblockPatchHierarchy<DIM>::Neighbor
MultiblockPatchHierarchy<DIM>::createNeighbor(
   int id,
   hier::BoxArray<DIM>& domain,
   const RotationIdentifier rotate,
   const hier::IntVector<DIM>& shift,
   const bool is_singularity) 
{
   Neighbor neighbor;
   neighbor.d_id = id;
   neighbor.d_translated_domain = domain;
   neighbor.d_rotation = rotate;
   neighbor.d_shift = shift;
   neighbor.d_is_singularity = is_singularity;

   return(neighbor);
}

/*
 * ************************************************************************
 *                                                                        *
 * Get a pointer to the specified heirarchy in the multiblock's array.    *
 *                                                                        *
 * ************************************************************************
 */
template<int DIM> tbox::Pointer< hier::PatchHierarchy<DIM> >&
MultiblockPatchHierarchy<DIM>::getHierarchy(const int hiera_num)
{
   return (d_hierarchies[hiera_num]);
}

/*
 * ************************************************************************
 *                                                                        *
 * Set block to be the domain of translated_block in the index space of   *
 * base_block.                                                            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void
MultiblockPatchHierarchy<DIM>::getTranslatedBlock(hier::BoxArray<DIM>& block,
                                     const int base_block,
                                     const int translated_block)
{
   for (typename tbox::List<Neighbor>::Iterator
        ni(d_block_neighbors[base_block]); ni; ni++) {
      if (ni().d_id == translated_block) {
         block = ni().d_translated_domain;
         break;
      }
   }
}

template<int DIM> typename MultiblockPatchHierarchy<DIM>::RotationIdentifier
MultiblockPatchHierarchy<DIM>::getRotationIdentifier(
   const tbox::Array<std::string>& rotation_string) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(rotation_string.getSize() == DIM);
#endif

   RotationIdentifier id = NO_ROTATE;
   bool is_error = false;

   if (DIM == 2) {
      if (rotation_string[0] == "I_UP") {
         if (rotation_string[1] == "J_UP") {
            id = IUP_JUP; //0;
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "I_DOWN") {
         if (rotation_string[1] == "J_DOWN") {
            id = IDOWN_JDOWN; //2;
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_UP") {
         if (rotation_string[1] == "I_DOWN") {
            id = JUP_IDOWN; //1;
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_DOWN") {
         if (rotation_string[1] == "I_UP") {
            id = JDOWN_IUP; //3;
         } else {
            is_error = true;
         }
      }
   } else if (DIM == 3) {
      if (rotation_string[0] == "I_UP") {
         if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "K_UP") {
               id = IUP_JUP_KUP; //0;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "K_DOWN") {
               id = IUP_JDOWN_KDOWN; //20;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "J_DOWN") {
               id = IUP_KUP_JDOWN; //13;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "J_UP") {
               id = IUP_KDOWN_JUP; //7;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "I_DOWN") {
         if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "K_DOWN") {
               id = IDOWN_JUP_KDOWN; //10;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "K_UP") {
               id = IDOWN_JDOWN_KUP; //16;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "J_UP") {
               id = IDOWN_KUP_JUP; //3;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "J_DOWN") {
               id = IDOWN_KDOWN_JDOWN; //23;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_UP") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "K_DOWN") {
               id = JUP_IUP_KDOWN; //8;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "K_UP") {
               id = JUP_IDOWN_KUP; //5;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "I_UP") {
               id = JUP_KUP_IUP; //2;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "I_DOWN") {
               id = JUP_KDOWN_IDOWN; //11;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "J_DOWN") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "K_UP") {
               id = JDOWN_IUP_KUP; //12;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "K_DOWN") {
               id = JDOWN_IDOWN_KDOWN; //21;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_UP") {
            if (rotation_string[2] == "I_DOWN") {
               id = JDOWN_KUP_IDOWN; //15;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "K_DOWN") {
            if (rotation_string[2] == "I_UP") {
               id = JDOWN_KDOWN_IUP; //18;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "K_UP") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "J_UP") {
               id = KUP_IUP_JUP; //1;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "J_DOWN") {
               id = KUP_IDOWN_JDOWN; //17;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "I_DOWN") {
               id = KUP_JUP_IDOWN; //4;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "I_UP") {
               id = KUP_JDOWN_IUP; //14;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else if (rotation_string[0] == "K_DOWN") {
         if (rotation_string[1] == "I_UP") {
            if (rotation_string[2] == "J_DOWN") {
               id = KDOWN_IUP_JDOWN; //19;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "I_DOWN") {
            if (rotation_string[2] == "J_UP") {
               id = KDOWN_IDOWN_JUP; //9;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_UP") {
            if (rotation_string[2] == "I_UP") {
               id = KDOWN_JUP_IUP; //6;
            } else {
               is_error = true;
            }
         } else if (rotation_string[1] == "J_DOWN") {
            if (rotation_string[2] == "I_DOWN") {
               id = KDOWN_JDOWN_IDOWN; //22;
            } else {
               is_error = true;
            }
         } else {
            is_error = true;
         }
      } else {
         is_error = true;
      }

      if (is_error) {
         TBOX_ERROR("Rotation_input " << rotation_string[0] << " " <<
                    rotation_string[1] << " " << rotation_string[2] <<
                    " is invalid.\n");
      }
   }

   return (id);
}

template<int DIM>
typename MultiblockPatchHierarchy<DIM>::RotationIdentifier
MultiblockPatchHierarchy<DIM>::getReverseRotationIdentifier(
   const RotationIdentifier rotation)
{
   RotationIdentifier reverse_id;

   switch (rotation) {

      case IUP_JUP_KUP:
         reverse_id = IUP_JUP_KUP;
         break;

      case KUP_IUP_JUP:
         reverse_id = JUP_KUP_IUP;
         break;

      case JUP_KUP_IUP:
         reverse_id = KUP_IUP_JUP;
         break;

      case IDOWN_KUP_JUP:
         reverse_id = IDOWN_KUP_JUP;
         break;

      case KUP_JUP_IDOWN:
         reverse_id = KDOWN_JUP_IUP;
         break;

      case JUP_IDOWN_KUP:
         reverse_id = JDOWN_IUP_KUP;
         break;

      case KDOWN_JUP_IUP:
         reverse_id = KUP_JUP_IDOWN;
         break;

      case IUP_KDOWN_JUP:
         reverse_id = IUP_KUP_JDOWN;
         break;

      case JUP_IUP_KDOWN:
         reverse_id = JUP_IUP_KDOWN;
         break;

      case KDOWN_IDOWN_JUP:
         reverse_id = JDOWN_KUP_IDOWN;
         break;

      case IDOWN_JUP_KDOWN:
         reverse_id = IDOWN_JUP_KDOWN;
         break;

      case JUP_KDOWN_IDOWN:
         reverse_id = KDOWN_IUP_JDOWN;
         break;

      case JDOWN_IUP_KUP:
         reverse_id = JUP_IDOWN_KUP;
         break;

      case IUP_KUP_JDOWN:
         reverse_id = IUP_KDOWN_JUP;
         break;

      case KUP_JDOWN_IUP:
         reverse_id = KUP_JDOWN_IUP;
         break;

      case JDOWN_KUP_IDOWN:
         reverse_id = KDOWN_IDOWN_JUP;
         break;

      case IDOWN_JDOWN_KUP:
         reverse_id = IDOWN_JDOWN_KUP;
         break;

      case KUP_IDOWN_JDOWN:
         reverse_id = JDOWN_KDOWN_IUP;
         break;

      case JDOWN_KDOWN_IUP:
         reverse_id = KUP_IDOWN_JDOWN;
         break;

      case KDOWN_IUP_JDOWN:
         reverse_id = JUP_KDOWN_IDOWN;
         break;

      case IUP_JDOWN_KDOWN:
         reverse_id = IUP_JDOWN_KDOWN;
         break;

      case JDOWN_IDOWN_KDOWN:
         reverse_id = JDOWN_IDOWN_KDOWN;
         break;

      case KDOWN_JDOWN_IDOWN:
         reverse_id = KDOWN_JDOWN_IDOWN;
         break;

      case IDOWN_KDOWN_JDOWN:
         reverse_id = IDOWN_KDOWN_JDOWN;
         break;

      default:
          TBOX_ERROR("MultiblockPatchHierarchy<DIM>::getReverseRotationIdentifier error...\n"
          << " Invalid RotationIdentifier value given" << std::endl);

         reverse_id = IUP_JUP_KUP;
         break;
   }

   return(reverse_id);
}

/*
 * ************************************************************************
 *                                                                        *
 * Rotate and shift the boxes in the given array according to the         *
 * rotation and shift that is used to translated the index space of       *
 * translated_block into the index space of base_block.                   *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void
MultiblockPatchHierarchy<DIM>::translateBoxArray(hier::BoxArray<DIM>& boxes,
                                    const hier::IntVector<DIM>& ratio,
                                    const int base_block,
                                    const int translated_block)
{
   for (typename tbox::List<Neighbor>::Iterator
        ni(d_block_neighbors[base_block]); ni; ni++) {
      if (ni().d_id == translated_block) {
         hier::IntVector<DIM> refined_shift = (ni().d_shift) * (ratio);
         boxes.rotate(ni().d_rotation);
         boxes.shift(refined_shift); 
      }
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Return a multiblock level that represents all of the patch levels in   *
 * the multiblock domain that have the given level number.                *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> tbox::Pointer< MultiblockPatchLevel<DIM> >
MultiblockPatchHierarchy<DIM>::getMultiblockPatchLevel(
   const int level_num) const
{
   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > >
      level_array(d_number_blocks);

   for (int b = 0; b < d_number_blocks; b++) {
      if (d_hierarchies[b]->getNumberLevels() > level_num) {
         level_array[b] = d_hierarchies[b]->getPatchLevel(level_num);
      } else {
         level_array[b].setNull();
      }
   }

   tbox::Pointer< MultiblockPatchLevel<DIM> > mb_level =
      new MultiblockPatchLevel<DIM>(level_array);

   return (mb_level);
}

/*
 * ************************************************************************
 *                                                                        *
 * Query to find if the two blocks represented by the integer identifiers *
 * are neighbors.                                                         *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
bool MultiblockPatchHierarchy<DIM>::areNeighbors(const int a, const int b)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(a != b);
#endif

   bool are_neighbors = false;
   for (typename tbox::List<Neighbor>::Iterator
        ni(d_block_neighbors[a]); ni; ni++) {
      if (ni().d_id == b) {
         are_neighbors = true;
         break;
      }
   }

   return (are_neighbors);
}

/*
 * ************************************************************************
 *                                                                        *
 * Returns the number of blocks in the multiblock domain.                 *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
int MultiblockPatchHierarchy<DIM>::getNumberBlocks() const
{
   return (d_number_blocks);
}

template<int DIM>
void MultiblockPatchHierarchy<DIM>::adjustMultiblockPatchLevelBoundaries(
   tbox::Pointer< MultiblockPatchLevel<DIM> > mb_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!mb_level.isNull());
   TBOX_ASSERT(mb_level->getNumberBlocks() == d_number_blocks);
#endif

   for (int mb = 0; mb < d_number_blocks; mb++) {
      tbox::Pointer< hier::PatchLevel<DIM> > patch_level =
          mb_level->getPatchLevelForBlock(mb);

      if (!patch_level.isNull()) {

         hier::IntVector<DIM> gcw =
           patch_level->getPatchDescriptor()->getMaxGhostWidth();
         hier::BoxList<DIM> singularity(getSingularityBoxList(mb));
                                                                                
         singularity.refine(mb_level->getRatio());

         hier::BoxList<DIM> pseudo_domain;

         tbox::List< typename MultiblockPatchHierarchy<DIM>::Neighbor> neighbors =
            getNeighbors(mb);
    
         for (typename tbox::List<Neighbor>::Iterator
              nei(neighbors); nei; nei++) {
            pseudo_domain.unionBoxes(nei().d_translated_domain);
         }

         pseudo_domain.refine(patch_level->getRatio());

         pseudo_domain.unionBoxes(patch_level->getPhysicalDomain());
         pseudo_domain.unionBoxes(singularity);

         pseudo_domain.coalesceBoxes();

         for (typename hier::PatchLevel<DIM>::Iterator sp(patch_level);
              sp; sp++) {
            tbox::Pointer< hier::Patch<DIM> > patch =
               patch_level->getPatch(sp());

            adjustBoundaryBoxesOnPatch(*patch,
                                       patch_level->getGridGeometry(),
                                       pseudo_domain,
                                       gcw,
                                       singularity);

         }
      }
   }
}

template<int DIM>
void MultiblockPatchHierarchy<DIM>::adjustBoundaryBoxesOnPatch(
   const hier::Patch<DIM>& patch,
   const tbox::Pointer< hier::GridGeometry<DIM> > grid_geometry,
   const hier::BoxList<DIM>& pseudo_domain,
   const hier::IntVector<DIM>& gcw,
   const hier::BoxList<DIM>& singularity)
{

   /*
    * Avoid adjusting boundary boxes for the case where we just use
    * a single block, since this is equivalent to not using multiblocks
    * at all.
    */
   if (d_number_blocks > 1) {
      tbox::Array< hier::BoundaryBox<DIM> > boundaries[DIM];
      
      grid_geometry->getBoundaryBoxes(boundaries,
                                      patch.getBox(),
                                      pseudo_domain,
                                      gcw,
                                      hier::IntVector<DIM>(0));

      tbox::Array< hier::BoundaryBox<DIM> > codim_boundaries[DIM];
      tbox::List<int> boundaries_in_sing[DIM];
      for (int codim = 2; codim <= DIM; codim++) {

         codim_boundaries[codim-1] =
            patch.getPatchGeometry()->getCodimensionBoundary(codim);

         int num_boxes = codim_boundaries[codim-1].size();

         for (int n = 0; n < num_boxes; n++) {
            hier::Box<DIM> border_box(codim_boundaries[codim-1][n].getBox());
            hier::BoxList<DIM> sing_test_list(singularity);
            sing_test_list.intersectBoxes(border_box);
            if (sing_test_list.size() != 0) {
               boundaries_in_sing[codim-1].addItem(n);
            }
         }
      }
   
      for (int i = 0; i < DIM; i++) {
         if (boundaries_in_sing[i].size() != 0) {
            int old_size = boundaries[i].size();
            boundaries[i].resizeArray(old_size + boundaries_in_sing[i].size());
            int nb = 0;
            for (tbox::List<int>::Iterator b(boundaries_in_sing[i]); b; b++) {
               boundaries[i][old_size+nb] = codim_boundaries[i][b()];
               boundaries[i][old_size+nb].setIsMultiblockSingularity(true);
               nb++;
            }
         }
         patch.getPatchGeometry()->setCodimensionBoundary(boundaries[i],i+1);
      }

   }

}


/* ************************************************************************
 *                                                                        *
 * Returns the number of neighbors a specific block of the multiblock     *
 * domain has.                                                            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
int MultiblockPatchHierarchy<DIM>::getNumberNeighbors(const int block_number)
{
   return (d_block_neighbors[block_number].getNumberItems());
}

/*
 * ************************************************************************
 *                                                                        *
 * Return a list of Neighbor objects describing all of the neighbors      *
 * of the block indicated by the block_number.                            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
tbox::List< typename MultiblockPatchHierarchy<DIM>::Neighbor>&
MultiblockPatchHierarchy<DIM>::getNeighbors(const int block_number)
{
   return (d_block_neighbors[block_number]);
}


/*
 * ************************************************************************
 *                                                                        *
 * Return a BoxList that describes all of the singularity points touched  *
 * by the block indicated by block_number.                                *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
hier::BoxList<DIM>&
MultiblockPatchHierarchy<DIM>::getSingularityBoxList(const int block_number)
{
   return (d_singularity[block_number]);
}

template<int DIM>
int MultiblockPatchHierarchy<DIM>::getFinestLevelNumber() const
{
   int finest_level_number = -1;
   for (int i = 0; i < d_number_blocks; i++) {
      finest_level_number =
         tbox::MathUtilities<int>::Max(finest_level_number,
                               d_hierarchies[i]->getFinestLevelNumber());
   }

   return (finest_level_number);
}

template<int DIM>
int MultiblockPatchHierarchy<DIM>::getNumberLevels() const
{
   return (getFinestLevelNumber() + 1);
}

template<int DIM>
tbox::Pointer< hier::BasePatchLevel<DIM> >
MultiblockPatchHierarchy<DIM>::getPatchLevel(const int l) const
{
   tbox::Pointer< MultiblockPatchLevel<DIM> > mlevel =
      getMultiblockPatchLevel(l);
   tbox::Pointer< hier::BasePatchLevel<DIM> > level = mlevel;
   return(level);
}

template<int DIM>
bool MultiblockPatchHierarchy<DIM>::finerLevelExists(const int l) const
{
   bool finer_level_exists = false;
   for (int i = 0; i < d_number_blocks; i++) {
      if (d_hierarchies[i]->finerLevelExists(l)) {
         finer_level_exists = true;
         break;
      }
   }

   return (finer_level_exists);
}

template<int DIM> void
MultiblockPatchHierarchy<DIM>::getDomainOutsideBlock(
   hier::BoxList<DIM>& domain_outside_block,
   const int block_number)
{
   for (typename tbox::List<Neighbor>::Iterator
        nei(d_block_neighbors[block_number]); nei; nei++) {
      domain_outside_block.unionBoxes(nei().d_translated_domain);
   }
}

template<int DIM> bool
MultiblockPatchHierarchy<DIM>::reducedConnectivityExists(
   const int block_number) const
{
   return (d_reduced_connect[block_number]);
}

/*
*************************************************************************
*                                                                       *
* Writes out the class version number and the number of levels in the   *
* hierarchy and has each patch_level write itself out.                  *
* The database keys for the patch levels are given by                   *
* "level#" where # is the level number for the patch_level.             *
* The patchdata that are written to the database are determined by      *
* which those bits in the hier::VariableDatabase<DIM>                   *
* d_patchdata_restart_table that are set.                               *
*                                                                       *
* Asserts that the database pointer passed in is not NULL.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MultiblockPatchHierarchy<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("MBLK_PATCH_HIERARCHY_VERSION",
                         MBLK_PATCH_HIERARCHY_VERSION);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( BLOCK_NAME_BUF_SIZE > (5 + 1 + 4 + 1) );
#endif
   char block_name[BLOCK_NAME_BUF_SIZE];

   for (int nb = 0; nb < d_number_blocks; nb++) {
      sprintf(block_name,"block_%04d",nb);

      tbox::Pointer<tbox::Database> block_database =
         database->putDatabase(block_name);

      d_hierarchies[nb]->putToDatabase(block_database);
   }
}

/*
*************************************************************************
*                                                                       *
* Writes out the class version number and the number of levels in the   *
* hierarchy and has each patch_level write itself out.                  *
* The database keys for the patch levels are given by                   *
* "level#" where # is the level number for the patch_level.             *
* The patchdata that are written to the database are determined by      *
* which those bits in the specified hier::ComponentSelector that are     *
* set.                                                                  *
*                                                                       *
* Asserts that the database pointer passed in is not NULL.              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void MultiblockPatchHierarchy<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> database,
   const hier::ComponentSelector& patchdata_write_table)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!database.isNull());
#endif

   database->putInteger("MBLK_PATCH_HIERARCHY_VERSION",
                         MBLK_PATCH_HIERARCHY_VERSION);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( BLOCK_NAME_BUF_SIZE > (5 + 1 + 4 + 1) );
#endif
   char block_name[BLOCK_NAME_BUF_SIZE];

   for (int nb = 0; nb < d_number_blocks; nb++) {
      sprintf(block_name,"block_%04d",nb);

      tbox::Pointer<tbox::Database> block_database =
         database->putDatabase(block_name);

      d_hierarchies[nb]->putToDatabase(block_database, patchdata_write_table);
   }
}

/*
*************************************************************************
*                                                                       *
* Gets the database in the root database that corresponds to the object *
* name.  This method then checks the class version against restart      *
* file version.  If they match, it creates each hierarchy level and     *
* reads in the level data.   The number of levels read from restart is  *
* the minimum of the argument max levels and the number of levels in    *
* the restart file.                                                     *
*                                                                       *
*************************************************************************
*/
template<int DIM> void MultiblockPatchHierarchy<DIM>::getFromRestart(
   const int max_levels)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_levels > 0);
#endif

   tbox::Pointer<tbox::Database> restart_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> database;

   if ( restart_db->isDatabase(d_object_name) ) {
      database = restart_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("MultiblockPatchHierarchy<DIM>::getFromRestart() error...\n"
              << "   Restart database with name "
              << d_object_name << " not found in restart file" << std::endl);
   }

   int ver = database->getInteger("MBLK_PATCH_HIERARCHY_VERSION");
   if (ver != MBLK_PATCH_HIERARCHY_VERSION) {
      TBOX_ERROR("MultiblockPatchHierarchy<DIM>::getFromRestart error...\n"
          << "  object name = " << d_object_name
          << " : Restart file version different than class version" << std::endl);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( BLOCK_NAME_BUF_SIZE > (5 + 1 + 4 + 1) );
#endif

   char block_name[BLOCK_NAME_BUF_SIZE];

   for (int nb = 0; nb < d_number_blocks; nb++) {
      sprintf(block_name,"block_%04d",nb);

      tbox::Pointer<tbox::Database> block_database = 
         database->getDatabase(block_name);
  
      d_hierarchies[nb]->getFromDatabase(
         block_database,
         hier::VariableDatabase<DIM>::getDatabase()->getPatchDataRestartTable(),
         max_levels);
   }
}

}
}


#endif
