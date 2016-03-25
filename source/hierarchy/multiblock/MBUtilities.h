//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/multiblock/MBUtilities.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	utility functions for multiblock
//

#ifndef included_hier_MBUtilities
#define included_hier_MBUtilities

#include "SAMRAI_config.h"
#include "MultiblockPatchHierarchy.h"

namespace SAMRAI {
    namespace hier {


/*!
 * @brief Class MBUtilities contains utility functions
 * related to multiblock functionality.
 *
 * @see hier::MultiblockPatchHierarchy
 */
template<int DIM>
class MBUtilities
{
public:

   /*!
    * Constructor
    */
   MBUtilities<DIM>();

   /*!
    * Virtual destructor does nothing
    */
   virtual ~MBUtilities<DIM>();

   /*!
    * @brief Copy patch data from src to dst using the shift and rotate
    * arguments.
    *
    * @param dst_patch   destination data
    * @param dst_id      destination id
    * @param src_patch   source data
    * @param src_id      source id
    * @param shift       the shift needed after rotation
    * @param rotate      identifier of the rotation between index spaces
    */
   static
   void translateAndCopyData(
      hier::Patch<DIM>& dst_patch,
      const int dst_id,
      const hier::Patch<DIM>& src_patch,
      const int src_id,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   /*!
    * @brief Fill patch data from src to dst using the shift and rotate
    * arguments.
    *
    * This is an empty virtual function that allows for a user-defined
    * implementation of the filling of destination patch data.
    *
    * @param dst_patch   destination data
    * @param dst_id      destination id
    * @param src_patch   source data
    * @param src_id      source id
    * @param shift       the shift needed after rotation
    * @param rotate      identifier of the rotation between index spaces
    */
   static
   void translateAndFillData(
      hier::Patch<DIM>& dst_patch,
      const int dst_id,
      const hier::Patch<DIM>& src_patch,
      const int src_id,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   /*!
    * @brief rotate an index from one index space to another
    *
    * The parameter index is an int pointer with points to an array of
    * int data, length DIM.  It signifies an ijk location in an index
    * space.  According to the rotation number, the location will be
    * rotated around the origin, with the new values overwriting the original
    * values in the array pointed to by index.
    *
    * @param index array identifying a point in index space
    * @param rotation    identifier of the rotation that will be applied
    *                        to index
    */
   static
   void rotateIndex(
      int* index,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotation);

private:
   /*!
    * @brief private routine to rotate an index around an axis
    *
    * In 3D, rotation of an index about the origin is decomposed into a
    * series of rotations about an axis.  This function performs one such
    * rotation.
    *
    * @param index array identifying a point in index space
    * @param axis axis around which index will be rotated
    * @param num_rotations number of 90-degree rotations around the axis
    */
   static
   void rotateAboutAxis(int* index,
                        const int axis,
                        const int num_rotations);


};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MBUtilities.C"
#endif

#endif
