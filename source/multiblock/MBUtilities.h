//
// File:	MBUtilities.h
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 476 $
// Modified:	$Date: 2005-07-05 15:42:59 -0700 (Tue, 05 Jul 2005) $
// Description:	utility functions for multiblock
//

#ifndef included_mblk_MBUtilities
#define included_mblk_MBUtilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_mblk_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif

namespace SAMRAI {
    namespace mblk {


/*!
 * @brief Class MBUtilities is a utility class to hold static functions
 * related to multiblock functionality.
 *
 * @see mblk_MultiblockPatchHierarchy
 */
template<int DIM>
class MBUtilities
{
public:
   /*! 
    * Empty constructor and destructor.
    */
   MBUtilities();

   virtual ~MBUtilities<DIM>();

   /*!
    * @brief Copy patch data from src to dst using the shift and rotate
    * arguments.
    *
    * @param dst destination data
    * @param src source data
    * @param shift the shift needed after rotation
    * @param rotate identifier of the rotation between index spaces
    */
   static
   void translateAndCopyData(
      hier::PatchData<DIM>& dst,
      const hier::PatchData<DIM>& src,
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
    * @param rotation_number identifier of the rotation that will be applied
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
