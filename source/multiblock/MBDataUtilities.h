//
// File:	MBDataUtilities.h
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 649 $
// Modified:	$Date: 2005-10-05 11:04:02 -0700 (Wed, 05 Oct 2005) $
// Description:	Templated operations copying patch data.
//

#ifndef included_mblk_MBDataUtilities
#define included_mblk_MBDataUtilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_pdat_CellData
#include "CellData.h"
#endif
#ifndef included_pdat_EdgeData
#include "EdgeData.h"
#endif
#ifndef included_pdat_FaceData
#include "FaceData.h"
#endif
#ifndef included_mblk_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif
#ifndef included_pdat_NodeData
#include "NodeData.h"
#endif
#ifndef included_pdat_SideData
#include "SideData.h"
#endif

namespace SAMRAI {
    namespace mblk {

/*!
 * @brief Class MBDataUtilities<DIM,TYPE> is a templated utilitiy class that
 * contains a set of static member functions that can be used to copy
 * patch data between index spaces that are not necessarily aligned
 * on the same axes.
 *
 * This class currently contains functions to copy cell, edge, node, face,
 * and side-centered data, as well as array data.
 *
 * @see hier::PatchData
 * @see mblk::MultiblockPatchHierarchy
 * @see mblk::MBUtilities
 */

template<int DIM, class TYPE>
class MBDataUtilities
{
public:

   /*! 
    * Empty constructor and destructor.
    */
   MBDataUtilities();

   virtual ~MBDataUtilities<DIM,TYPE>();

   /*!
    * @brief Translate and copy cell data from src to dst according the shift
    * and rotation.
    *
    * @param dst destination data
    * @param src source data
    * @param shift shift needed after rotation
    * @param rotate identifier of the rotation between two index spaces 
    */
   static
   void translateAndCopyCellData(
      pdat::CellData<DIM,TYPE>& dst,
      const pdat::CellData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   /*!
    * @brief Translate and copy node data from src to dst according the shift
    * and rotation.
    *
    * @param dst destination data
    * @param src source data
    * @param shift shift needed after rotation
    * @param rotate identifier of the rotation between two index spaces
    */
   static
   void translateAndCopyNodeData(
      pdat::NodeData<DIM,TYPE>& dst,
      const pdat::NodeData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   /*!
    * @brief Translate and copy face data from src to dst according the shift
    * and rotation.
    *
    * @param dst destination data
    * @param src source data
    * @param shift shift needed after rotation
    * @param rotate identifier of the rotation between two index spaces
    */
   static
   void translateAndCopyFaceData(
      pdat::FaceData<DIM,TYPE>& dst,
      const pdat::FaceData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   /*!
    * @brief Translate and copy side data from src to dst according the shift
    * and rotation.
    *
    * @param dst destination data
    * @param src source data
    * @param shift shift needed after rotation
    * @param rotate identifier of the rotation between two index spaces
    */
   static
   void translateAndCopySideData(
      pdat::SideData<DIM,TYPE>& dst,
      const pdat::SideData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);


   /*!
    * @brief Translate and copy edge data from src to dst according the shift
    * and rotation.
    *
    * @param dst destination data
    * @param src source data
    * @param shift shift needed after rotation
    * @param rotate identifier of the rotation between two index spaces
    */
   static
   void translateAndCopyEdgeData(
      pdat::EdgeData<DIM,TYPE>& dst,
      const pdat::EdgeData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   /*!
    * @brief Translate and copy array data from src to dst according the shift
    * and rotation.
    *
    * @param dst destination data
    * @param src source data
    * @param shift shift needed after rotation
    * @param rotate identifier of the rotation between two index spaces
    */
   static
   void translateAndCopyArrayData(
      pdat::ArrayData<DIM,TYPE>& dst,
      const pdat::ArrayData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

private:

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MBDataUtilities.C"
#endif

#endif
