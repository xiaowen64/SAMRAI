//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/multiblock/MBDataUtilities.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated operations copying patch data.
//

#ifndef included_pdat_MBDataUtilities
#define included_pdat_MBDataUtilities

#include "SAMRAI_config.h"
#include "CellData.h"
#include "EdgeData.h"
#include "FaceData.h"
#include "MultiblockPatchHierarchy.h"
#include "NodeData.h"
#include "SideData.h"

namespace SAMRAI {
    namespace pdat {

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
 * @see hier::MultiblockPatchHierarchy
 * @see hier::MBUtilities
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
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

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
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

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
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

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
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);


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
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

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
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

private:

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MBDataUtilities.C"
#endif

#endif
