//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/multiblock/MBDataUtilities.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated operations copying patch data.
//

#ifndef included_pdat_MBDataUtilities
#define included_pdat_MBDataUtilities

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
#ifndef included_hier_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif
#ifndef included_pdat_NodeData
#include "NodeData.h"
#endif
#ifndef included_pdat_SideData
#include "SideData.h"
#endif

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
