//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/patchdata/boxgeometry/MultiblockCellDataTranslator.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1318 $
// Modified:	$LastChangedDate: 2006-12-04 12:10:56 -0800 (Mon, 04 Dec 2006) $
// Description:	hier::Box geometry information for cell centered objects
//

#ifndef included_pdat_MultiblockCellDataTranslator
#define included_pdat_MultiblockCellDataTranslator

#include "SAMRAI_config.h"
#include "MultiblockDataTranslator.h"
#include "ArrayData.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * Class MultiblockCellDataTranslator<DIM>
 */

template<int DIM, class TYPE> class MultiblockCellDataTranslator :
public hier::MultiblockDataTranslator<DIM>
{
public:
   /*!
    * @brief Constructor
    */
   MultiblockCellDataTranslator<DIM,TYPE>(); 

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockCellDataTranslator<DIM,TYPE>();

   virtual void translateAndCopyData(
      hier::Patch<DIM>& dst_patch,
      const int dst_id,
      const hier::Patch<DIM>& src_patch,
      const int src_id,
      const hier::IntVector<DIM>& shift,
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);

   virtual void translateAndFillData(
      hier::Patch<DIM>& dst_patch,
      const int dst_id,
      const hier::Patch<DIM>& src_patch,
      const int src_id,
      const hier::IntVector<DIM>& shift,
      const typename hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
{
   (void) dst_patch;
   (void) dst_id;
   (void) src_patch;
   (void) src_id;
   (void) shift;
   (void) rotate;
}

private:

   void translateAndCopyArrayData(
      ArrayData<DIM,TYPE>& dst,
      const ArrayData<DIM,TYPE>& src,
      const hier::IntVector<DIM>& shift,
      const typename
         hier::MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate);


};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockCellDataTranslator.C"
#endif
