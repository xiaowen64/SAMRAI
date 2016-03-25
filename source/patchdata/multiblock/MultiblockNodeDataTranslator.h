//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/patchdata/boxgeometry/MultiblockNodeDataTranslator.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1318 $
// Modified:	$LastChangedDate: 2006-12-04 12:10:56 -0800 (Mon, 04 Dec 2006) $
// Description:	hier::Box geometry information for cell centered objects
//

#ifndef included_pdat_MultiblockNodeDataTranslator
#define included_pdat_MultiblockNodeDataTranslator

#include "SAMRAI_config.h"
#include "MultiblockDataTranslator.h"


namespace SAMRAI {
    namespace pdat {

/*!
 * Class MultiblockNodeDataTranslator<DIM>
 */

template<int DIM, class TYPE> class MultiblockNodeDataTranslator :
public hier::MultiblockDataTranslator<DIM>
{
public:
   /*!
    * @brief Constructor
    */
   MultiblockNodeDataTranslator<DIM,TYPE>(); 

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockNodeDataTranslator<DIM,TYPE>();

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

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockNodeDataTranslator.C"
#endif
