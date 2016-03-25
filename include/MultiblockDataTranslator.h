//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/multiblock/MultiblockGridGeometry.h $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description:	data translator for Multiblock.
//
 
#ifndef included_hier_MultiblockDataTranslator
#define included_hier_MultiblockDataTranslator

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_MultiblockPatchHierarchy
#include "MultiblockPatchHierarchy.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace hier {

/*!
 * @brief Class MultiblockDataTranslator<DIM>
 */

template<int DIM>
class MultiblockDataTranslator
{
public:
   /*!
    * @brief Default constructor
    */
   MultiblockDataTranslator();
 
   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockDataTranslator<DIM>();

   virtual void translateAndCopyData(
      Patch<DIM>& dst_patch,
      const int dst_id,
      const Patch<DIM>& src_patch,
      const int src_id,
      const IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
      = 0;

   virtual void translateAndFillData(
      Patch<DIM>& dst_patch,
      const int dst_id,
      const Patch<DIM>& src_patch,
      const int src_id,
      const IntVector<DIM>& shift,
      const typename MultiblockPatchHierarchy<DIM>::RotationIdentifier rotate)
       = 0;

private:

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockDataTranslator.C"
#endif

#endif
