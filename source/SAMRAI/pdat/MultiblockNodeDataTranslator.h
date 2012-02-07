/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_MultiblockNodeDataTranslator
#define included_pdat_MultiblockNodeDataTranslator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MultiblockDataTranslator.h"

namespace SAMRAI {
namespace pdat {

/*!
 * Class MultiblockNodeDataTranslator<DIM>
 */

template<class TYPE>
class MultiblockNodeDataTranslator:
   public hier::MultiblockDataTranslator
{
public:
   /*!
    * @brief Constructor
    */
   MultiblockNodeDataTranslator<TYPE>();

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockNodeDataTranslator<TYPE>();

   virtual void
   translateAndCopyData(
      hier::Patch& dst_patch,
      const int dst_id,
      const hier::Patch& src_patch,
      const int src_id,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate);

   virtual void translateAndFillData(
      hier::Patch& dst_patch,
      const int dst_id,
      const hier::Patch& src_patch,
      const int src_id,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate)
   {
      NULL_USE(dst_patch);
      NULL_USE(dst_id);
      NULL_USE(src_patch);
      NULL_USE(src_id);
      NULL_USE(shift);
      NULL_USE(rotate);
   }

private:
};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SAMRAI/pdat/MultiblockNodeDataTranslator.C"
#endif

#endif
