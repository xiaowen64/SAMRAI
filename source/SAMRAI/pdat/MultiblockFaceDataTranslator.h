/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier 
 *
 ************************************************************************/

#ifndef included_pdat_MultiblockFaceDataTranslator
#define included_pdat_MultiblockFaceDataTranslator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/MultiblockDataTranslator.h"

namespace SAMRAI {
namespace pdat {

/*!
 * Class MultiblockFaceDataTranslator<DIM>
 */

template<class TYPE>
class MultiblockFaceDataTranslator:
   public hier::MultiblockDataTranslator
{
public:
   /*!
    * @brief Constructor
    */
   MultiblockFaceDataTranslator<TYPE>();

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockFaceDataTranslator<TYPE>();

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
      (void)dst_patch;
      (void)dst_id;
      (void)src_patch;
      (void)src_id;
      (void)shift;
      (void)rotate;
   }

private:
   void
   translateAndCopyArrayData(
      ArrayData<TYPE>& dst,
      const ArrayData<TYPE>& src,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate);

};

}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SAMRAI/pdat/MultiblockFaceDataTranslator.C"
#endif

#endif
