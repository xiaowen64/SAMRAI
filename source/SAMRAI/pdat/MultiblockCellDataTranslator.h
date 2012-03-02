/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_MultiblockCellDataTranslator
#define included_pdat_MultiblockCellDataTranslator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MultiblockDataTranslator.h"
#include "SAMRAI/pdat/ArrayData.h"

namespace SAMRAI {
namespace pdat {

/*!
 * Class MultiblockCellDataTranslator<DIM>
 */

template<class TYPE>
class MultiblockCellDataTranslator:
   public hier::MultiblockDataTranslator
{
public:
   /*!
    * @brief Constructor
    */
   MultiblockCellDataTranslator<TYPE>();

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~MultiblockCellDataTranslator<TYPE>();

   virtual void
   translateAndCopyData(
      hier::Patch& dst_patch,
      const int dst_id,
      const hier::Patch& src_patch,
      const int src_id,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate);

   virtual void
   translateAndFillData(
      hier::Patch& dst_patch,
      const int dst_id,
      const hier::Patch& src_patch,
      const int src_id,
      const hier::IntVector& shift,
      const hier::Transformation::RotationIdentifier rotate);

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

#include "SAMRAI/pdat/MultiblockCellDataTranslator.C"

#endif
