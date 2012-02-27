/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   data translator for Multiblock.
 *
 ************************************************************************/

#ifndef included_hier_MultiblockDataTranslator
#define included_hier_MultiblockDataTranslator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class MultiblockDataTranslator
 */

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
   virtual ~MultiblockDataTranslator();

   virtual void
   translateAndCopyData(
      Patch& dst_patch,
      const int dst_id,
      const Patch& src_patch,
      const int src_id,
      const IntVector& shift,
      const Transformation::RotationIdentifier rotate) = 0;

   virtual void
   translateAndFillData(
      Patch& dst_patch,
      const int dst_id,
      const Patch& src_patch,
      const int src_id,
      const IntVector& shift,
      const Transformation::RotationIdentifier rotate) = 0;

private:
};

}
}

#endif
