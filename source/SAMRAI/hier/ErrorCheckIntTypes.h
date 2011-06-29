/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Trivial struct types for enforcing type checking in
 *                multiple int function arguments. 
 *
 ************************************************************************/

#ifndef included_hier_ErrorCheckIntTypes
#define included_hier_ErrorCheckIntTypes

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace hier {

/*!
 * The PatchNumber struct associates a type to an integral
 * patch number to prevent errors in function arguments.
 */
struct PatchNumber {
   int pn;
   explicit PatchNumber(
      int patch_num):pn(patch_num)
   {
      TBOX_ASSERT(patch_num >= 0);
   }
};

/*!
 * The LevelNumber struct associates a type to an integral
 * level number to prevent errors in function arguments.
 */
struct LevelNumber {
   int ln;
   explicit LevelNumber(
      int level_num):ln(level_num) {
   }
};

/*!
 * The PatchDataId struct associates a type to an integral
 * patch data index to prevent errors in function arguments.
 */
struct PatchDataId {
   int pd;
   explicit PatchDataId(
      int patch_data_id):pd(patch_data_id) {
   }
};

}
}

#endif
