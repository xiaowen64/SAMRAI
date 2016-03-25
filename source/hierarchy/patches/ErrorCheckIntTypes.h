//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/ErrorCheckIntTypes.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2132 $
// Modified:    $LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Trivial struct types for enforcing type checking in
//              multiple int function arguments.
//

#ifndef included_hier_ErrorCheckIntTypes
#define included_hier_ErrorCheckIntTypes

#include "SAMRAI_config.h"


namespace SAMRAI {
    namespace hier {

   /*!
    * The PatchNumber struct associates a type to an integral
    * patch number to prevent errors in function arguments.
    */
   struct PatchNumber {
      int pn;
      explicit PatchNumber(int patch_num) : pn(patch_num) 
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(patch_num >= 0);
#endif
      }
   };

   /*!
    * The LevelNumber struct associates a type to an integral
    * level number to prevent errors in function arguments.
    */
   struct LevelNumber {
      int ln;
      explicit LevelNumber(int level_num) : ln(level_num) { }
   };

   /*!
    * The PatchDataId struct associates a type to an integral
    * patch data index to prevent errors in function arguments.
    */
   struct PatchDataId {
      int pd;
      explicit PatchDataId(int patch_data_id) : pd(patch_data_id) { }
   };


}
}

#endif


