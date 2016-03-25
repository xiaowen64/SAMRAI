//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/ErrorCheckIntTypes.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Trivial struct types for enforcing type checking in
//              multiple int function arguments.
//

#ifndef included_hier_ErrorCheckIntTypes
#define included_hier_ErrorCheckIntTypes

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif


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


