/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for error reporting, file manipulation, etc.
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 * Routine to rename a file.
 */
SAMRAI_INLINE_KEYWORD
void
Utilities::renameFile(
   const std::string& old_filename,
   const std::string& new_filename)
{
   TBOX_ASSERT(!old_filename.empty());
   TBOX_ASSERT(!new_filename.empty());

   rename(old_filename.c_str(), new_filename.c_str());
}

SAMRAI_INLINE_KEYWORD
std::string
Utilities::nodeToString(
   int num)
{
   return intToString(num, s_node_width);
}

SAMRAI_INLINE_KEYWORD
std::string
Utilities::processorToString(
   int num)
{
   return intToString(num, s_processor_width);
}

SAMRAI_INLINE_KEYWORD
std::string
Utilities::patchToString(
   int num) {
   return intToString(num, s_patch_width);
}

SAMRAI_INLINE_KEYWORD
std::string
Utilities::levelToString(
   int num) {
   return intToString(num, s_level_width);
}

SAMRAI_INLINE_KEYWORD
std::string
Utilities::blockToString(
   int num) {
   return intToString(num, s_block_width);
}

}
}
