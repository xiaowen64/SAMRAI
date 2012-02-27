/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility class to read and write boxes to an HDF database.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * Return number of levels in the Database.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
BoxIOUtility::getNumberOfLevels()
{
   return d_level_boxes.getSize();
}

/*
 *************************************************************************
 *
 * Return number of entries for a particular level.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int
BoxIOUtility::getNumberOfEntries(
   const int level_number)
{
   if (level_number + 1 > d_level_boxes.getSize()) {
      TBOX_ERROR("BoxIOUtility::getNumberOfEntries() error: "
         << "invalid level number. "
         << "\n The level boxes database holds data only up "
         << "\n to level " << d_level_boxes.getSize() - 1
         << "\n You requested data from level " << level_number
         << std::endl);
   }

   return d_level_boxes[level_number].getSize();
}

}
}
