/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A database structure that stores HDF5 format data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Public method to return the group_id so VisIt can access an
 * object's HDF database.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
hid_t
HDFDatabase::getGroupId()
{
   return d_group_id;
}

}
}
