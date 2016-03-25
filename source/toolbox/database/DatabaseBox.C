//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/database/DatabaseBox.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	A box structure representing a portion of the AMR index space
//

#include "tbox/DatabaseBox.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/DatabaseBox.I"
#endif

namespace SAMRAI {
   namespace tbox {


DatabaseBox::DatabaseBox(const int dimension, const int *lower, const int *upper)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((dimension >= 0) && (dimension <= DatabaseBox_MAX_DIM));
#endif
   d_data.d_dimension = dimension;
   for (int i = 0; i < d_data.d_dimension; i++) {
      d_data.d_lo[i] = lower[i];
      d_data.d_hi[i] = upper[i];
   }
   for (int j = d_data.d_dimension; j < DatabaseBox_MAX_DIM ; j++) {
      d_data.d_lo[j] = 0;
      d_data.d_hi[j] = 0;
   }
}

bool DatabaseBox::empty() const
{
   bool is_empty = (d_data.d_dimension == 0 ? true : false);
   for (int i = 0; i < d_data.d_dimension; i++) {
      if (d_data.d_hi[i] < d_data.d_lo[i]) is_empty = true;
   }
   return(is_empty);
}

int DatabaseBox::operator==(const DatabaseBox& box) const
{
   bool equals = (d_data.d_dimension == box.d_data.d_dimension);
   for (int i = 0; i < d_data.d_dimension; i++) {
      if (d_data.d_lo[i] != box.d_data.d_lo[i]) equals = false;
      if (d_data.d_hi[i] != box.d_data.d_hi[i]) equals = false;
   }
   return(equals);
}


}
}

