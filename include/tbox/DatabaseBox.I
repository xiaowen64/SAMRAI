//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/database/DatabaseBox.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	A box structure representing a portion of the AMR index space
//

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
DatabaseBox::DatabaseBox()
{
   d_data.d_dimension = 0;
   d_data.d_lo[0] = d_data.d_hi[0] = 0;
   d_data.d_lo[1] = d_data.d_hi[1] = 0;
   d_data.d_lo[2] = d_data.d_hi[2] = 0;
}

inline
DatabaseBox::DatabaseBox(const DatabaseBox& box)
{
   d_data.d_dimension = box.d_data.d_dimension;
   d_data.d_lo[0]     = box.d_data.d_lo[0];
   d_data.d_lo[1]     = box.d_data.d_lo[1];
   d_data.d_lo[2]     = box.d_data.d_lo[2];
   d_data.d_hi[0]     = box.d_data.d_hi[0];
   d_data.d_hi[1]     = box.d_data.d_hi[1];
   d_data.d_hi[2]     = box.d_data.d_hi[2];
}

inline
DatabaseBox& DatabaseBox::operator=(const DatabaseBox& box)
{
   d_data.d_dimension = box.d_data.d_dimension;
   d_data.d_lo[0]     = box.d_data.d_lo[0];
   d_data.d_lo[1]     = box.d_data.d_lo[1];
   d_data.d_lo[2]     = box.d_data.d_lo[2];
   d_data.d_hi[0]     = box.d_data.d_hi[0];
   d_data.d_hi[1]     = box.d_data.d_hi[1];
   d_data.d_hi[2]     = box.d_data.d_hi[2];
   return(*this);
}

inline
DatabaseBox::~DatabaseBox()
{
}

inline
int DatabaseBox::getDimension() const
{
   return(d_data.d_dimension);
}

inline
void DatabaseBox::setDimension(const int dimension)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((dimension >= 0) && (dimension <= DatabaseBox_MAX_DIM));
#endif
   d_data.d_dimension = dimension;
}

inline
int& DatabaseBox::lower(const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
#endif
   return(d_data.d_lo[i]);
}

inline
int& DatabaseBox::upper(const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
#endif
   return(d_data.d_hi[i]);
}

inline
int DatabaseBox::lower(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
#endif
   return(d_data.d_lo[i]);
}

inline
int DatabaseBox::upper(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
#endif
   return(d_data.d_hi[i]);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
