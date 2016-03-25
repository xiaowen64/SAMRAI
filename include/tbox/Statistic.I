//
// File:        Statistic.I
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    \f$       \f$
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Class to record statistics during program execution.
//

#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*                                                                       *
* Utility functions to retrieve statistic name and type information.    *
*                                                                       *
*************************************************************************
*/

inline string Statistic::getName() const
{
   return(d_object_name);
}

inline string Statistic::getType() const
{
   return( (d_stat_type == PROC_STAT) ? "PROC_STAT" : "PATCH_STAT" );
}

inline int Statistic::getInstanceId() const
{
   return(d_instance_id);
}

inline int Statistic::getStatSequenceLength() const
{
   return(d_seq_counter);
}

inline void Statistic::reset() {
   d_proc_array.resizeArray(0);
   d_patch_array.resizeArray(0);
}

/*
*************************************************************************
*                                                                       *
* Protected utility functions to retrieve statistic record data.        *
*                                                                       *
*************************************************************************
*/

inline const Array<Statistic::ProcStat>&
Statistic::getProcStatSeqArray() const
{
   return(d_proc_array);
}

inline const Array<Statistic::PatchStat>&
Statistic::getPatchStatSeqArray() const
{
   return(d_patch_array);
}

/*
*************************************************************************
*                                                                       *
* Return false indicating that statistic data size cannot be computed   *
* on each processor a priori.  This can probably be optimized later.    *
*                                                                       *
*************************************************************************
*/

inline bool Statistic::canEstimateDataStreamSize()
{ 
   return(false);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
