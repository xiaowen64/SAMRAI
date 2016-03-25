//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/timers/Statistic.I $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    \f$       \f$
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

inline std::string Statistic::getName() const
{
   return(d_object_name);
}

inline std::string Statistic::getType() const
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
