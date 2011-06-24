/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution. 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *                                                                       *
 * Utility functions to retrieve statistic name and type information.    *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD const std::string& Statistic::getName() const
{
   return d_object_name;
}

SAMRAI_INLINE_KEYWORD std::string Statistic::getType() const
{
   return (d_stat_type == PROC_STAT) ? "PROC_STAT" : "PATCH_STAT";
}

SAMRAI_INLINE_KEYWORD int Statistic::getInstanceId() const
{
   return d_instance_id;
}

SAMRAI_INLINE_KEYWORD int Statistic::getStatSequenceLength() const
{
   return d_seq_counter;
}

SAMRAI_INLINE_KEYWORD void Statistic::reset() {
   d_proc_array.clear();
   d_patch_array.clear();
}

/*
 *************************************************************************
 *                                                                       *
 * Protected utility functions to retrieve statistic record data.        *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD const Array<Statistic::ProcStat>&
Statistic::getProcStatSeqArray() const
{
   return d_proc_array;
}

SAMRAI_INLINE_KEYWORD const Array<Statistic::PatchStat>&
Statistic::getPatchStatSeqArray() const
{
   return d_patch_array;
}

/*
 *************************************************************************
 *                                                                       *
 * Return false indicating that statistic data size cannot be computed   *
 * on each processor a priori.  This can probably be optimized later.    *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD bool Statistic::canEstimateDataStreamSize()
{
   return false;
}

}
}
