//
// File:	MultiblockCoarsenAlgorithm.C
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 713 $
// Modified:	$Date: 2005-11-08 13:48:27 -0800 (Tue, 08 Nov 2005) $
// Description:	Coarsening algorithm for data transfer between AMR levels
//

#ifndef included_mblk_MultiblockCoarsenAlgorithm_C
#define included_mblk_MultiblockCoarsenAlgorithm_C
 
#include "MultiblockCoarsenAlgorithm.h"

#include "VariableDatabase.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

namespace SAMRAI {
    namespace mblk {

/*
*************************************************************************
*                                                                       *
* The constructor creates a new xfer::CoarsenClasses<DIM> object             *
* and caches a boolean indiating whether to copy data to the            *
* destination space on the coarse level before coarsening.              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
MultiblockCoarsenAlgorithm<DIM>::MultiblockCoarsenAlgorithm(
   tbox::Pointer< MultiblockPatchHierarchy<DIM> > multiblock,
   bool fill_coarse_data)
{
   d_multiblock = multiblock;
   d_fill_coarse_data = fill_coarse_data;
   d_coarsen_classes = new xfer::CoarsenClasses<DIM>(d_fill_coarse_data);
}

/*
*************************************************************************
*									*
* The destructor implicitly deallocates the list data.                  *
*									*
*************************************************************************
*/

template<int DIM>
MultiblockCoarsenAlgorithm<DIM>::~MultiblockCoarsenAlgorithm()
{
}

/*
*************************************************************************
*                                                                       *
* Register a coarsening operation with the coarsening algorithm.        *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void MultiblockCoarsenAlgorithm<DIM>::registerCoarsen(
   const int dst,
   const int src,
   const tbox::Pointer< xfer::CoarsenOperator<DIM> > opcoarsen,
   const hier::IntVector<DIM>& gcw_to_coarsen)
{
   typename xfer::CoarsenClasses<DIM>::Data data;

   data.d_dst              = dst;
   data.d_src              = src;
   data.d_gcw_to_coarsen   = gcw_to_coarsen;
   data.d_opcoarsen        = opcoarsen;

   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   tbox::Pointer<hier::Variable<DIM> > var;
   if (!var_db->mapIndexToVariable(dst, var)) {
      TBOX_ERROR("MultiblockCoarsenAlgorithm<DIM>::registerCoarsen error..."
                 << "\nNo variable associated with dst patch data index." << endl);
   }

   data.d_fine_bdry_reps_var = var->fineBoundaryRepresentsVariable();

   d_coarsen_classes->insertEquivalenceClassItem(data);
}

/*
*************************************************************************
*                                                                       *
* Create a communication schedule that will coarsen data from fine      *
* patch level to the coarse patch level.                                *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer< MultiblockCoarsenSchedule<DIM> >
MultiblockCoarsenAlgorithm<DIM>::createSchedule(
   tbox::Pointer< MultiblockPatchLevel<DIM> > crse_level,
   tbox::Pointer< MultiblockPatchLevel<DIM> > fine_level,
   MultiblockCoarsenPatchStrategy<DIM>* patch_strategy,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy) const
{
   return(new MultiblockCoarsenSchedule<DIM>(crse_level,
                                             fine_level,
                                             d_coarsen_classes,
                                             d_multiblock,
                                             patch_strategy,
                                             refine_strategy,
                                             d_fill_coarse_data));
}


}
}
#endif
