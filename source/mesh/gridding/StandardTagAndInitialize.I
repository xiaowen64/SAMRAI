//
// File:        StandardTagAndInitialize.I
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 7 $
// Modified:    $Date: 2004-11-30 13:18:17 -0800 (Tue, 30 Nov 2004) $
// Description: Routines for performing cell-tagging and initializing 
//              a new level.
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace mesh {

template<int DIM> inline
bool StandardTagAndInitialize<DIM>::usesTimeIntegration() const
{
  return(d_use_richardson_extrapolation);
}

template<int DIM> inline
int StandardTagAndInitialize<DIM>::getErrorCoarsenRatio() const
{
   return(d_error_coarsen_ratio);
}

template<int DIM> inline
void StandardTagAndInitialize<DIM>::turnOnGradientDetector()
{
   d_use_gradient_detector = true;
}

template<int DIM> inline
void StandardTagAndInitialize<DIM>::turnOffGradientDetector() 
{
   d_use_gradient_detector = false;
}

template<int DIM> inline
void StandardTagAndInitialize<DIM>::turnOnRichardsonExtrapolation()
{
   d_use_richardson_extrapolation = true;
}

template<int DIM> inline
void StandardTagAndInitialize<DIM>::turnOffRichardsonExtrapolation()
{
   d_use_richardson_extrapolation = false;
}

template<int DIM> inline
void StandardTagAndInitialize<DIM>::turnOnRefineBoxes()
{
   d_use_refine_boxes = true;
}

template<int DIM> inline
void StandardTagAndInitialize<DIM>::turnOffRefineBoxes()
{
   d_use_refine_boxes = false;
}

template<int DIM> inline
bool StandardTagAndInitialize<DIM>::refineUserBoxInputOnly() const
{
   bool use_only_refine_boxes = false;
   if (d_use_refine_boxes) {
      use_only_refine_boxes = true;
      if (d_use_gradient_detector || d_use_richardson_extrapolation) {
         use_only_refine_boxes = false;
      }
   }
   return(use_only_refine_boxes);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
