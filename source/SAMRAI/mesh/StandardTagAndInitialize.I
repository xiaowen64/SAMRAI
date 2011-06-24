/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Routines for performing cell-tagging and initializing
 *                a new level. 
 *
 ************************************************************************/

namespace SAMRAI {
namespace mesh {

SAMRAI_INLINE_KEYWORD
bool StandardTagAndInitialize::usesTimeIntegration() const
{
   return d_use_richardson_extrapolation;
}

SAMRAI_INLINE_KEYWORD
int StandardTagAndInitialize::getErrorCoarsenRatio() const
{
   return d_error_coarsen_ratio;
}

SAMRAI_INLINE_KEYWORD
void StandardTagAndInitialize::turnOnGradientDetector()
{
   d_use_gradient_detector = true;
}

SAMRAI_INLINE_KEYWORD
void StandardTagAndInitialize::turnOffGradientDetector()
{
   d_use_gradient_detector = false;
}

SAMRAI_INLINE_KEYWORD
void StandardTagAndInitialize::turnOnRichardsonExtrapolation()
{
   d_use_richardson_extrapolation = true;
}

SAMRAI_INLINE_KEYWORD
void StandardTagAndInitialize::turnOffRichardsonExtrapolation()
{
   d_use_richardson_extrapolation = false;
}

SAMRAI_INLINE_KEYWORD
void StandardTagAndInitialize::turnOnRefineBoxes()
{
   d_use_refine_boxes = true;
}

SAMRAI_INLINE_KEYWORD
void StandardTagAndInitialize::turnOffRefineBoxes()
{
   d_use_refine_boxes = false;
}

SAMRAI_INLINE_KEYWORD
bool StandardTagAndInitialize::refineUserBoxInputOnly() const
{
   bool use_only_refine_boxes = false;
   if (d_use_refine_boxes) {
      use_only_refine_boxes = true;
      if (d_use_gradient_detector || d_use_richardson_extrapolation) {
         use_only_refine_boxes = false;
      }
   }
   return use_only_refine_boxes;
}

SAMRAI_INLINE_KEYWORD
const StandardTagAndInitializeConnectorWidthRequestor
&StandardTagAndInitialize::getConnectorWidthRequestor() const
{
   return d_staicwri;
}

}
}
