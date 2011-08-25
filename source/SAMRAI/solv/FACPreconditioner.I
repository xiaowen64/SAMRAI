/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   FAC algorithm for solving linear equations on a hierarchy
 *
 ************************************************************************/

#include <algorithm>
#include <cctype>

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *                                                                       *
 * Utility routines to get information about solver.                     *
 *                                                                       *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
int FACPreconditioner::getNumberOfIterations() const
{
   return d_number_iterations;
}

SAMRAI_INLINE_KEYWORD
double FACPreconditioner::getResidualNorm() const
{
   return d_residual_norm;
}

SAMRAI_INLINE_KEYWORD
void FACPreconditioner::getConvergenceFactors(
   double& avg_factor,
   double& final_factor) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_number_iterations <= 0) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                               << "a solve is invalid.\n");
   }
#endif
   avg_factor = d_avg_convergence_factor;
   final_factor = d_convergence_factor[d_number_iterations - 1];
}

SAMRAI_INLINE_KEYWORD
double FACPreconditioner::getNetConvergenceFactor() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_number_iterations <= 0) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                               << "a solve is invalid.\n");
   }
#endif
   return d_net_convergence_factor;
}

SAMRAI_INLINE_KEYWORD
double FACPreconditioner::getAvgConvergenceFactor() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_number_iterations <= 0) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                               << "a solve is invalid.\n");
   }
#endif
   return d_avg_convergence_factor;
}

SAMRAI_INLINE_KEYWORD
double FACPreconditioner::getFinalConvergenceFactor() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_number_iterations <= 0) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                               << "a solve is invalid.\n");
   }
#endif
   return d_convergence_factor[d_number_iterations - 1];
}

SAMRAI_INLINE_KEYWORD
void FACPreconditioner::setPresmoothingSweeps(
   int num_pre_sweeps) {
   d_presmoothing_sweeps = num_pre_sweeps;
}

SAMRAI_INLINE_KEYWORD
void FACPreconditioner::setPostsmoothingSweeps(
   int num_post_sweeps) {
   d_postsmoothing_sweeps = num_post_sweeps;
}

SAMRAI_INLINE_KEYWORD
void FACPreconditioner::setMaxCycles(
   int max_cycles) {
   d_max_iterations = max_cycles;
}

SAMRAI_INLINE_KEYWORD
void FACPreconditioner::setResidualTolerance(
   double residual_tol,
   double relative_residual_tol) {
   d_residual_tolerance = residual_tol;
   d_relative_residual_tolerance = relative_residual_tol;
}

SAMRAI_INLINE_KEYWORD
void FACPreconditioner::setAlgorithmChoice(
   const std::string& choice)
{
   /* This ptr_function helps resolve to the correct tolower method */
   int (* ptr_function)(
      int) = std::tolower;
   std::string lower = choice;
   std::transform(lower.begin(),
      lower.end(),
      lower.begin(),
      ptr_function);
#ifdef DEBUG_CHECK_ASSERTIONS
   if (lower != "default"               /* Recursive from BTNG */
       && lower != "mccormick-s4.3"     /* McCormick's section 4.3 */
       && lower != "pernice"            /* Translation of Pernice's */
       ) {
      TBOX_ERROR(
         d_object_name << ": algorithm should be set to one of\n"
                       << "'default' (recommended), 'mccormick-s4.3' or 'pernice'\n");
   }
#endif
   d_algorithm_choice = lower;
}

/*
 *************************************************************************
 *                                                                       *
 * Functions to set logging parameters.                                  *
 *                                                                       *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void FACPreconditioner::enableLogging(
   bool enabled) {
   d_do_log = enabled;
}

SAMRAI_INLINE_KEYWORD
const std::string& FACPreconditioner::getObjectName() const {
   return d_object_name;
}

}
}
