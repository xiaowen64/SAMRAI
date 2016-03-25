//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/FAC/FACPreconditioner.I $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	FAC algorithm for solving linear equations on a hierarchy
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace solv {


/*
*************************************************************************
*                                                                       *
* Utility routines to get information about solver.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> inline
int FACPreconditioner<DIM>::getNumberIterations() const
{ 
   return(d_number_iterations);
}

template<int DIM> inline
double FACPreconditioner<DIM>::getResidualNorm() const
{ 
   return(d_residual_norm);
}

template<int DIM> inline
void FACPreconditioner<DIM>::getConvergenceFactors(
   double &avg_factor, 
   double &final_factor) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if( d_number_iterations <= 0 ) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
		 << "a solve is invalid.\n");
   }
#endif
   avg_factor = d_avg_convergence_factor;
   final_factor = d_convergence_factor[d_number_iterations-1];
   return;
}


template<int DIM> inline
double FACPreconditioner<DIM>::getNetConvergenceFactor() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if( d_number_iterations <= 0 ) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
		 << "a solve is invalid.\n");
   }
#endif
   return d_net_convergence_factor;
}


template<int DIM> inline
double FACPreconditioner<DIM>::getAvgConvergenceFactor() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if( d_number_iterations <= 0 ) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
		 << "a solve is invalid.\n");
   }
#endif
   return d_avg_convergence_factor;
}


template<int DIM> inline
double FACPreconditioner<DIM>::getFinalConvergenceFactor() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if( d_number_iterations <= 0 ) {
      TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
		 << "a solve is invalid.\n");
   }
#endif
   return d_convergence_factor[d_number_iterations-1];
}


template<int DIM> inline
void FACPreconditioner<DIM>::setPresmoothingSweeps( int num_pre_sweeps ) {
   d_presmoothing_sweeps = num_pre_sweeps;
   return;
}

template<int DIM> inline
void FACPreconditioner<DIM>::setPostsmoothingSweeps( int num_post_sweeps ) {
   d_postsmoothing_sweeps = num_post_sweeps;
   return;
}

template<int DIM> inline
void FACPreconditioner<DIM>::setMaxCycles( int max_cycles ) {
   d_max_iterations = max_cycles;
   return;
}

template<int DIM> inline
void FACPreconditioner<DIM>::setResidualTolerance( double residual_tol,
                                                   double relative_residual_tol ) {
   d_residual_tolerance = residual_tol;
   d_relative_residual_tolerance = relative_residual_tol;
   return;
}




template<int DIM> inline
void FACPreconditioner<DIM>::setAlgorithmChoice( const std::string &choice )
{
   std::string lower = choice;
   std::string::iterator c;
   for ( c=lower.begin(); c<lower.end(); ++c ) {
      *c = tolower(*c);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (    lower != "default"           /* Recursive from BTNG */
	&& lower != "mccormick-s4.3"	/* McCormick's secion 4.3 */
	&& lower != "pernice"           /* Translation of Pernice's */
        ) {
      TBOX_ERROR(d_object_name << ": algorithm should be set to one of\n"
		 << "'default' (recommended), 'mccormick-s4.3' or 'pernice'\n");
   }
#endif
   d_algorithm_choice = lower;
   return;
}



/*
*************************************************************************
*                                                                       *
* Functions to set logging parameters.                                  *
*                                                                       *
*************************************************************************
*/
template<int DIM> inline
void FACPreconditioner<DIM>::enableLogging(
   bool enabled
) {
   d_do_log = enabled;
   return;
}


}
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
