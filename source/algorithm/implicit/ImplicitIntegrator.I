//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/implicit/ImplicitIntegrator.I $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Implicit time integration manager class for nonlinear problems.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace algs {

template<int DIM> inline
double ImplicitIntegrator<DIM>::getCurrentDt() const
{
   return(d_current_dt);
}

template<int DIM> inline
double ImplicitIntegrator<DIM>::getInitialTime() const
{
   return(d_initial_time);
}

template<int DIM> inline
double ImplicitIntegrator<DIM>::getFinalTime() const
{
   return(d_final_time);
}

template<int DIM> inline
double ImplicitIntegrator<DIM>::getCurrentTime() const
{
   return (d_current_time);
}

template<int DIM> inline
bool ImplicitIntegrator<DIM>::stepsRemaining() const
{
   return(d_integrator_step < d_max_integrator_steps);
}

template<int DIM> inline 
int ImplicitIntegrator<DIM>::getIntegratorStep() const
{
   return(d_integrator_step);
}

template<int DIM> inline 
int ImplicitIntegrator<DIM>::getMaxIntegratorSteps() const
{
   return(d_max_integrator_steps);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
