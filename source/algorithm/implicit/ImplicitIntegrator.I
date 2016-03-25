//
// File:        ImplicitIntegrator.I
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 7 $
// Modified:    $Date: 2004-11-30 13:18:17 -0800 (Tue, 30 Nov 2004) $
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
