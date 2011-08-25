/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Implicit time integration manager class for nonlinear problems.
 *
 ************************************************************************/

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
double ImplicitIntegrator::getCurrentDt() const
{
   return d_current_dt;
}

SAMRAI_INLINE_KEYWORD
double ImplicitIntegrator::getInitialTime() const
{
   return d_initial_time;
}

SAMRAI_INLINE_KEYWORD
double ImplicitIntegrator::getFinalTime() const
{
   return d_final_time;
}

SAMRAI_INLINE_KEYWORD
double ImplicitIntegrator::getCurrentTime() const
{
   return d_current_time;
}

SAMRAI_INLINE_KEYWORD
bool ImplicitIntegrator::stepsRemaining() const
{
   return d_integrator_step < d_max_integrator_steps;
}

SAMRAI_INLINE_KEYWORD
int ImplicitIntegrator::getIntegratorStep() const
{
   return d_integrator_step;
}

SAMRAI_INLINE_KEYWORD
int ImplicitIntegrator::getMaxIntegratorSteps() const
{
   return d_max_integrator_steps;
}

SAMRAI_INLINE_KEYWORD
const std::string& ImplicitIntegrator::getObjectName() const
{
   return d_object_name;
}

}
}
