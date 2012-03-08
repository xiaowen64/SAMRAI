//
// File:        $URL$
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Integration routines for single level in AMR hierarchy
//              (basic hyperbolic systems).
//

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
const std::string&
HyperbolicLevelIntegrator::getObjectName() const
{
   return d_object_name;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
HyperbolicLevelIntegrator::getCurrentContext() const
{
   return d_current;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
HyperbolicLevelIntegrator::getNewContext() const
{
   return d_new;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
HyperbolicLevelIntegrator::getOldContext() const
{
   return d_old;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
HyperbolicLevelIntegrator::getScratchContext() const
{
   return d_scratch;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
HyperbolicLevelIntegrator::getPlotContext() const
{
   return d_plot_context;
}

SAMRAI_INLINE_KEYWORD
bool
HyperbolicLevelIntegrator::usingRefinedTimestepping() const
{
   return d_use_time_refinement;
}

}
}
