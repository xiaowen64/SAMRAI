/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   RefineSchedule's implementation of PatchHierarchy
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 ***********************************************************************
 * Register the static ConnectorWidthRequestor to be
 * automatically used by PatchHierarchy objects.
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineScheduleConnectorWidthRequestor::initializeCallback()
{
   hier::PatchHierarchy::registerAutoConnectorWidthRequestorStrategy(
      s_auto_registered_connector_width_requestor);
}

/*
 ***************************************************************************
 ***************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineScheduleConnectorWidthRequestor::finalizeCallback()
{
}

}
}
