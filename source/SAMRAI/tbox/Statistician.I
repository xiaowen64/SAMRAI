//
// File:        $URL$
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
// Revision:    \f$       \f$
// Modified:    $LastChangedDate$
// Description: Singleton manager class for statistic objects.
//

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
const std::string& StatisticRestartDatabase::getObjectName() const
{
   return d_object_name;
}

}
}
