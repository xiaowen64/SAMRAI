/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Blueprint utilities.
 *
 ************************************************************************/
#ifndef included_hier_BlueprintUtils
#define included_hier_BlueprintUtils

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlueprintUtilsStrategy.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for performing simple common tasks on a container
 * of Boxes.

Combined with a BlueprintUtilsStrategy, this loops over a hierarchy,
fills Blueprint topology and coordinate entries, calls back to user code
for specific coord choices:  uniform, rectilinear, explicit
 */

class PatchHierarchy;

class BlueprintUtils
{

public:

   BlueprintUtils(BlueprintUtilsStrategy* strategy);

   virtual ~BlueprintUtils();

   void putTopologyAndCoordinatesToDatabase(
      const std::shared_ptr<tbox::Database>& database,
      const PatchHierarchy& hierarchy) const;

private:

   BlueprintUtilsStrategy* d_strategy;

};

}
}

#endif  // included_hier_BlueprintUtils
