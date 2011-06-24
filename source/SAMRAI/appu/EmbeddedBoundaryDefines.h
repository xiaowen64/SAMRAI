/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Definitions of cells and nodes on embedded boundaries. 
 *
 ************************************************************************/

#ifndef included_appu_EmbeddedBoundaryDefines
#define included_appu_EmbeddedBoundaryDefines

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace appu {

/*!
 * @brief Class EmbeddedBoundaryDefines sets the enumerated types used
 * to define cells and nodes on the embedded boundary level.
 *
 * @see appu::EmbeddedBoundaryGeometry
 * @see appu::CutCell
 * @see appu::BoundaryNode
 */

class EmbeddedBoundaryDefines
{
public:
   /*!
    * Enumerated type for the different cell classifications.
    *
    * - \b SOLID          {Cell is located in the "solid" region.}
    * - \b CUT            {Cell is cut, meaning a CutCell data
    *                      structure will be maintained at this cell.}
    * - \b BORDER         {Cell neighbors a cut cell, in the "flow" region.}
    * - \b FLOW           {Cell is located in the "flow" region.}
    */
   enum CELL_TYPE { SOLID = 0,
                    CUT = 1,
                    BORDER = 2,
                    FLOW = 3 };

   /*!
    * Enumerated type for inside/outside node classification.
    *
    * - \b INSIDE         {Node is located "inside" the prescribed geometry.}
    * - \b OUTSIDE        {Node is outside the geometry.}
    * - \b BOUNDARY       {Node is on the boundary of the geometry.  That is
    *                      it is the first one "inside" the geometry.}
    * - \b ONBOUNDARY     {Node is located exactly on the boundary of the
    *                      geometry (used to avoid divide-by-zero problems)
    *                      in numerical operations at embedded boundary.}
    */
   enum NODE_TYPE { OUTSIDE = 0,
                    INSIDE = 1,
                    BOUNDARY = 2,
                    ONBOUNDARY = -9 };
};

}
}
#endif
