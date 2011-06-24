//
// File:        $URL$
// Package:     SAMRAI
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2010 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Construction and management of embedded boundary data
//              on an AMR hierarchy.
//
//

namespace SAMRAI {
namespace appu {

SAMRAI_INLINE_KEYWORD
const std::string& EmbeddedBoundaryGeometry::getObjectName() const
{
   return d_object_name;
}

} // namespace appu
} // namespace SAMRAI
