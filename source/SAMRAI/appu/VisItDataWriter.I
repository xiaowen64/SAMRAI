//
// File:        $URL$
// Package:     SAMRAI
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2011 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision$
// Modified:    $LastChangedDate$
// Description: Writes data files for visualization by VisIt
//
//

namespace SAMRAI {
namespace appu {

SAMRAI_INLINE_KEYWORD
const std::string& VisItDataWriter::getObjectName() const
{
   return d_object_name;
}

} // namespace appu
} // namespace SAMRAI
