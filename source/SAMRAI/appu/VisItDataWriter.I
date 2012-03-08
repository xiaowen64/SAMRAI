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
const std::string&
VisItDataWriter::getObjectName() const
{
   return d_object_name;
}

/*
 *************************************************************************
 *
 * Set default derived data writer.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
VisItDataWriter::setDefaultDerivedDataWriter(
   VisDerivedDataStrategy* derived_writer)
{
   TBOX_ASSERT(derived_writer != (VisDerivedDataStrategy *)NULL);

   d_default_derived_writer = derived_writer;
}

/*
 *************************************************************************
 *
 * Set materials data writer object.
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
VisItDataWriter::setMaterialsDataWriter(
   VisMaterialsDataStrategy* materials_writer)
{
   TBOX_ASSERT(materials_writer != (VisMaterialsDataStrategy *)NULL);

   d_materials_writer = materials_writer;
}

/*
 *************************************************************************
 *
 * Write plot data from given hierarchy to HDF file
 *
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
VisItDataWriter::setSummaryFilename(
   std::string& filename)
{
   TBOX_ASSERT(!filename.empty());

   d_summary_filename = filename + ".samrai";
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
VisItDataWriter::initializeCallback()
{
   t_write_plot_data =
      tbox::TimerManager::getManager()->getTimer("appu:VisItDataWriter::writePlotData()");
}

/*
 *************************************************************************
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
VisItDataWriter::finalizeCallback()
{
   t_write_plot_data.reset();
}

}
}
