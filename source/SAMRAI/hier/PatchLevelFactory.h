/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch level objects
 *
 ************************************************************************/

#ifndef included_hier_PatchLevelFactory
#define included_hier_PatchLevelFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchFactory.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/DescribedClass.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Factory used to create new patch levels.
 *
 * New types of patch level objects can be introduced into SAMRAI by deriving
 * from PatchLevelFactory and re-defining allocate.
 *
 * @see hier::PatchLevel
 */
class PatchLevelFactory:public tbox::DescribedClass
{
public:
   /*!
    * @brief Construct a patch level factory object.
    */
   PatchLevelFactory();

   /*!
    * @brief Virtual destructor for patch level factory objects.
    */
   virtual ~PatchLevelFactory();

   /*!
    * @brief Allocate a patch level with the specified boxes and processor mappings.
    *
    * Redefine this function to change the method for creating patch levels.
    *
    * @return A Pointer to the newly created PatchLevel.
    *
    * @param[in]  mapped_box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  factory @b Default: a Pointer to the standard PatchFactory
    */
   virtual tbox::Pointer<PatchLevel>
   allocate(
      const MappedBoxLevel& mapped_box_level,
      const tbox::Pointer<GridGeometry> grid_geometry,
      const tbox::Pointer<PatchDescriptor> descriptor,
      tbox::Pointer<PatchFactory> factory =
         tbox::Pointer<PatchFactory>(NULL)) const;

   /*!
    * @brief Allocate a patch level using the data from the database to
    * initialize it.
    *
    * The component_selector argument is used to specify which patch data
    * components to allocate and read in from the database.
    * @note
    * If desired, pass a ComponentSelector with all bits set to false to
    * indicate that no patch data components are read/allocated.
    *
    * Redefine this function to change the method for creating
    * patch levels from a database.
    *
    * @return A Pointer to the newly created PatchLevel.
    *
    * @param[in]  database
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  component_selector
    * @param[in]  factory @b Default: a Pointer to the standard PatchFactory
    * @param[in]  defer_boundary_box_creation @b Default: false
    */
   virtual tbox::Pointer<PatchLevel>
   allocate(
      tbox::Pointer<tbox::Database> database,
      const tbox::Pointer<GridGeometry> grid_geometry,
      const tbox::Pointer<PatchDescriptor> descriptor,
      const ComponentSelector& component_selector,
      tbox::Pointer<PatchFactory> factory = tbox::Pointer<PatchFactory>(NULL),
      const bool defer_boundary_box_creation = false) const;

private:
   /*
    * Copy constructor and assignment are not implemented.
    */
   PatchLevelFactory(
      const PatchLevelFactory&);
   void
   operator = (
      const PatchLevelFactory&);

};

}
}
#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/PatchLevelFactory.I"
#endif
#endif
