/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Utility class to read and write boxes to an HDF database.
 *
 ************************************************************************/

#ifndef included_hier_BoxIOUtility
#define included_hier_BoxIOUtility

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>

namespace SAMRAI {
namespace hier {

/**
 * Class BoxIOUtility supports writing and reading box information
 * to an HDF file.
 *
 * TODO: This class in used only in the application tests.  Should it
 * be moved there?
 */

class BoxIOUtility
{
public:
   /**
    * Enumerated type for specification of whether to read or write
    * data.
    *
    * - \b READ           { read from HDF database}
    * - \b WRITE          { write to HDF database};
    *
    */
   enum IOTYPE { READ = 0,
                 WRITE = 1 };
   /**
    * The constructor requires the name of the HDF database
    * to write or read to, and the IOTYPE.
    */
   BoxIOUtility(
      const std::string& filename,
      const IOTYPE iotype);

   /**
    * Destructor.
    */
   ~BoxIOUtility();

   /**
    * Pulls refinement boxes corresponding to the provided level and
    * entry number from storage array - returns
    * a boxlist with the corresponding refine boxes.
    */
   void
   getLevelBoxes(
      BoxContainer& level_boxes,
      const int level_number,
      const int entry_number);

   /**
    * Puts new refinement boxes corresponding to the provided level and
    * entry number into storage arrays.
    */
   void
   putLevelBoxes(
      const BoxContainer& level_boxes,
      const int level_number,
      const int entry_number);

   /**
    * Returns the number of levels in the database.
    */
   int
   getNumberOfLevels();

   /**
    * Returns the number of entries in the database for the specified
    * level.
    */
   int
   getNumberOfEntries(
      const int level_number);

   /**
    * Opens and writes to an HDF database directory with the prescribed
    * name a set of refinement boxes used during the run.
    */
   void
   writeLevelBoxesDatabase();

   /**
    * Print the boxes stored in the database to the specified IO stream.
    */
   void
   printBoxes(
      std::ostream& os);

private:
   /*
    * Opens and reads from an HDF database with the prescribed directory
    * the set of refine boxes for controlled gridding operations.
    */
   void
   readLevelBoxesDatabase();

   /*
    * HDF database name.
    */
   std::string d_hdf_filename;

   /*
    * The IO type - read or write
    */
   IOTYPE d_iotype;

   /*
    * tbox::Array to store the level boxes.
    */
   tbox::Array<tbox::Array<BoxContainer> > d_level_boxes;

};

}
}

#endif
