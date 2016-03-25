//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxIOUtility.h $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Utility class to read and write boxes to an HDF database.
//

#ifndef included_hier_BoxIOUtility
#define included_hier_BoxIOUtility

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_hier_BoxArray
#include "BoxArray.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <std::string>
#define included_String
#endif


namespace SAMRAI {
   namespace hier {

/**
 * Class BoxIOUtility<DIM> supports writing and reading box information
 * to an HDF file. 
 */

template<int DIM> class BoxIOUtility 
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
   enum IOTYPE{ READ = 0,
                WRITE = 1 };
   /**
    * The constructor requires the name of the HDF database
    * to write or read to, and the IOTYPE.
    */
   BoxIOUtility(const std::string& dirname,
                      const IOTYPE iotype);

   /**
    * Virtual destructor.
    */
   virtual ~BoxIOUtility<DIM>();

   /**
    * Pulls refinement boxes corresponding to the provided level and 
    * entry number from storage array - returns
    * a boxlist with the corresponding refine boxes.
    */
   void getLevelBoxes(BoxArray<DIM>& level_boxes,
                      const int level_number,
                      const int entry_number);

   /**
    * Puts new refinement boxes corresponding to the provided level and 
    * entry number into storage arrays.
    */
   void putLevelBoxes(const BoxArray<DIM>& level_boxes,
                      const int level_number,
                      const int entry_number);

   /**
    * Returns the number of levels in the database. 
    */
   int getNumberOfLevels();

   /**
    * Returns the number of entries in the database for the specified
    * level.
    */
   int getNumberOfEntries(const int level_number);

   /**
    * Opens and writes to an HDF database directory with the prescribed
    * name a set of refinement boxes used during the run.
    */
   void writeLevelBoxesDatabase();

   /**
    * Print the boxes stored in the database to the specified IO stream.
    */
   void printBoxes(std::ostream& os);

private:

   /*
    * Opens and reads from an HDF database with the prescribed directory 
    * the set of refine boxes for controlled gridding operations. 
    */
   void readLevelBoxesDatabase();

   /*
    * HDF database name.
    */
   std::string d_hdf_dirname;

   /*
    * The IO type - read or write
    */
   IOTYPE d_iotype;

   /*
    * tbox::Array to store the level boxes.
    */
   tbox::Array< tbox::Array< BoxArray<DIM> > > d_level_boxes;

};


}
}

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxIOUtility.C"
#endif

#endif

