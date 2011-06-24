/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Utility class to read and write boxes to an HDF database. 
 *
 ************************************************************************/

#ifndef included_hier_BoxIOUtility_C
#define included_hier_BoxIOUtility_C

#include "SAMRAI/hier/BoxIOUtility.h"

#ifdef HAVE_HDF5
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#endif
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <fstream>

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *                                                                       *
 * Constructor and destructor for BoxIOUtility.                    *
 *                                                                       *
 *************************************************************************
 */

BoxIOUtility::BoxIOUtility(
   const std::string& filename,
   const IOTYPE iotype)
{
   TBOX_ASSERT(!filename.empty());

#ifndef HAVE_HDF5
   TBOX_ERROR("BoxIOUtility constructor error"
      << "\n HDF Library not included in SAMRAI configure - "
      << "cannot read or write box information" << std::endl);
#endif

   d_hdf_filename = filename;
   d_iotype = iotype;
   if (d_iotype == READ) {
      readLevelBoxesDatabase();
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Destructor writes refine boxes to HDF database.                       *
 *                                                                       *
 *************************************************************************
 */
BoxIOUtility::~BoxIOUtility()
{
}

/*
 *************************************************************************
 *                                                                       *
 * Read information from the storage arrays.                             *
 *                                                                       *
 *************************************************************************
 */
void BoxIOUtility::getLevelBoxes(
   BoxList& level_boxes,
   const int level_number,
   const int entry_number)
{
   tbox::plog << "Reading boxes for level: " << level_number
              << " entry number: " << entry_number << std::endl;

   /*
    * Make sure the data we are reading valid data from the database.
    */
   if (level_number + 1 > d_level_boxes.getSize()) {
      TBOX_ERROR("BoxIOUtility::getLevelBoxes() error:  invalid level number. "
         << "\n The level boxes database holds data only up "
         << "\n to level " << d_level_boxes.getSize() - 1
         << "\n You requested data from level " << level_number
         << std::endl);
   }

   if (entry_number + 1 > d_level_boxes[level_number].getSize()) {
      TBOX_ERROR("BoxIOUtility::getLevelBoxes() error:  invalid entry number. "
         << "\n The level boxes database holds entries up "
         << "\n to size " << d_level_boxes[level_number].getSize()
         << "\n You requested data for entry " << entry_number
         << std::endl);
   }

   tbox::plog << "Returning BoxList containing "
   << d_level_boxes[level_number][entry_number].getNumberOfBoxes()
   << " boxes. " << std::endl;

   /*
    * Return BoxList entry for this step.
    */
   level_boxes = d_level_boxes[level_number][entry_number];
}

/*
 *************************************************************************
 *                                                                       *
 * Pack information into storage arrays.                                 *
 *                                                                       *
 *************************************************************************
 */
void BoxIOUtility::putLevelBoxes(
   const BoxList& level_boxes,
   const int level_number,
   const int entry_number)
{
   const tbox::Dimension dim(level_boxes.getDim());

   tbox::plog << "Writing boxes for level: " << level_number
              << " entry number: " << entry_number << std::endl;

   /*
    * Reset dimension of arrays, if necessary.
    */
   if (d_level_boxes.getSize() < level_number + 1) {
      d_level_boxes.resizeArray(level_number + 1);
   }

   /*
    * Increment entry size by 10 to buffer the resize call.
    */
   int level_entry_size = d_level_boxes[level_number].getSize();
   if (level_entry_size < entry_number + 1) {
      d_level_boxes[level_number].resizeArray(level_entry_size + 10,
         hier::BoxList(dim));
   }

   /*
    * Pack array
    */
   d_level_boxes[level_number][entry_number] = level_boxes;
}

/*
 *************************************************************************
 *                                                                       *
 * Return number of levels in the Database.                              *
 *                                                                       *
 *************************************************************************
 */
int BoxIOUtility::getNumberOfLevels()
{
   return d_level_boxes.getSize();
}

/*
 *************************************************************************
 *                                                                       *
 * Return number of entries for a particular level.                      *
 *                                                                       *
 *************************************************************************
 */
int BoxIOUtility::getNumberOfEntries(
   const int level_number)
{
   if (level_number + 1 > d_level_boxes.getSize()) {
      TBOX_ERROR("BoxIOUtility::getNumberOfEntries() error: "
         << "invalid level number. "
         << "\n The level boxes database holds data only up "
         << "\n to level " << d_level_boxes.getSize() - 1
         << "\n You requested data from level " << level_number
         << std::endl);
   }

   return d_level_boxes[level_number].getSize();
}

/*
 *************************************************************************
 *                                                                       *
 * Opens an HDF5 database which we will either write to or read from.    *
 *                                                                       *
 *************************************************************************
 */
void BoxIOUtility::readLevelBoxesDatabase()
{

#ifdef HAVE_HDF5
   /*
    * Open the HDF5 database.
    */
   tbox::Pointer<tbox::HDFDatabase> db(new tbox::HDFDatabase("root"));
   int stat = db->open(d_hdf_filename);
   if (stat < 0) {
      TBOX_ERROR("BoxIOUtility::readLevelBoxesDatabase() error: "
         "\n Error opening HDF database: " << d_hdf_filename << std::endl);
   }

   /*
    * Read number of levels.
    */
   int num_levels = db->getInteger("num_levels");

   d_level_boxes.resizeArray(num_levels);

   int i;

   /*
    * Cycle through the levels and read the number of entries
    * for each level, followed by the box arrays.
    *
    * Format:  nboxes[i]      - number of boxes for each entry
    *          BoxList[i]     - box arrays for each of the entries
    */

   for (int ln = 0; ln < num_levels; ln++) {

      /*
       * Read number of boxes for each entry of the level.
       */
      std::string s1 = "nboxes[" + tbox::Utilities::intToString(ln) + "]";
      tbox::Array<int> number_of_boxes = db->getIntegerArray(s1);

      /*
       * Read box array for each entry of the level
       */
      int nentries = number_of_boxes.getSize();

      for (i = 0; i < nentries; i++) {
         std::string s2 = "BoxList[" + tbox::Utilities::intToString(ln) + "]["
            + tbox::Utilities::intToString(i) + "]";
         if (number_of_boxes[i] > 0) {
            tbox::Array<tbox::DatabaseBox> db_array = db->getDatabaseBoxArray(
                  s2);
            /*
             * Do not have DIM so we need to resize inside the loop after getting
             * a box.
             */
            if (d_level_boxes[ln].size() < nentries) {
               d_level_boxes[ln].resizeArray(nentries,
                  hier::BoxList(db_array[0].getDim()));
            }

            d_level_boxes[ln][i] = db_array;
         }
      }
   }

   db->close();
#endif
}

/*
 *************************************************************************
 *                                                                       *
 * Writes level box information to HDF5 database                         *
 *                                                                       *
 *************************************************************************
 */
void BoxIOUtility::writeLevelBoxesDatabase()
{
#ifdef HAVE_HDF5
   /*
    * Open the HDF5 database.
    */
   tbox::Pointer<tbox::Database> db(new tbox::HDFDatabase("root"));
   int stat = db->create(d_hdf_filename);
   if (stat < 0) {
      TBOX_ERROR("BoxIOUtility::writeLevelBoxesDatabase() error"
         << "\n Error opening HDF database: " << d_hdf_filename << std::endl);
   }

   /*
    * Cycle through the levels and write the number of entries
    * followed by the box array.
    *
    * Format:  num_levels     - number of levels in problem
    *          nboxes[i]      - number of boxes for each entry of the level
    *          BoxList[i]     - box arrays for each of the entries
    */

   /*
    * Write number of levels.
    */
   int num_levels = d_level_boxes.getSize();
   db->putInteger("num_levels", num_levels);

   int i;

   for (int ln = 0; ln < num_levels; ln++) {

      /*
       * Write nboxes[i].
       */
      int nentries = d_level_boxes[ln].getSize();
      tbox::Array<int> number_of_boxes(nentries);
      for (i = 0; i < nentries; i++) {
         number_of_boxes[i] = d_level_boxes[ln][i].getNumberOfBoxes();
      }
      std::string s1 = "nboxes[" + tbox::Utilities::intToString(ln) + "]";
      db->putIntegerArray(s1, number_of_boxes);

      /*
       * Write box array[i]
       */
      for (i = 0; i < nentries; i++) {
         std::string s2 = "BoxList[" + tbox::Utilities::intToString(ln)
            + "][" + tbox::Utilities::intToString(i) + "]";
         if (number_of_boxes[i] > 0) {
            db->putDatabaseBoxArray(s2, d_level_boxes[ln][i]);
         }
      }
   }
   db->close();

   d_level_boxes.resizeArray(0);
#endif
}

/*
 *************************************************************************
 *                                                                       *
 * Print the boxes stored om the database to the specified IO stream.    *
 *                                                                       *
 *************************************************************************
 */
void BoxIOUtility::printBoxes(
   std::ostream& os)
{
   os << "\n\n---------------------------------------------" << std::endl;
   os << "Boxes stored in database:" << std::endl;

   int nlevels = d_level_boxes.getSize();

   for (int ln = 0; ln < nlevels; ln++) {
      os << "Level " << ln << ":" << std::endl;
      os << "-------" << std::endl;

      for (int i = 0; i < d_level_boxes[ln].getSize(); i++) {

         os << "   Entry " << i << ": " << std::endl;
         d_level_boxes[ln][i].print(os);
      }
      os << "\n" << std::endl;
   }
   os << "---------------------------------------------\n\n" << std::endl;

}

}
}

#endif
