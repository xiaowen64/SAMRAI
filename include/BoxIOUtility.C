//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxIOUtility.C $
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Utility class to read and write boxes to an HDF database.
//

#ifndef included_hier_BoxIOUtility_C
#define included_hier_BoxIOUtility_C

#include "BoxIOUtility.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#ifdef HAVE_HDF5
#include "tbox/HDFDatabase.h"
#endif 
#include "tbox/List.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"


#define MESH_REFINE_BOX_IO_UTILITY (1)

namespace SAMRAI {
   namespace hier {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for BoxIOUtility<DIM>.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM>  BoxIOUtility<DIM>::BoxIOUtility(
   const std::string& dirname, 
   const IOTYPE iotype) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!dirname.empty());
#endif

#ifndef HAVE_HDF5
    TBOX_ERROR("BoxIOUtility<DIM> constructor error"
               << "\n HDF Library not included in SAMRAI configure - "
               << "cannot read or write box information" << std::endl);
#endif


    d_hdf_dirname = dirname;
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
template<int DIM>  BoxIOUtility<DIM>::~BoxIOUtility()
{
   if (d_iotype == WRITE) {
      writeLevelBoxesDatabase();
   }
}

/*
*************************************************************************
*                                                                       *
* Read information from the storage arrays.                             *
*                                                                       *
*************************************************************************
*/
template<int DIM> void BoxIOUtility<DIM>::getLevelBoxes(
   BoxArray<DIM>& level_boxes,
   const int level_number,
   const int entry_number)
{
   tbox::plog << "Reading boxes for level: " << level_number
        << " entry number: " << entry_number << std::endl;

   /*
    * Make sure the data we are reading valid data from the database.
    */
   if (level_number+1 > d_level_boxes.getSize()) {
      TBOX_ERROR("BoxIOUtility::getLevelBoxes() error:  invalid level number. "
                 << "\n The level boxes database holds data only up "
                 << "\n to level " << d_level_boxes.getSize()-1 
                 << "\n You requested data from level " << level_number
                 << std::endl);
   }


   if (entry_number+1 > d_level_boxes[level_number].getSize()) {
      TBOX_ERROR("BoxIOUtility::getLevelBoxes() error:  invalid entry number. "
                 << "\n The level boxes database holds entries up "
                 << "\n to size " << d_level_boxes[level_number].getSize()
                 << "\n You requested data for entry " << entry_number
                 << std::endl);      
   } 

   tbox::plog << "Returning BoxArray containing " 
        << d_level_boxes[level_number][entry_number].getNumberOfBoxes() 
        << " boxes. " << std::endl;
   
   /*
    * Return BoxArray entry for this step.
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
template<int DIM> void BoxIOUtility<DIM>::putLevelBoxes(
   const BoxArray<DIM>& level_boxes, 
   const int level_number,
   const int entry_number)
{  

   tbox::plog << "Writing boxes for level: " << level_number
        << " entry number: " << entry_number << std::endl;

   /*
    * Reset dimension of arrays, if necessary.
    */
   if ( d_level_boxes.getSize() < level_number+1) {
      d_level_boxes.resizeArray(level_number+1);
   }
   
   /* 
    * Increment entry size by 10 to buffer the resize call.
    */
   int level_entry_size = d_level_boxes[level_number].getSize();
   if (level_entry_size < entry_number+1) {
      d_level_boxes[level_number].resizeArray(level_entry_size+10);
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
template<int DIM> int BoxIOUtility<DIM>::getNumberOfLevels()
{
   return(d_level_boxes.getSize());
}

/*
*************************************************************************
*                                                                       *
* Return number of entries for a particular level.                      *
*                                                                       *
*************************************************************************
*/
template<int DIM> int BoxIOUtility<DIM>::getNumberOfEntries(
   const int level_number)
{
   if (level_number+1 > d_level_boxes.getSize()) {
      TBOX_ERROR("BoxIOUtility::getNumberOfEntries() error: " 
                 << "invalid level number. "
                 << "\n The level boxes database holds data only up "
                 << "\n to level " << d_level_boxes.getSize()-1 
                 << "\n You requested data from level " << level_number
                 << std::endl);
   }

   return(d_level_boxes[level_number].getSize());
}


/*
*************************************************************************
*                                                                       *
* Opens an HDF5 database which we will either write to or read from.    *
*                                                                       *
*************************************************************************
*/
template<int DIM> void BoxIOUtility<DIM>::readLevelBoxesDatabase() 
{

#ifdef HAVE_HDF5
   /*
    * Open the HDF5 database.
    */
   tbox::Pointer<tbox::HDFDatabase> db = new tbox::HDFDatabase("root");
   int stat = db->mount(d_hdf_dirname,"R");
   if (stat < 0) {
     TBOX_ERROR("BoxIOUtility<DIM>::readLevelBoxesDatabase() error: "
                "\n Error opening HDF database: " << d_hdf_dirname << std::endl);
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
    *          BoxArray[i]    - box arrays for each of the entries
    */     

  for (int ln = 0; ln < num_levels; ln++) {
      
      /*
       * Read number of boxes for each entry of the level.
       */
      char *buffer1 = new char[16];
      sprintf(buffer1, "nboxes[%d]", ln);
      std::string s1(buffer1);
      delete [] buffer1;
      tbox::Array<int> number_of_boxes = db->getIntegerArray(s1);

      /*
       * Read box array for each entry of the level
       */
      int nentries = number_of_boxes.getSize();
      d_level_boxes[ln].resizeArray(nentries);
      for (i = 0; i < nentries; i++) {
         char *buffer2 = new char[16];
         sprintf(buffer2, "BoxArray[%d][%d]", ln, i);
         std::string s2(buffer2);
         delete [] buffer2;
         if (number_of_boxes[i] > 0) {
            d_level_boxes[ln][i] = db->getDatabaseBoxArray(s2);
         }
      }
  }
#endif
}

/*
*************************************************************************
*                                                                       *
* Writes level box information to HDF5 database                         *
*                                                                       *
*************************************************************************
*/
template<int DIM> void BoxIOUtility<DIM>::writeLevelBoxesDatabase()
{       
#ifdef HAVE_HDF5
   /*
    * Open the HDF5 database.
    */
   tbox::Pointer<tbox::HDFDatabase> db = new tbox::HDFDatabase("root");
   int stat = db->mount(d_hdf_dirname,"W");
   if (stat < 0) {
     TBOX_ERROR("BoxIOUtility<DIM>::writeLevelBoxesDatabase() error" 
                << "\n Error opening HDF database: " << d_hdf_dirname << std::endl);
   }

   /*
    * Cycle through the levels and write the number of entries
    * followed by the box array.
    *
    * Format:  num_levels     - number of levels in problem
    *          nboxes[i]      - number of boxes for each entry of the level
    *          BoxArray[i]    - box arrays for each of the entries
    */          
     
   /*
    * Write number of levels.
    */
   int num_levels =  d_level_boxes.getSize();
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
      char *buffer1 = new char[16];
      sprintf(buffer1, "nboxes[%d]", ln);
      std::string s1(buffer1);
      delete buffer1;
      db->putIntegerArray(s1, number_of_boxes );

      /*
       * Write box array[i]
       */
      for (i = 0; i < nentries; i++) {
         char *buffer2 = new char[16];
         sprintf(buffer2, "BoxArray[%d][%d]", ln, i);
         std::string s2(buffer2);
         delete buffer2;
         if (number_of_boxes[i] > 0) {
            db->putDatabaseBoxArray(s2, d_level_boxes[ln][i]);
         }
      }
   }
   db->unmount();

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
template<int DIM> void BoxIOUtility<DIM>::printBoxes(std::ostream& os)
{  
   os << "\n\n---------------------------------------------" << std::endl;
   os << "Boxes stored in database:"  << std::endl;

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
