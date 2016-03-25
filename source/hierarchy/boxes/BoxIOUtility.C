//
// File:        BoxIOUtility.C
// Package:     SAMRAI hierarchy
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Utility class to read and write boxes to an HDF database.
//

#ifndef included_hier_BoxIOUtility_C
#define included_hier_BoxIOUtility_C

#include "BoxIOUtility.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

#ifdef HAVE_HDF5
#include "tbox/HDFDatabase.h"
#endif 
#include "tbox/List.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

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
   const string& dirname, 
   const IOTYPE iotype) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!dirname.empty());
#endif

#ifndef HAVE_HDF5
    TBOX_ERROR("BoxIOUtility<DIM> constructor error"
               << "\n HDF Library not included in SAMRAI configure - "
               << "cannot read or write box information" << endl);
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
        << " entry number: " << entry_number << endl;

   /*
    * Make sure the data we are reading valid data from the database.
    */
   if (level_number+1 > d_level_boxes.getSize()) {
      TBOX_ERROR("BoxIOUtility::getLevelBoxes() error:  invalid level number. "
                 << "\n The level boxes database holds data only up "
                 << "\n to level " << d_level_boxes.getSize()-1 
                 << "\n You requested data from level " << level_number
                 << endl);
   }


   if (entry_number+1 > d_level_boxes[level_number].getSize()) {
      TBOX_ERROR("BoxIOUtility::getLevelBoxes() error:  invalid entry number. "
                 << "\n The level boxes database holds entries up "
                 << "\n to size " << d_level_boxes[level_number].getSize()
                 << "\n You requested data for entry " << entry_number
                 << endl);      
   } 

   tbox::plog << "Returning BoxArray containing " 
        << d_level_boxes[level_number][entry_number].getNumberOfBoxes() 
        << " boxes. " << endl;
   
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
        << " entry number: " << entry_number << endl;

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
                 << endl);
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
                "\n Error opening HDF database: " << d_hdf_dirname << endl);
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
      string s1(buffer1);
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
         string s2(buffer2);
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
                << "\n Error opening HDF database: " << d_hdf_dirname << endl);
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
      string s1(buffer1);
      delete buffer1;
      db->putIntegerArray(s1, number_of_boxes );

      /*
       * Write box array[i]
       */
      for (i = 0; i < nentries; i++) {
         char *buffer2 = new char[16];
         sprintf(buffer2, "BoxArray[%d][%d]", ln, i);
         string s2(buffer2);
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
template<int DIM> void BoxIOUtility<DIM>::printBoxes(ostream& os)
{  
   os << "\n\n---------------------------------------------" << endl;
   os << "Boxes stored in database:"  << endl;

   int nlevels = d_level_boxes.getSize();
   
   for (int ln = 0; ln < nlevels; ln++) {
      os << "Level " << ln << ":" << endl;
      os << "-------" << endl;

      for (int i = 0; i < d_level_boxes[ln].getSize(); i++) {

         os << "   Entry " << i << ": " << endl;
         d_level_boxes[ln][i].print(os);
      }
      os << "\n" << endl;
   }
   os << "---------------------------------------------\n\n" << endl;

}
       
}
}

#endif
