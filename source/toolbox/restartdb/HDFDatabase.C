//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/toolbox/restartdb/HDFDatabase.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1846 $
// Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description: A database structure that stores HDF5 format data.
//

#include "tbox/HDFDatabase.h"

#ifdef HAVE_HDF5

#include "tbox/IOStream.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


/*
*************************************************************************
*                                                                       *
* Integer keys for identifying types in HDF5 database.  Negative        *
* entries are used to distinguish arrays from scalars when printing     *
* key information.                                                      *
*                                                                       *
*************************************************************************
*/
#define KEY_DATABASE        (0)
#define KEY_BOOL_ARRAY      (1)
#define KEY_BOX_ARRAY       (2)
#define KEY_CHAR_ARRAY      (3)
#define KEY_COMPLEX_ARRAY   (4)
#define KEY_DOUBLE_ARRAY    (5)
#define KEY_FLOAT_ARRAY     (6)
#define KEY_INT_ARRAY       (7)
#define KEY_STRING_ARRAY    (8)

#define KEY_BOOL_SCALAR     (-1)
#define KEY_BOX_SCALAR      (-2)
#define KEY_CHAR_SCALAR     (-3)
#define KEY_COMPLEX_SCALAR  (-4)
#define KEY_DOUBLE_SCALAR   (-5)
#define KEY_FLOAT_SCALAR    (-6)
#define KEY_INT_SCALAR      (-7)
#define KEY_STRING_SCALAR   (-8)


/*
  Macros starting with H5T_SAMRAI_ are for controlling the data
  type that is actually written to the file.  As long as
  these are not "native" types, the file should be portable.
*/

// Type used for writing simple (non-compound) data.
#define H5T_SAMRAI_INT      H5T_STD_I32BE
#define H5T_SAMRAI_FLOAT    H5T_IEEE_F32BE
#define H5T_SAMRAI_DOUBLE   H5T_IEEE_F64BE
#define H5T_SAMRAI_BOOL     H5T_STD_I8BE

// Type used for writing the data attribute key.
#define H5T_SAMRAI_ATTR H5T_STD_I8BE


/*
*************************************************************************
*                                                                       *
* Macros to suppress the HDF5 messages sent to standard i/o; handle     *
* errors explicity within this code.                                    *
*                                                                       *
*************************************************************************
*/

#define BEGIN_SUPPRESS_HDF5_WARNINGS                  \
{                                                     \
   herr_t (*H5E_saved_efunc) (void*) = NULL;          \
   void *H5E_saved_edata = NULL;                      \
   H5Eget_auto(&H5E_saved_efunc, &H5E_saved_edata);   \
   H5Eset_auto(NULL, NULL);                              

#define END_SUPPRESS_HDF5_WARNINGS                     \
   H5Eset_auto(H5E_saved_efunc, H5E_saved_edata);      \
}



/*
*************************************************************************
* We may wish to assert HDF5 return values regardless of debug modes.   *
*************************************************************************
*/
#define ASSERT_HDF5_RETURN_VALUES
#ifdef ASSERT_HDF5_RETURN_VALUES
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif



namespace SAMRAI {
   namespace tbox {

std::string HDFDatabase::s_top_level_search_group = std::string();
std::string HDFDatabase::s_group_to_search = std::string();
int HDFDatabase::s_still_searching = 0;
int HDFDatabase::s_found_group = 0;

/*
*************************************************************************
*                                                                       *
* Static member function to iterate through the hdf5 data file and      *
* assemble a list of desired (key, type) pairs.                         *
*                                                                       *
*************************************************************************
*/

herr_t HDFDatabase::iterateKeys(
   hid_t loc_id,
   const char *name,
   void *opdata)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(name != (char*)NULL);
#endif

   if (s_still_searching) {

     H5G_stat_t statbuf;
     int type_key;
     herr_t errf;

     errf = H5Gget_objinfo(loc_id, name, 0, &statbuf);
#ifdef ASSERT_HDF5_RETURN_VALUES
     TBOX_ASSERT( errf >= 0 );
#endif

     switch (statbuf.type) {
     case H5G_GROUP: {
       if (s_top_level_search_group == "/") {
         addKeyToList(name, KEY_DATABASE, opdata);
       } else if ( !strcmp(name, s_group_to_search.c_str()) ) {
         hid_t grp;
         grp = H5Gopen(loc_id, name);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( grp >= 0 );
#endif
         s_found_group = true;
         s_still_searching =
           H5Giterate(grp, ".", NULL,
                      HDFDatabase::iterateKeys, opdata);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( s_still_searching >= 0 );
#endif
         s_found_group = false;
       } else {
         hid_t grp;
         grp = H5Gopen(loc_id, name);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( grp >= 0 );
#endif
         if (s_found_group) {
           addKeyToList(name, KEY_DATABASE, opdata);
         } else {
           errf = H5Giterate(grp, ".", NULL,
                             HDFDatabase::iterateKeys, opdata);
#ifdef ASSERT_HDF5_RETURN_VALUES
           TBOX_ASSERT( errf >= 0 );
#endif
         }
       }
       break;
     }

     case H5G_DATASET: {
       if (s_still_searching && s_found_group) {
         hid_t this_set;
         BEGIN_SUPPRESS_HDF5_WARNINGS
         this_set = H5Dopen(loc_id, name);
         END_SUPPRESS_HDF5_WARNINGS
         if (this_set > 0) {
            hid_t attr = H5Aopen_name(this_set, "Type");
#ifdef ASSERT_HDF5_RETURN_VALUES
            TBOX_ASSERT( attr >= 0 );
#endif
            errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
#ifdef ASSERT_HDF5_RETURN_VALUES
            TBOX_ASSERT( errf >= 0 );
#endif
            hid_t this_space = H5Dget_space(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
            TBOX_ASSERT( this_space >= 0 );
#endif
            hsize_t nsel = H5Sget_select_npoints(this_space);
            int array_size = int(nsel); 
            addKeyToList(name,
                         (array_size == 1 ? -type_key : type_key),
                         opdata);
            errf = H5Sclose(this_space);
#ifdef ASSERT_HDF5_RETURN_VALUES
            TBOX_ASSERT( errf >= 0 );
#endif
            errf = H5Aclose(attr);
#ifdef ASSERT_HDF5_RETURN_VALUES
            TBOX_ASSERT( errf >= 0 );
#endif
            errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
            TBOX_ASSERT( errf >= 0 );
#endif
         }
       }
       break;
     }

     default: {
       TBOX_ERROR("HDFDatabase key search error....\n"
                  << "   Unable to identify key = " << name 
                  << " as a known group or dataset" << std::endl);
     }
     }

   }
   return 0;
}

/*
*************************************************************************
*                                                                       *
* Static member function to add key to list for database associated     *
* with void* argument.                                                  *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::addKeyToList(
   const char *name,
   int type,
   void* database) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(name != (char*)NULL);
   TBOX_ASSERT(database != NULL);
#endif

   KeyData key_item;
   key_item.d_key  = name; 
   key_item.d_type = type; 

   ((HDFDatabase*)database)->d_keydata.appendItem(key_item);
}

/*
*************************************************************************
*                                                                       *
* Public HDF database constructor creates an empty database with the    *
* specified name.  It sets the group_ID to a default value of -1.       *
* This data is used by member functions to track parent databases.      *
*                                                                       *
*************************************************************************
*/

HDFDatabase::HDFDatabase(const std::string& name) :
   d_is_file(false),
   d_file_id(-1),
   d_group_id(-1),
   d_database_name(name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
#endif
   d_keydata.clearItems(); 
}

/*
*************************************************************************
*                                                                       *
* Private HDF database constructor creates an empty database with the   *
* specified name.  The group_ID is used privately within                *
* the member functions to track parent databases.                       *
*                                                                       *
*************************************************************************
*/

HDFDatabase::HDFDatabase(
   const std::string& name, 
   hid_t group_ID) :
   d_is_file(false),
   d_file_id(-1),
   d_group_id(group_ID),
   d_database_name(name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
#endif
   d_keydata.clearItems(); 
}

/*
*************************************************************************
*                                                                       *
* The database destructor closes the opened file or group.              *
*                                                                       *
*************************************************************************
*/

HDFDatabase::~HDFDatabase()
{
   herr_t errf;

   NULL_USE(errf);

   if (d_is_file) {
      unmount();
   } 

   if ( d_group_id != -1 ) {
      errf = H5Gclose(d_group_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT( errf >= 0 );
#endif

   }

}

/*
*************************************************************************
*                                                                       *
* Return true if the key exists within the database; false otherwise.   *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::keyExists(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   bool key_exists = false;
   herr_t errf;
   
   hid_t this_set;
   BEGIN_SUPPRESS_HDF5_WARNINGS
   this_set = H5Dopen(d_group_id, key.c_str());
   END_SUPPRESS_HDF5_WARNINGS
   if (this_set > 0) {
      key_exists = true;
      errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }
   if (!key_exists) {
      hid_t this_group;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_group = H5Gopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_group > 0) {
         key_exists = true;
         errf = H5Gclose(this_group);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return key_exists;
}

/*
*************************************************************************
*                                                                       *
* Return all keys in the database.                                      *
*                                                                       *
*************************************************************************
*/

Array<std::string> HDFDatabase::getAllKeys()
{
   performKeySearch();

   Array<std::string> tmp_keys(d_keydata.getNumberItems());

   int k = 0;
   for (List<KeyData>::Iterator i(d_keydata); i; i++) {
      tmp_keys[k] = i().d_key; 
      k++;
   }

   cleanupKeySearch();

   return(tmp_keys);
}

/*
*************************************************************************
*                                                                       *
* Return the size of the array associated with the key.  If the key     *
* does not exist, then zero is returned.                                *
* Array size is set based on the number of elements (points) within     *
* the dataspace defined by the named dataset (or key).                  *
*                                                                       *
*************************************************************************
*/

int HDFDatabase::getArraySize(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   herr_t errf;
   int array_size  = 0;

   hid_t this_set;
   BEGIN_SUPPRESS_HDF5_WARNINGS
   this_set = H5Dopen(d_group_id, key.c_str());
   END_SUPPRESS_HDF5_WARNINGS
   if (this_set > 0) {
      hid_t this_space = H5Dget_space(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( this_space >= 0 );
#endif
      hsize_t nsel = H5Sget_select_npoints(this_space);
      array_size = int(nsel);
      errf = H5Sclose(this_space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   return array_size;
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a database entry.  If the key does not exist, then false   *
* is returned.  The key represents a database (or hdf group) if the     *
* H5Gopen function on the key is successful.                            *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isDatabase(const std::string& key)
{
   bool is_database = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_group;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_group = H5Gopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_group > 0) {
         is_database = true;
         errf = H5Gclose(this_group);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_database;
}

/*
*************************************************************************
*                                                                       *
* Create a new database with the specified key name.                    *
*                                                                       *
*************************************************************************
*/

Pointer<Database> 
HDFDatabase::putDatabase(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif

   std::string parent_name = d_database_name;
   hid_t  parent_id   = d_group_id;

   hid_t this_group = H5Gcreate(parent_id, key.c_str(), 0);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( this_group >= 0 );

#endif
   Pointer<Database> new_database = 
      new HDFDatabase(key, this_group);

   return(new_database);
}

/*
************************************************************************
*                                                                      *
* Get the database with the specified key name.                        *
*                                                                      *
************************************************************************
*/

Pointer<Database> 
HDFDatabase::getDatabase(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if (!isDatabase(key)) {
      TBOX_ERROR("HDFDatabase::getDatabase() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a database." << std::endl);
   }

   hid_t this_group = H5Gopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( this_group >= 0 );

#endif
   Pointer<Database> database = 
      new HDFDatabase(key, this_group);

   return(database);
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a boolean entry.  If the key does not exist, then false    *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isBool(const std::string& key)
{
   bool is_boolean  = false;
   herr_t errf;
   
   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_BOOL_ARRAY) {
            is_boolean = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_boolean;
}

/*
*************************************************************************
*                                                                       *
* Create a boolean scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putBool(
   const std::string& key, 
   const bool& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putBoolArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a boolean array entry in the database with the specified       *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putBoolArray(
   const std::string& key, 
   const Array<bool>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putBoolArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putBoolArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a boolean array entry in the database with the specified       *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putBoolArray(
   const std::string& key, 
   const bool* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (bool*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[1] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( space >= 0 );
#endif

      /*
        We cannot be sure exactly what bool is because it is
        represented differently on different platforms, and
        it may have been redefined, i.e., by the Boolean
        type.  We are unsure what the bool is so we convert it
        to the native int type (H5T_NATIVE_INT) before giving
        it to HDF.  When we write a bool, we write it the
        shortest integer type we can find, the H5T_SAMRAI_BOOL
        type.
      */
      Array<int> data1( nelements );
      for ( int i=0; i<nelements; ++i ) data1[i] = data[i];

      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_BOOL,
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &data1[0]);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_BOOL_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putBoolArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get boolean scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a boolean type.                                                 *
*                                                                      *
************************************************************************
*/

bool HDFDatabase::getBool(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   bool ret_val;
   getBoolArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get boolean scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a boolean type.                                                 *
*                                                                      *
************************************************************************
*/

bool HDFDatabase::getBoolWithDefault(
   const std::string& key, 
   const bool& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<bool> local_bool = getBoolArray(key);
   bool *locptr = local_bool.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get boolean arrays from the database with the        *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a boolean type.                   *
*                                                                      *
************************************************************************
*/

Array<bool> HDFDatabase::getBoolArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if (!isBool(key)) {
      TBOX_ERROR("HDFDatabase::getBoolArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a bool array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;
   herr_t errf;

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   Array<bool> bool_array(nsel);

   if (nsel > 0) {
      /*
        We cannot be sure exactly what bool is because it is
        represented differently on different platforms, and
        it may have been redefined, i.e., by the Boolean
        type.  So we read bools into native integer memory
        then convert.
      */
      Array<int> data1( nsel );
      int* locPtr = data1.getPointer();
      errf = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      // Convert what was just read in.
      for ( size_t i=0; i<nsel; ++i ) bool_array[i] = data1[i];
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );

#endif
   return bool_array;
}

void HDFDatabase::getBoolArray(
   const std::string& key,
   bool* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<bool> tmp = getBoolArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getBoolArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a box entry.  If the key does not exist, then false        *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isDatabaseBox(const std::string& key)
{
   bool is_box  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_BOX_ARRAY) {
            is_box = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_box;
}

/*
*************************************************************************
*                                                                       *
* Create a box entry in the database with the specified                 *
* key name.  A box entry is an array of one.                            *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDatabaseBox(
   const std::string& key, 
   const DatabaseBox& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putDatabaseBoxArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a box array entry in the database with the specified key name. *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDatabaseBoxArray(
   const std::string& key, 
   const Array<DatabaseBox>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putDatabaseBoxArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putDatabaseBoxArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a box array entry in the database with the specified key name. *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (DatabaseBox*)NULL);
#endif
   if (nelements > 0) {

      herr_t errf;

      // Memory type
      hid_t mtype = createCompoundDatabaseBox('n');
      // Storage type
      hid_t stype = createCompoundDatabaseBox('s');

      hsize_t length = nelements;
      hid_t space = H5Screate_simple(1, &length, NULL);

      hid_t dataset =
         H5Dcreate( d_group_id, key.c_str(), stype, space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );

#endif
      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_BOX_ARRAY, dataset );

      errf = H5Tclose(mtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Tclose(stype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );

#endif
   } else {
      TBOX_ERROR("HDFDatabase::putDatabaseBoxArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get box scalar entry from the database with the specified key        *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a box type.                                                     *
*                                                                      *
************************************************************************
*/

DatabaseBox HDFDatabase::getDatabaseBox(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   DatabaseBox ret_val;
   getDatabaseBoxArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get box scalar entry from the database with the specified key        *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a box type.                                                     *
*                                                                      *
************************************************************************
*/

DatabaseBox HDFDatabase::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<DatabaseBox> local_box = getDatabaseBoxArray(key);
   DatabaseBox *locptr = local_box.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get box arrays from the database with the            *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a box type.                       *
*                                                                      *
************************************************************************
*/

Array<DatabaseBox> HDFDatabase::getDatabaseBoxArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( !isDatabaseBox(key) ) {
      TBOX_ERROR("HDFDatabase::getDatabaseBoxArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a box array." << std::endl);
   }

   hid_t   dset, dspace;
   hsize_t nsel;
   herr_t errf;

   // Memory type
   hid_t mtype = createCompoundDatabaseBox('n');

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   Array<DatabaseBox> boxArray(nsel);

   if (nsel > 0) {
      DatabaseBox* locPtr = boxArray.getPointer();
      errf = H5Dread(dset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   errf = H5Tclose(mtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );

#endif
   return boxArray;
}

void HDFDatabase::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<DatabaseBox> tmp = getDatabaseBoxArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getDatabaseBoxArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

hid_t HDFDatabase::createCompoundDatabaseBox( char type_spec ) const {
   herr_t errf;
   hid_t int_type_spec=H5T_SAMRAI_INT;
   switch (type_spec) {
   case 'n':
      // Use native type specs.
      int_type_spec = H5T_NATIVE_INT;
      break;
   case 's':
      // Use storage type specs.
      int_type_spec = H5T_SAMRAI_INT;
      break;
   default:
      TBOX_ERROR("HDFDatabase::createCompundDatabaseBox() error in database "
         << d_database_name
         << "\n    Unknown type specifier found. " << std::endl);
   }
   hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(DatabaseBox));
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT(type >= 0);
#endif
   errf = H5Tinsert(type, "dim", HOFFSET(DatabaseBox_POD,d_dimension),
                    int_type_spec);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0);
#endif
   const hsize_t box_dim = DatabaseBox_MAX_DIM /* defined in DatabaseBox.h */;
   insertArray(type, "lo", HOFFSET(DatabaseBox_POD,d_lo), 1, &box_dim,
               NULL, int_type_spec);
   insertArray(type, "hi", HOFFSET(DatabaseBox_POD,d_hi), 1, &box_dim,
               NULL, int_type_spec);
   return type;
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a char entry.  If the key does not exist, then false       *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isChar(const std::string& key)
{
   bool is_char  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_CHAR_ARRAY) {
            is_char = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_char;
}

/*
*************************************************************************
*                                                                       *
* Create a char scalar entry in the database with the specified         *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putChar(
   const std::string& key, 
   const char& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putCharArray(key, &data, 1);

}

/*
*************************************************************************
*                                                                       *
* Create a char array entry in the database with the specified          *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putCharArray(
   const std::string& key, 
   const Array<char>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putCharArray(key, data.getPointer(), data.getSize());
   } else { 
      TBOX_ERROR("HDFDatabase::putCharArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a char array entry in the database with the specified          *
* key name. The charentry is defined by the hdf type H5T_C_S1.          *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putCharArray(
   const std::string& key,
   const char* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (char*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hid_t atype, space, dataset;

      char* local_buf = new char[nelements];

      atype = H5Tcopy(H5T_C_S1);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( atype >= 0 );
#endif
      errf = H5Tset_size(atype, nelements);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Tset_strpad(atype, H5T_STR_NULLTERM);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      for (int i = 0; i < nelements; i++) {
         local_buf[i] = data[i];
      }

      space = H5Screate(H5S_SCALAR);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( space >= 0 );
#endif

      dataset = 
         H5Dcreate(d_group_id, key.c_str(), atype, space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif

      errf = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, local_buf);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_CHAR_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Tclose(atype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      delete[] local_buf;

   } else {
      TBOX_ERROR("HDFDatabase::putCharArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get char scalar entry from the database with the specified key       *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a char type.                                                    *
*                                                                      *
************************************************************************
*/

char HDFDatabase::getChar(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   char ret_val;
   getCharArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get char scalar entry from the database with the specified key       *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a char type.                                                    *
*                                                                      *
************************************************************************
*/

char HDFDatabase::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<char> local_char = getCharArray(key);
   char *locptr = local_char.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get char arrays from the database with the           *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a char type.                      *
*                                                                      *
************************************************************************
*/

Array<char> HDFDatabase::getCharArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( !isChar(key) ) {
      TBOX_ERROR("HDFDatabase::getCharArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a char array." << std::endl);
   } 

   hid_t   dset, dspace, dtype;
   size_t  nsel = 0;
   herr_t errf;

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   dtype  = H5Dget_type(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dtype >= 0 );
#endif
   nsel   = H5Tget_size(dtype);

   Array<char> charArray(nsel);

   if (nsel > 0) {
      char* locPtr = charArray.getPointer();
      errf = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Tclose(dtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif

   return charArray;
}

void HDFDatabase::getCharArray(
   const std::string& key,
   char* data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<char> tmp = getCharArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getCharArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a complex entry.  If the key does not exist, then false    *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isComplex(const std::string& key)
{
   bool is_complex  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set = H5Dopen(d_group_id, key.c_str());
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_COMPLEX_ARRAY) {
            is_complex = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_complex;
}

/*
*************************************************************************
*                                                                       *
* Create a complex scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putComplex(
   const std::string& key,
   const dcomplex& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putComplexArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a complex array entry in the database with the specified       *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putComplexArray(
   const std::string& key,
   const Array<dcomplex>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putComplexArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putComplexArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a complex array entry in the database with the specified       *
* key name.  The complex array is a compound type based on the hdf      *
* type H5T_NATIVE_DOUBLE (for real and imag parts).                     * 
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putComplexArray(
   const std::string& key,
   const dcomplex* const data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (dcomplex*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hid_t space, dataset;

      // Memory type
      hid_t mtype = createCompoundComplex('n');
      // Storage type
      hid_t stype = createCompoundComplex('s');

      hsize_t dim[] = {nelements};
      space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( space >= 0 );
#endif
   
      dataset = 
         H5Dcreate( d_group_id, key.c_str(), stype, space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_COMPLEX_ARRAY, dataset );

      errf = H5Tclose(mtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Tclose(stype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putComplexArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get complex scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a complex type.                                                 *
*                                                                      *
************************************************************************
*/

dcomplex HDFDatabase::getComplex(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   dcomplex ret_val;
   getComplexArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get complex scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a complex type.                                                 *
*                                                                      *
************************************************************************
*/

dcomplex HDFDatabase::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<dcomplex> local_dcomplex = getComplexArray(key);
   dcomplex *locptr = local_dcomplex.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get complex arrays from the database with the        *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a complex type.                   *
*                                                                      *
************************************************************************
*/

Array<dcomplex> HDFDatabase::getComplexArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   herr_t errf;
   if ( !isComplex(key) ) {
      TBOX_ERROR("HDFDatabase::getComplexArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a complex array." << std::endl);
   }

   hid_t   dset, dspace;
   hsize_t nsel;

   // Memory type
   hid_t mtype = createCompoundComplex('n');

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   Array<dcomplex> complexArray(nsel);

   if (nsel > 0) {
      dcomplex* locPtr = complexArray.getPointer();
      errf = H5Dread(dset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   errf = H5Tclose(mtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );

#endif
   return complexArray;
}

void HDFDatabase::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<dcomplex> tmp = getComplexArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getComplexArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}


hid_t HDFDatabase::createCompoundComplex( char type_spec ) const {
   herr_t errf;
   hid_t double_type_spec=H5T_SAMRAI_DOUBLE;
   switch (type_spec) {
   case 'n':
      // Use native type specs.
      double_type_spec = H5T_NATIVE_DOUBLE;
      break;
   case 's':
      // Use storage type specs.
      double_type_spec = H5T_SAMRAI_DOUBLE;
      break;
   default:
      TBOX_ERROR("HDFDatabase::createCompundComplex() error in database "
         << d_database_name
         << "\n    Unknown type specifier found. " << std::endl);
   }
   hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(dcomplex));
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( type >= 0 );
#endif
   errf = H5Tinsert(type, "real", HOFFSET(hdf_complex,re), double_type_spec);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Tinsert(type, "imag", HOFFSET(hdf_complex,im), double_type_spec);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   return type;
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a double entry.  If the key does not exist, then false     *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isDouble(const std::string& key)
{
   bool is_double = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_DOUBLE_ARRAY) {
            is_double = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_double;
}

/*
*************************************************************************
*                                                                       *
* Create a double scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDouble(
   const std::string& key, 
   const double& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putDoubleArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDoubleArray(
   const std::string& key,
   const Array<double>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putDoubleArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putDoubleArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.  The array type is based on the hdf type H5T_NATIVE_HDOUBLE.*
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDoubleArray(
   const std::string& key,
   const double* const data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (double*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( space >= 0 );

#endif
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_DOUBLE, 
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );

#endif
      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_DOUBLE_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );

#endif
   } else {
      TBOX_ERROR("HDFDatabase::putDoubleArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   } 
}

/*
************************************************************************
*                                                                      *
* Get double scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a double type.                                                  *
*                                                                      *
************************************************************************
*/

double HDFDatabase::getDouble(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   double ret_val;
   getDoubleArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get double scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a double type.                                                  *
*                                                                      *
************************************************************************
*/

double HDFDatabase::getDoubleWithDefault(
   const std::string& key, 
   const double& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<double> local_double = getDoubleArray(key);
   double *locptr = local_double.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get double arrays from the database with the         *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a double type.                    *
*                                                                      *
************************************************************************
*/

Array<double> HDFDatabase::getDoubleArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   herr_t errf;
   if (!isDouble(key)) {
     TBOX_ERROR("HDFDatabase::getDoubleArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a double array." << std::endl);
   }

   hid_t   dset, dspace;
   hsize_t nsel;

   dset = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   Array<double> doubleArray(nsel);

   if (nsel > 0) {
      double* locPtr = doubleArray.getPointer();
      errf = H5Dread(dset, H5T_NATIVE_DOUBLE, 
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif

   return doubleArray;
}

void HDFDatabase::getDoubleArray(
   const std::string& key, 
   double* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<double> tmp = getDoubleArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getDoubleArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a float entry.  If the key does not exist, then false      *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isFloat(const std::string& key)
{
   bool is_float  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_FLOAT_ARRAY) {
            is_float = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_float;
}

/*
*************************************************************************
*                                                                       *
* Create a float scalar entry in the database with the specified        *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putFloat(
   const std::string& key, 
   const float& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putFloatArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a float array entry in the database with the specified         *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putFloatArray(
   const std::string& key,
   const Array<float>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putFloatArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putFloatArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
  
}

/*
*************************************************************************
*                                                                       *
* Create a float array entry in the database with the specified         *
* key name.  The array type is based on the hdf type H5T_NATIVE_HFLOAT. *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putFloatArray(
   const std::string& key,
   const float* const data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (float*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( space >= 0 );

#endif
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_FLOAT, 
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_FLOAT_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putFloatArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get float scalar entry from the database with the specified key      *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a float type.                                                   *
*                                                                      *
************************************************************************
*/

float HDFDatabase::getFloat(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   float ret_val;
   getFloatArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get float scalar entry from the database with the specified key      *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a float type.                                                   *
*                                                                      *
************************************************************************
*/

float HDFDatabase::getFloatWithDefault(
   const std::string& key, 
   const float& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<float> local_float = getFloatArray(key);
   float *locptr = local_float.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get float arrays from the database with the          *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a float type.                     *
*                                                                      *
************************************************************************
*/

Array<float> HDFDatabase::getFloatArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   herr_t errf;
   if (!isFloat(key)) {
      TBOX_ERROR("HDFDatabase::getFloatArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a float array." << std::endl);
   }
 
   hid_t   dset, dspace;
   hsize_t nsel;

   dset = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   Array<float> floatArray(nsel);
 
   if (nsel > 0) {
      float* locPtr = floatArray.getPointer();
      errf = H5Dread(dset, H5T_NATIVE_FLOAT, 
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif

   return floatArray;

}

void HDFDatabase::getFloatArray(
   const std::string& key,
   float* data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<float> tmp = getFloatArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getFloatArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a integer entry.  If the key does not exist, then false    *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isInteger(const std::string& key)
{
   bool is_int  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_INT_ARRAY) {
            is_int = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_int;
}

/*
*************************************************************************
*                                                                       *
* Create a integer scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putInteger(
   const std::string& key,
   const int& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putIntegerArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create an integer array entry in the database with the specified      *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putIntegerArray(
   const std::string& key, 
   const Array<int>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putIntegerArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putIntegerArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create an integer array entry in the database with the specified      *
* key name.  The array type is based on the hdf type H5T_NATIVE_HINT.   *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putIntegerArray(
   const std::string& key, 
   const int* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (int*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT(space >= 0);
#endif

      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_INT, 
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT(errf >= 0);
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_INT_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putIntegerArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get integer scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a integer type.                                                 *
*                                                                      *
************************************************************************
*/

int HDFDatabase::getInteger(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   int ret_val;
   getIntegerArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get integer scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a integer type.                                                 *
*                                                                      *
************************************************************************
*/

int HDFDatabase::getIntegerWithDefault(
   const std::string& key, 
   const int& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<int> local_int = getIntegerArray(key);
   int *locptr = local_int.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get integer arrays from the database with the        *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a integer type.                   *
*                                                                      *
************************************************************************
*/

Array<int> HDFDatabase::getIntegerArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   herr_t errf;
   if ( !isInteger(key) ) {
      TBOX_ERROR("HDFDatabase::getIntegerArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not an integer array." << std::endl);
   }

   hid_t   dset, dspace;
   hsize_t nsel;

   dset = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   Array<int> intArray(nsel);

   if (nsel > 0) {
      int* locPtr = intArray.getPointer();
      errf = H5Dread(dset, H5T_NATIVE_INT, 
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );

#endif
   return intArray;
}

void HDFDatabase::getIntegerArray(
   const std::string& key, 
   int* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<int> tmp = getIntegerArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getIntegerArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a string entry.  If the key does not exist, then false     *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isString(const std::string& key)
{
   bool is_string  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_STRING_ARRAY) {
            is_string = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( errf >= 0 );
#endif
      }
   }

   return is_string;
}

/*
*************************************************************************
*                                                                       *
* Create a string scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putString(
   const std::string& key, 
   const std::string& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putStringArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a string array entry in the database with the specified        *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putStringArray(
   const std::string& key, 
   const Array<std::string>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putStringArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("HDFDatabase::putStringArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.  The array type is based on the hdf type H5T_C_S1.          *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putStringArray(
   const std::string& key, 
   const std::string* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != (std::string*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      int maxlen = 0;
      int current, data_size;
      int i;
      for (i = 0; i < nelements; i++) {
         current = data[i].size(); 
         if ( current > maxlen ) maxlen = current;
      }

      char* local_buf = new char[nelements*(maxlen+1)];
      for (i = 0; i < nelements; i++) {
         strcpy(&local_buf[i*(maxlen+1)], data[i].c_str());
         data_size = data[i].size();
         if (data_size < maxlen) {
            memset(&local_buf[i*(maxlen+1)] + data_size + 1, 0, 
                   maxlen - data_size);
         }
      }

      hid_t atype = H5Tcopy(H5T_C_S1);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( atype >= 0 );
#endif
      errf = H5Tset_size(atype, maxlen+1);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Tset_strpad(atype, H5T_STR_NULLTERM);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( space >= 0 );
#endif

      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), 
                                atype, space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( dataset >= 0 );
#endif

      errf = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, local_buf);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_STRING_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Tclose(atype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      delete[] local_buf;

   } else {
      TBOX_ERROR("HDFDatabase::putStringArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get string scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a string type.                                                  *
*                                                                      *
************************************************************************
*/

std::string HDFDatabase::getString(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   std::string ret_val;
   getStringArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Get string scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a string type.                                                  *
*                                                                      *
************************************************************************
*/

std::string HDFDatabase::getStringWithDefault(
   const std::string& key, 
   const std::string& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<std::string> local_string = getStringArray(key);
   std::string *locptr = local_string.getPointer();
   return ( locptr == NULL ? defaultvalue : *locptr);
}

/*
************************************************************************
*                                                                      *
* Two routines to get string arrays from the database with the         *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a string type.                    *
*                                                                      *
************************************************************************
*/

Array<std::string> HDFDatabase::getStringArray(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   herr_t errf;
   if (!isString(key)) {
      TBOX_ERROR("HDFDatabase::getStringArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a string array." << std::endl);
   }

   hsize_t nsel;
   size_t  dsize;
   hid_t   dset, dspace, dtype;
   char*   local_buf;

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dspace >= 0 );
#endif
   dtype  = H5Dget_type(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( dtype >= 0 );
#endif
   dsize  = H5Tget_size(dtype);
   nsel   = H5Sget_select_npoints(dspace);

   local_buf = new char[nsel*dsize];

   errf = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, local_buf);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif

   Array<std::string> stringArray(nsel);

   for (int i = 0; i < (int)nsel; i++) {
      std::string* locPtr = stringArray.getPointer(i);
      *locPtr = &local_buf[i*dsize];
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Tclose(dtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif

   delete[] local_buf;
   return stringArray;
}

void HDFDatabase::getStringArray(
   const std::string& key, 
   std::string* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   Array<std::string> tmp = getStringArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getStringArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}


void HDFDatabase::writeAttribute( int type_key,
                                       hid_t dataset_id
                                       )
{
   herr_t errf;
   hid_t attr_id = H5Screate(H5S_SCALAR);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( attr_id >= 0 );
#endif
   hid_t attr = H5Acreate(dataset_id, "Type", H5T_SAMRAI_ATTR, 
                          attr_id, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( attr >= 0 );
#endif
   errf = H5Awrite(attr, H5T_NATIVE_INT, &type_key);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Aclose(attr);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Sclose(attr_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
}


int HDFDatabase::readAttribute( hid_t dataset_id )
{
   herr_t errf;
   hid_t attr = H5Aopen_name(dataset_id, "Type");
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( attr >= 0 );
#endif
   int type_key;
   errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Aclose(attr);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   return type_key;
}


/*
*************************************************************************
*                                                                       *
* Print contents of current database to the specified output stream.    *
* Note that contents of subdatabases will not be printed.  This must    *
* be done by iterating through all the subdatabases individually.       * 
*                                                                       *
*************************************************************************
*/

void HDFDatabase::printClassData(std::ostream& os)
{

   performKeySearch();

   if (d_keydata.getNumberItems() == 0) {
      os << "Database named `"<< d_database_name 
         << "' has zero keys..." << std::endl;
   } else {
      os << "Printing contents of database named `" 
         << d_database_name << "'..." << std::endl;
   }

   for (List<KeyData>::Iterator i(d_keydata); i; i++) {
      int t = i().d_type; 
      switch ( tbox::MathUtilities<int>::Abs(t) ) {
         case KEY_DATABASE: {
            os << "   Data entry `"<< i().d_key << "' is"
               << " a database" << std::endl;   
            break;
         }
         case KEY_BOOL_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a boolean ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_BOX_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a box ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_CHAR_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a char ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_COMPLEX_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a complex ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_DOUBLE_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a double ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_FLOAT_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a float ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_INT_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " an integer ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         case KEY_STRING_ARRAY: {
            os << "   Data entry `"<< i().d_key << "' is" << " a string ";
            os << ( (t < 0) ? "scalar" : "array") << std::endl;   
            break;
         }
         default: {
            TBOX_ERROR("HDFDatabase::printClassData error....\n"
               << "   Unable to identify key = " << i().d_key
               << " as a known group or dataset" << std::endl);
         }
      }
   }

   cleanupKeySearch();

}

/*
*************************************************************************
*                                                                       *
* Open a HDF5 data-base file.  Flags can be "R" for read only or "W"    *
* for write.  If the file does not exist when read only is specified,   *
* an error status is returned (<0).  If the file exists when opened     *
* for write, an error status is returned (<0).                          * 
*                                                                       *
*************************************************************************
*/ 

int HDFDatabase::mount(
   const std::string& file_name, 
   const std::string& flags)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!file_name.empty());
   TBOX_ASSERT(!flags.empty());
#endif

   int status = 1;

   hid_t file_id = 0;

   if (flags == "R") {
     file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
     if( file_id < 0 ) {
        TBOX_ERROR("Unable to open HDF5 file " << file_name << "\n");
     }
   } else if (flags == "W") {
     file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, 
                         H5P_DEFAULT, H5P_DEFAULT);
     if( file_id < 0 ) {
        TBOX_ERROR("Unable to open HDF5 file " << file_name << "\n");
     }
   } else {
     TBOX_ERROR("HDFDatabase::mount error...\n"
                << "   database name is " << d_database_name
                << "\n    unrecognized flag = " << flags << std::endl);  
   }

   if (file_id < 0) {
     status = (int)file_id;
   } else {
     d_is_file  = true;
     d_group_id = file_id;
     d_file_id  = file_id;
   }

   return status;

}

/*
*************************************************************************
*                                                                       *
* Close the open HDF data file specified by d_file_id.                  *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::unmount()
{
   herr_t errf;
   if (d_is_file) {
      errf = H5Fclose(d_file_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      if ( d_group_id == d_file_id ) d_group_id = -1;
      d_file_id = -1;
      d_is_file = false;
   }
}

/*
*************************************************************************
*                                                                       *
* Private helper function for writing arrays in HDF5.  This function    *
* was deprecated in HDF5 1.4.  We replicate it here since it makes      *
* arrays easier to use in this database class.                          *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::insertArray(
   hid_t parent_id, 
   const char *name, 
   hsize_t offset, 
   int ndims, 
   const hsize_t dim[/*ndims*/], 
   const int *perm, 
   hid_t member_id) const
{
   herr_t errf;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 2))

   hid_t array = H5Tarray_create(member_id, ndims, dim, perm);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( array >= 0 );
#endif
   errf = H5Tinsert(parent_id, name, offset, array);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
   errf = H5Tclose(array);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
#else
   size_t newdim[H5S_MAX_RANK];
   for(int i = 0; i < ndims; i++) {
     newdim[i] = dim[i];
   }
    
   errf = H5Tinsert_array(parent_id, name, offset, ndims, newdim, perm, member_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
#endif
}

/*
*************************************************************************
*                                                                       *
* Private helper function for searching database keys.                  *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::performKeySearch()
{
   herr_t errf;
   if (d_is_file) {
      s_group_to_search = "/";
      s_top_level_search_group = "/";
      s_found_group = 1;
   } else {
      s_group_to_search = d_database_name;
      s_top_level_search_group = std::string();
      s_found_group = 0;
   }

   s_still_searching = 1;

   errf = H5Giterate(d_group_id, "/", NULL,
                     HDFDatabase::iterateKeys, (void*)this);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
}

void HDFDatabase::cleanupKeySearch()
{
   s_top_level_search_group = std::string();
   s_group_to_search = std::string();
   s_still_searching = 0;
   s_found_group = 0;

   d_keydata.clearItems();
}

/*
*************************************************************************
*                                                                       *
* Public method to return the group_id so VisIt can access an           *
* object's HDF database.                                                *
*                                                                       *
*************************************************************************
*/
hid_t HDFDatabase::getGroupId()
{
  return (d_group_id);
}

}
}

#endif
