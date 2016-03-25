//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/restartdb/HDFDatabase.C $
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2620 $
// Modified:    $LastChangedDate: 2008-11-19 14:24:28 -0800 (Wed, 19 Nov 2008) $
// Description: A database structure that stores HDF5 format data.
//

#include <cstring>

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

/*
 * SGS Note:  Can the new HDF5 stack stuff be a better solution to this?
 */
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
#define BEGIN_SUPPRESS_HDF5_WARNINGS                  \
{                                                     \
   herr_t (*H5E_saved_efunc) (hid_t, void*) = NULL;   \
   void *H5E_saved_edata = NULL;                      \
   H5Eget_auto(H5E_DEFAULT, &H5E_saved_efunc, &H5E_saved_edata); \
   H5Eset_auto(H5E_DEFAULT, NULL, NULL);                              

#define END_SUPPRESS_HDF5_WARNINGS                     \
   H5Eset_auto(H5E_DEFAULT, H5E_saved_efunc, H5E_saved_edata);	\
}
#else
#define BEGIN_SUPPRESS_HDF5_WARNINGS                  \
{                                                     \
   herr_t (*H5E_saved_efunc) (void*) = NULL;          \
   void *H5E_saved_edata = NULL;                      \
   H5Eget_auto(&H5E_saved_efunc, &H5E_saved_edata);   \
   H5Eset_auto(NULL, NULL);                              

#define END_SUPPRESS_HDF5_WARNINGS                     \
   H5Eset_auto(H5E_saved_efunc, H5E_saved_edata);      \
}
#endif



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
   void *void_database)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(name != (char*)NULL);
#endif

   HDFDatabase *database = (HDFDatabase *)(void_database);

   if (database -> d_still_searching) {

     H5G_stat_t statbuf;
     int type_key;
     herr_t errf;

     errf = H5Gget_objinfo(loc_id, name, 0, &statbuf);
#ifdef ASSERT_HDF5_RETURN_VALUES
     TBOX_ASSERT( errf >= 0 );
#endif

     switch (statbuf.type) {
     case H5G_GROUP: {
       if (database -> d_top_level_search_group == "/") {
         addKeyToList(name, KEY_DATABASE, void_database);
       } else if ( !strcmp(name, database -> d_group_to_search.c_str()) ) {
         hid_t grp;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
	 grp = H5Gopen(loc_id, name, H5P_DEFAULT);
#else	 
         grp = H5Gopen(loc_id, name);
#endif

#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( grp >= 0 );
#endif
         database -> d_found_group = true;
         database -> d_still_searching =
           H5Giterate(grp, ".", NULL,
                      HDFDatabase::iterateKeys, void_database);
#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( database -> d_still_searching >= 0 );
#endif
         database -> d_found_group = false;
       } else {
         hid_t grp;

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
	 grp = H5Gopen(loc_id, name, H5P_DEFAULT);
#else	 
         grp = H5Gopen(loc_id, name);
#endif

#ifdef ASSERT_HDF5_RETURN_VALUES
         TBOX_ASSERT( grp >= 0 );
#endif
         if (database -> d_found_group) {
           addKeyToList(name, KEY_DATABASE, void_database);
         } else {
           errf = H5Giterate(grp, ".", NULL,
                             HDFDatabase::iterateKeys, void_database);
#ifdef ASSERT_HDF5_RETURN_VALUES
           TBOX_ASSERT( errf >= 0 );
#endif
         }
       }
       break;
     }

     case H5G_DATASET: {
       if (database -> d_still_searching && database -> d_found_group) {
         hid_t this_set;
         BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(loc_id, name, H5P_DEFAULT);
#else	 
      this_set = H5Dopen(loc_id, name);
#endif
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
                         void_database);
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
   d_still_searching(0),
   d_found_group(0),
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
      close();
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
   BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	
   this_set = H5Dopen(d_group_id, key.c_str()); 
#endif
   END_SUPPRESS_HDF5_WARNINGS;
   if (this_set > 0) {
      key_exists = true;
      errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
   }
   if (!key_exists) {
      hid_t this_group;
      BEGIN_SUPPRESS_HDF5_WARNINGS;

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_group = H5Gopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_group = H5Gopen(d_group_id, key.c_str());
#endif

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

   Array<std::string> tmp_keys(d_keydata.getNumberOfItems());

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
*									*
* Get the type of the array entry associated with the specified key	*
*									*
*************************************************************************
*/
enum Database::DataType HDFDatabase::getArrayType(const std::string& key)
{

   enum Database::DataType type = Database::SAMRAI_INVALID;

   herr_t errf;

   if (!key.empty()) {

      if(isDatabase(key)) {
	 type = Database::SAMRAI_DATABASE;
      } else {

	 hid_t this_set;
	 BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
	 this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
	 this_set = H5Dopen(d_group_id, key.c_str());
#endif
	 END_SUPPRESS_HDF5_WARNINGS;
	    if (this_set > 0) {
	       int type_key = readAttribute(this_set);
	       
	       switch (type_key) {
		  case KEY_DATABASE:
		     type = Database::SAMRAI_DATABASE;
		     break;
		  case KEY_BOOL_ARRAY:
		     type = Database::SAMRAI_BOOL;
		     break;
		  case KEY_BOX_ARRAY:
		     type = Database::SAMRAI_BOX;
		     break;
		  case KEY_CHAR_ARRAY:
		     type = Database::SAMRAI_CHAR;
		     break;
		  case KEY_COMPLEX_ARRAY:
		     type = Database::SAMRAI_COMPLEX;
		     break;
		  case KEY_DOUBLE_ARRAY:
		     type = Database::SAMRAI_DOUBLE;
		     break;
		  case KEY_FLOAT_ARRAY:
		     type = Database::SAMRAI_FLOAT;
		     break;
		  case KEY_INT_ARRAY:
		     type = Database::SAMRAI_INT;
		     break;
		  case KEY_STRING_ARRAY:
		     type = Database::SAMRAI_STRING;
		     break;
		  case KEY_BOOL_SCALAR:
		     type = Database::SAMRAI_BOOL;
		     break;
		  case KEY_BOX_SCALAR:
		     type = Database::SAMRAI_BOX;
		     break;
		  case KEY_CHAR_SCALAR:
		     type = Database::SAMRAI_CHAR;
		     break;
		  case KEY_COMPLEX_SCALAR:
		     type = Database::SAMRAI_COMPLEX;
		     break;
		  case KEY_DOUBLE_SCALAR:
		     type = Database::SAMRAI_DOUBLE;
		     break;
		  case KEY_FLOAT_SCALAR:
		     type = Database::SAMRAI_FLOAT;
		     break;
		  case KEY_INT_SCALAR:
		     type = Database::SAMRAI_INT;
		     break;
		  case KEY_STRING_SCALAR:
		     type = Database::SAMRAI_STRING;
		     break;
	       }

	       errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
	       TBOX_ASSERT( errf >= 0 );
#endif
	    }
      }
   }
   return type;
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
   BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   this_set = H5Dopen(d_group_id, key.c_str());
#endif
   END_SUPPRESS_HDF5_WARNINGS;
   if (this_set > 0) {
      hid_t this_space = H5Dget_space(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( this_space >= 0 );
#endif

      hsize_t nsel;
      if(readAttribute(this_set) == KEY_CHAR_ARRAY) {
	 hid_t dtype  = H5Dget_type(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
	 TBOX_ASSERT( dtype >= 0 );
#endif
	 nsel   = H5Tget_size(dtype);
      } else {
	 nsel = H5Sget_select_npoints(this_space);
      }
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
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif

   bool is_database = false;
   herr_t errf;

   hid_t this_group;
   BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   this_group = H5Gopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   this_group = H5Gopen(d_group_id, key.c_str());
#endif
   END_SUPPRESS_HDF5_WARNINGS;
   if (this_group > 0) {
      is_database = true;
      errf = H5Gclose(this_group);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   hid_t this_group = H5Gcreate(d_group_id, key.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT);
#else	
   hid_t this_group = H5Gcreate(d_group_id, key.c_str(), 0);
#endif


#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( this_group >= 0 );
#endif

   Pointer<Database> new_database = new HDFDatabase(key, this_group);

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   hid_t this_group = H5Gopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	
   hid_t this_group = H5Gopen(d_group_id, key.c_str()); 
#endif
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
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_BOOL,
                                space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_BOOL,
                                space, H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset   = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset   = H5Dopen(d_group_id, key.c_str());
#endif

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
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      hid_t dataset =
         H5Dcreate( d_group_id, key.c_str(), stype, space, H5P_DEFAULT, 
		    H5P_DEFAULT, H5P_DEFAULT);
#else	
      hid_t dataset =
         H5Dcreate( d_group_id, key.c_str(), stype, space, H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset   = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset   = H5Dopen(d_group_id, key.c_str());
#endif

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
      BEGIN_SUPPRESS_HDF5_WARNINGS;

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      dataset = H5Dcreate(d_group_id, key.c_str(), atype, space,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      dataset = H5Dcreate(d_group_id, key.c_str(), atype, space, 
			  H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset   = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset   = H5Dopen(d_group_id, key.c_str());
#endif
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
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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
   
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      dataset = H5Dcreate( d_group_id, key.c_str(), stype, space, 
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      dataset = H5Dcreate( d_group_id, key.c_str(), stype, space, 
			   H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset   = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset   = H5Dopen(d_group_id, key.c_str());
#endif

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
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_DOUBLE, 
                                space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_DOUBLE, 
                                space, H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset = H5Dopen(d_group_id, key.c_str());
#endif
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
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_FLOAT, 
                                space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_FLOAT, 
                                space, H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset = H5Dopen(d_group_id, key.c_str());
#endif

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
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_INT, 
                                space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_INT, 
                                space, H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	
   dset = H5Dopen(d_group_id, key.c_str()); 
#endif

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
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
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
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), 
                                atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else	
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), 
                                atype, space, H5P_DEFAULT);
#endif

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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   dset   = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else	 
   dset   = H5Dopen(d_group_id, key.c_str());
#endif

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

void HDFDatabase::writeAttribute( int type_key,
                                       hid_t dataset_id
                                       )
{
   herr_t errf;
   hid_t attr_id = H5Screate(H5S_SCALAR);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( attr_id >= 0 );
#endif
#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   hid_t attr = H5Acreate(dataset_id, "Type", H5T_SAMRAI_ATTR, 
                          attr_id, H5P_DEFAULT, H5P_DEFAULT);
#else	
   hid_t attr = H5Acreate(dataset_id, "Type", H5T_SAMRAI_ATTR, 
                          attr_id, H5P_DEFAULT);
#endif

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

   if (d_keydata.getNumberOfItems() == 0) {
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
* Create HDF data file specified by name.                               *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::create(const std::string& name) {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
#endif
   bool status = false;

   hid_t file_id = 0;

   file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, 
		       H5P_DEFAULT, H5P_DEFAULT);
   if( file_id < 0 ) {
      TBOX_ERROR("Unable to open HDF5 file " << name << "\n");
      status = false;
   } else {
      status = true;
      d_is_file  = true;
      d_group_id = file_id;
      d_file_id  = file_id;
   }

   return status;
}

/*
*************************************************************************
*                                                                       *
* Open HDF data file specified by name                                  *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::open(const std::string& name) {
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!name.empty());
#endif
   bool status = false;

   hid_t file_id = 0;

   file_id = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   if (file_id < 0) {
      TBOX_ERROR("Unable to open HDF5 file " << name << "\n");
      status = false;
   } else {
      status = true;
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

bool HDFDatabase::close()
{
   herr_t errf =0;
   if (d_is_file) {
      errf = H5Fclose(d_file_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
      TBOX_ASSERT( errf >= 0 );
#endif
      if ( d_group_id == d_file_id ) {
	 d_group_id = -1;
      }
      d_file_id = -1;
      d_is_file = false;
   }

   if(errf >= 0) {
      return true;
   }
   else {
      return false;
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

#if (H5_VERS_MAJOR>1) || ((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR > 6))
   hid_t array = H5Tarray_create(member_id, ndims, dim);
#else	
   /*
    * Note that perm is NOT used by HDF, see HDF documentation.
    */
   hid_t array = H5Tarray_create(member_id, ndims, dim, perm);
#endif


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
      d_group_to_search = "/";
      d_top_level_search_group = "/";
      d_found_group = 1;
   } else {
      d_group_to_search = d_database_name;
      d_top_level_search_group = std::string();
      d_found_group = 0;
   }

   d_still_searching = 1;

   errf = H5Giterate(d_group_id, "/", NULL,
                     HDFDatabase::iterateKeys, (void*)this);
#ifdef ASSERT_HDF5_RETURN_VALUES
   TBOX_ASSERT( errf >= 0 );
#endif
}

void HDFDatabase::cleanupKeySearch()
{
   d_top_level_search_group = std::string();
   d_group_to_search = std::string();
   d_still_searching = 0;
   d_found_group = 0;

   d_keydata.clearItems();
}

/*
*************************************************************************
* Attach to an already created HDF file.
*************************************************************************
*/
bool HDFDatabase::attachToFile(hid_t group_id)
{
   bool status = false;

   if(group_id > 0) {
      status = true;
      d_is_file=false;
      d_file_id = -1;
      d_group_id = group_id; 
   } else {      
      TBOX_ERROR("HDFDatabase: Invalid fileid supplied to attachToFile" 
		 << std::endl);
      status = false;
   }

   return status;
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

std::string HDFDatabase::getName(void)
{
   return d_database_name;
}

}
}

#endif
