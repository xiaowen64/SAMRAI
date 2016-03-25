//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/database/Database.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	An abstract base class for the SAMRAI database objects
//

#ifndef included_tbox_Database
#define included_tbox_Database

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_tbox_DatabaseBox
#include "tbox/DatabaseBox.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif


namespace SAMRAI {
   namespace tbox {


/**
 * @brief Class Database is an abstract base class for the input, restart,
 * and visualization databases.  
 *
 * SAMRAI databases store (key,value) pairs in a hierarchical
 * database.  Each value may be another database or a boolean, box,
 * character, double complex, double, float, integer, or string.
 * DatabaseBoxes are stored using the toolbox box structure.
 *
 * Data is entered into the database through methods of the general form
 * putTYPE(key, TYPE) or putTYPEArray(key, TYPE array), where TYPE is the
 * type of value created.  If the specified key already exists in the
 * database, then the existing key is silently deleted.
 *
 * Data is extracted from the database through methods of the general form
 * TYPE = getTYPE(key), where TYPE is the type of value to be returned
 * from the database.  There are two general lookup methods.  In the first,
 * a default value is provided (for scalars only).  If the specified key is
 * not found in the database, then the specified default is returned.  In
 * the second form, no default is provided, and the database exists with
 * an error message and program exits if the key is not found.  The array
 * version of getTYPE() works in a similar fashion.
 */

class Database : public DescribedClass
{
public:
   /**
    * The constructor for the database base class does nothing interesting.
    */
   Database();

   /**
    * The virtual destructor for the database base class does nothing
    * interesting.
    */
   virtual ~Database();

   /**
    * Return true if the specified key exists in the database and false
    * otherwise.
    *
    * @param key Key name to lookup.
    */
   virtual bool keyExists(const std::string& key) = 0;

   /**
    * Return all keys in the database.
    */
   virtual Array<std::string> getAllKeys() = 0;

   /**
    * Return the size of the array associated with the key.  If the key
    * does not exist, then zero is returned.
    *
    * @param key Key name in database.
    */
   virtual int getArraySize(const std::string& key) = 0;

   /**
    * Return whether the specified key represents a database entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isDatabase(const std::string& key) = 0;

   /**
    * Create a new database with the specified key name.  If the key already
    * exists in the database, then the old key record is deleted and the new
    * one is silently created in its place.
    *
    * @param key Key name in database.
    */
   virtual Pointer<Database> putDatabase(const std::string& key) = 0;

   /**
    * Get the database with the specified key name.  If the specified
    * key does not exist in the database or it is not a database, then
    * an error message is printed and the program exits.
    *
    * @param key Key name in database.
    */
   virtual Pointer<Database> getDatabase(const std::string& key) = 0;

   /**
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isBool(const std::string& key) = 0;

   /**
    * Create a boolean scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putBool(const std::string& key, const bool& data) = 0;

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putBoolArray(
      const std::string& key, const Array<bool>& data) = 0;

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putBoolArray(
      const std::string& key, const bool* const data, const int nelements) = 0;

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual bool getBool(const std::string& key) = 0;

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a boolean scalar,
    * then an error message is printed and the program exits.
    */
   virtual bool getBoolWithDefault(
      const std::string& key, const bool& defaultvalue) = 0;

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<bool> getBoolArray(const std::string& key) = 0;

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getBoolArray(
      const std::string& key, bool* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents a box entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isDatabaseBox(const std::string& key) = 0;

   /**
    * Create a box scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Data to put into database.
    */
   virtual void putDatabaseBox(const std::string& key, const DatabaseBox& data) = 0;

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putDatabaseBoxArray(
      const std::string& key, const Array<DatabaseBox>& data) = 0;

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putDatabaseBoxArray(
      const std::string& key, const DatabaseBox* const data, const int nelements) = 0;

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * box scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual DatabaseBox getDatabaseBox(const std::string& key) = 0;

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a box scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual DatabaseBox getDatabaseBoxWithDefault(
      const std::string& key, const DatabaseBox& defaultvalue) = 0;

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<DatabaseBox> getDatabaseBoxArray(const std::string& key) = 0;

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getDatabaseBoxArray(
      const std::string& key, DatabaseBox* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isChar(const std::string& key) = 0;

   /**
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putChar(const std::string& key, const char& data) = 0;

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putCharArray(
      const std::string& key, const Array<char>& data) = 0;

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putCharArray(
      const std::string& key, const char* const data, const int nelements) = 0;

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual char getChar(const std::string& key) = 0;

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a character scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual char getCharWithDefault(
      const std::string& key, const char& defaultvalue) = 0;

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<char> getCharArray(const std::string& key) = 0;

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getCharArray(
      const std::string& key, char* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents a complex entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isComplex(const std::string& key) = 0;

   /**
    * Create a complex scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putComplex(const std::string& key, const dcomplex& data) = 0;

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putComplexArray(
      const std::string& key, const Array<dcomplex>& data) = 0;

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putComplexArray(
      const std::string& key, const dcomplex* const data, const int nelements) = 0;

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * complex scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual dcomplex getComplex(const std::string& key) = 0;

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a complex scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual dcomplex getComplexWithDefault(
      const std::string& key, const dcomplex& defaultvalue) = 0;

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<dcomplex> getComplexArray(const std::string& key) = 0;

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getComplexArray(
      const std::string& key, dcomplex* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isDouble(const std::string& key) = 0;

   /**
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putDouble(const std::string& key, const double& data) = 0;

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putDoubleArray(
      const std::string& key, const Array<double>& data) = 0;

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putDoubleArray(
      const std::string& key, const double* const data, const int nelements) = 0;

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual double getDouble(const std::string& key) = 0;

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a double scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual double getDoubleWithDefault(
      const std::string& key, const double& defaultvalue) = 0;

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<double> getDoubleArray(const std::string& key) = 0;

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getDoubleArray(
      const std::string& key, double* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isFloat(const std::string& key) = 0;

   /**
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putFloat(const std::string& key, const float& data) = 0;

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putFloatArray(const std::string& key, 
                              const Array<float>& data) = 0;

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putFloatArray(
      const std::string& key, const float* const data, const int nelements) = 0;

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual float getFloat(const std::string& key) = 0;

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a float scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual float getFloatWithDefault(
      const std::string& key, const float& defaultvalue) = 0;

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<float> getFloatArray(const std::string& key) = 0;

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getFloatArray(
      const std::string& key, float* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents an integer entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isInteger(const std::string& key) = 0;

   /**
    * Create an integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putInteger(const std::string& key, const int& data) = 0;

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putIntegerArray(const std::string& key, 
                                const Array<int>& data) = 0;

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putIntegerArray(
      const std::string& key, const int* const data, const int nelements) = 0;

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual int getInteger(const std::string& key) = 0;

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not an integer scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual int getIntegerWithDefault(
      const std::string& key, const int& defaultvalue) = 0;

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<int> getIntegerArray(const std::string& key) = 0;

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getIntegerArray(
      const std::string& key, int* data, const int nelements) = 0;

   /**
    * Return whether the specified key represents a std::string entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isString(const std::string& key) = 0;

   /**
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putString(const std::string& key, const std::string& data) = 0;

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Array with data to put into database.
    */
   virtual void putStringArray(const std::string& key, 
                               const Array<std::string>& data) = 0;

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void putStringArray(
      const std::string& key, const std::string* const data, const int nelements) = 0;

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * string scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual std::string getString(const std::string& key) = 0;

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a string scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    */
   virtual std::string getStringWithDefault(
      const std::string& key, const std::string& defaultvalue) = 0;

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual Array<std::string> getStringArray(const std::string& key) = 0;

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void getStringArray(
      const std::string& key, std::string* data, const int nelements) = 0;

   /**
    * Get a bool entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * bool scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void getScalar(const std::string& key, bool& scalar);

   /**
    * Get a bool entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a bool array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const std::string& key, const bool scalar);


   /**
    * Get a bool entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a bool array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void getArray(const std::string& key, Array<bool>& array);

   /**
    * Create an bool array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const std::string& key, const Array<bool> array);

   /**
    * Get a char entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * char scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void getScalar(const std::string& key, char& scalar);

   /**
    * Get a char entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a char array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const std::string& key, const char scalar);


   /**
    * Get a char entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a char array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void getArray(const std::string& key, Array<char>& array);

   /**
    * Create an char array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const std::string& key, const Array<char> array);


   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * complex scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void getScalar(const std::string& key, dcomplex& scalar);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const std::string& key, const dcomplex scalar);


   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void getArray(const std::string& key, Array<dcomplex>& array);

   /**
    * Create an complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const std::string& key, const Array<dcomplex> array);


   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void getScalar(const std::string& key, float& scalar);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const std::string& key, const float scalar);


   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void getArray(const std::string& key, Array<float>& array);

   /**
    * Create an float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const std::string& key, const Array<float> array);


   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void getScalar(const std::string& key, double& scalar);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const std::string& key, const double scalar);


   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void getArray(const std::string& key, Array<double>& array);

   /**
    * Create an double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const std::string& key, const Array<double> array);

   /**
    * Get a integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void getScalar(const std::string& key, int& scalar);

   /**
    * Get a integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const std::string& key, const int scalar);


   /**
    * Get a integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void getArray(const std::string& key, Array<int>& array);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const std::string& key, const Array<int> array);


   /**
    * Print the current database to the specified output stream.  If
    * no output stream is specified, then data is written to stream pout.
    * 
    * @param os Output stream.
    */
   virtual void printClassData(std::ostream& os = pout) = 0;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Database.I"
#endif
#endif
