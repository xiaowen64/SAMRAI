//
// File:	InputDatabase.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 431 $
// Modified:	$Date: 2005-06-10 16:51:19 -0700 (Fri, 10 Jun 2005) $
// Description:	An input database structure that stores (key,value) pairs
//

#ifndef included_tbox_InputDatabase
#define included_tbox_InputDatabase

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class InputDatabase stores (key,value) pairs in a hierarchical
 * database.  Each value may be another database, boolean, box, character,
 * complex, double, float, integer, or string.  Note that boxes are stored
 * using the toolbox box class that can store boxes of any dimension in the
 * same data structure.
 *
 * See the Database class documentation for a description of the
 * generic database interface.
 *
 * Note that the input database will attempt to promote numerical types
 * where appropriate.  The promotion chain is int -> float -> double ->
 * complex.  For example, an integer key will be promoted to a complex
 * value if isComplex() or getComplex() is called.  Double values will also
 * be truncated to floats (with loss of information) if a float call is
 * made on a double value.
 *
 * It is assumed that all processors will access the database in the same
 * manner.  Thus, all error messages are output to pout instead of perr.
 */

class InputDatabase : public Database
{
public:
   /**
    * The input database constructor creates an empty database with the
    * specified name.
    */
   InputDatabase(const string& name);

   /**
    * The input database destructor deallocates the data in the database.
    */
   virtual ~InputDatabase();

   /**
    * Return string name of input database object.
    */
   virtual string getName() const;

   /**
    * Return true if the specified key exists in the database and false
    * otherwise.
    */
   virtual bool keyExists(const string& key);

   /**
    * Return all keys in the database.
    */
   virtual Array<string> getAllKeys();

   /**
    * Return the size of the array associated with the key.  If the key
    * does not exist, then zero is returned.
    */
   virtual int getArraySize(const string& key);

   /**
    * Return whether the specified key represents a database entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool isDatabase(const string& key);

   /**
    * Create a new database with the specified key name.  If the key already
    * exists in the database, then the old key record is deleted and the new
    * one is silently created in its place.
    */
   virtual Pointer<Database> putDatabase(const string& key);

   /**
    * Get the database with the specified key name.  If the specified
    * key does not exist in the database or it is not a database, then
    * an error message is printed and the program exits.
    */
   virtual Pointer<Database> getDatabase(const string& key);

   /**
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool isBool(const string& key);

   /**
    * Create a boolean scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putBool(const string& key, const bool& data);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putBoolArray(const string& key, 
                             const Array<bool>& data);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putBoolArray(
      const string& key, const bool* const data, const int nelements);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    */
   virtual bool getBool(const string& key);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a boolean scalar,
    * then an error message is printed and the program exits.
    */
   virtual bool getBoolWithDefault(const string& key, const bool& defaultvalue);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.
    */
   virtual Array<bool> getBoolArray(const string& key);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void getBoolArray(
      const string& key, bool* data, const int nelements);

   /**
    * Return whether the specified key represents a box entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool isDatabaseBox(const string& key);

   /**
    * Create a box scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putDatabaseBox(const string& key, const DatabaseBox& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putDatabaseBoxArray(const string& key, 
                            const Array<DatabaseBox>& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putDatabaseBoxArray(
      const string& key, const DatabaseBox* const data, const int nelements);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * box scalar, then an error message is printed and the program
    * exits.
    */
   virtual DatabaseBox getDatabaseBox(const string& key);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a box scalar,
    * then an error message is printed and the program exits.
    */
   virtual DatabaseBox getDatabaseBoxWithDefault(
      const string& key, const DatabaseBox& defaultvalue);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.
    */
   virtual Array<DatabaseBox> getDatabaseBoxArray(const string& key);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void getDatabaseBoxArray(
      const string& key, DatabaseBox* data, const int nelements);

   /**
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool isChar(const string& key);

   /**
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putChar(const string& key, const char& data);

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putCharArray(const string& key, 
                             const Array<char>& data);

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putCharArray(
      const string& key, const char* const data, const int nelements);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    */
   virtual char getChar(const string& key);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a character scalar,
    * then an error message is printed and the program exits.
    */
   virtual char getCharWithDefault(const string& key, const char& defaultvalue);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.
    */
   virtual Array<char> getCharArray(const string& key);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void getCharArray(
      const string& key, char* data, const int nelements);

   /**
    * Return whether the specified key represents a complex entry.  If
    * the key does not exist, then false is returned.  Complex values
    * may be promoted from integers, floats, or doubles.
    */
   virtual bool isComplex(const string& key);

   /**
    * Create a complex scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putComplex(const string& key, const dcomplex& data);

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putComplexArray(const string& key, 
                                const Array<dcomplex>& data);

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putComplexArray(
      const string& key, const dcomplex* const data, const int nelements);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * complex scalar, then an error message is printed and the program
    * exits.  Complex values may be promoted from integers, floats, or
    * doubles.
    */
   virtual dcomplex getComplex(const string& key);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a complex scalar,
    * then an error message is printed and the program exits.  Complex
    * values may be promoted from integers, floats, or doubles.
    */
   virtual dcomplex getComplexWithDefault(
      const string& key, const dcomplex& defaultvalue);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  Complex values may be promoted from integers,
    * floats, or doubles.
    */
   virtual Array<dcomplex> getComplexArray(const string& key);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Complex values may be promoted from integers, floats, or doubles.
    */
   virtual void getComplexArray(
      const string& key, dcomplex* data, const int nelements);

   /**
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.  Double values
    * may be promoted from integers or floats.
    */
   virtual bool isDouble(const string& key);

   /**
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putDouble(const string& key, const double& data);

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putDoubleArray(const string& key, 
                               const Array<double>& data);

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putDoubleArray(
      const string& key, const double* const data, const int nelements);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.  Double values may be promoted from integers or floats.
    */
   virtual double getDouble(const string& key);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a double scalar, then
    * an error message is printed and the program exits.  Double values may
    * be promoted from integers or floats.
    */
   virtual double getDoubleWithDefault(
      const string& key, const double& defaultvalue);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  Double values may be promoted from integers
    * or floats.
    */
   virtual Array<double> getDoubleArray(const string& key);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Double values may be promoted from integers or floats.
    */
   virtual void getDoubleArray(
      const string& key, double* data, const int nelements);

   /**
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.  Float values
    * may be promoted from integers or silently truncated from doubles.
    */
   virtual bool isFloat(const string& key);

   /**
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putFloat(const string& key, const float& data);

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putFloatArray(const string& key, 
                              const Array<float>& data);

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putFloatArray(
      const string& key, const float* const data, const int nelements);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.  Float values may be promoted from integers or silently
    * truncated from doubles.
    */
   virtual float getFloat(const string& key);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a float scalar, then
    * an error message is printed and the program exits.  Float values may
    * be promoted from integers or silently truncated from doubles.
    */
   virtual float getFloatWithDefault(
      const string& key, const float& defaultvalue);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  Float values may be promoted from integers
    * or silently truncated from doubles.
    */
   virtual Array<float> getFloatArray(const string& key);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Float values may be promoted from integers or silently truncated
    * from doubles.
    */
   virtual void getFloatArray(
      const string& key, float* data, const int nelements);

   /**
    * Return whether the specified key represents an integer entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool isInteger(const string& key);

   /**
    * Create an integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putInteger(const string& key, const int& data);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putIntegerArray(const string& key, 
                                const Array<int>& data);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putIntegerArray(
      const string& key, const int* const data, const int nelements);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    */
   virtual int getInteger(const string& key);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not an integer scalar,
    * then an error message is printed and the program exits.
    */
   virtual int getIntegerWithDefault(
      const string& key, const int& defaultvalue);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.
    */
   virtual Array<int> getIntegerArray(const string& key);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void getIntegerArray(
      const string& key, int* data, const int nelements);

   /**
    * Return whether the specified key represents a string entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool isString(const string& key);

   /**
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putString(const string& key, const string& data);

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putStringArray(const string& key, 
                               const Array<string>& data);

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void putStringArray(
      const string& key, const string* const data, const int nelements);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * string scalar, then an error message is printed and the program
    * exits.
    */
   virtual string getString(const string& key);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a string scalar,
    * then an error message is printed and the program exits.
    */
   virtual string getStringWithDefault(
      const string& key, const string& defaultvalue);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.
    */
   virtual Array<string> getStringArray(const string& key);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void getStringArray(
      const string& key, string* data, const int nelements);

   /**
    * Return whether the specified key has been accessed by one of the
    * lookup member functions.  If the key does not exist in the database,
    * then false is returned.
    */
   bool keyAccessed(const string& key);

   /**
    * Print the current database to the specified output stream.  After
    * each key, print whether that key came from the input file and was
    * used (input), came from the input file but was not used (unused),
    * or came from a default key value (default).  If no output stream
    * is specified, then data is written to stream pout.
    *
    * NOTE:  under the g++ compiler libraries, printClassData has a 
    * maximum output of 4096 characters per line.
    */
   virtual void printClassData(ostream& os = pout);

   /**
    * Print the database keys that were not used to the specified output
    * stream.
    */
   void printUnusedKeys(ostream& os = pout) const;

   /**
    * Print the database keys that were set via default calls to the specified
    * output stream.
    */
   void printDefaultKeys(ostream& os = pout) const;

private:
   InputDatabase(const InputDatabase&);	// not implemented
   void operator=(const InputDatabase&);		// not implemented

   /*
    * The following structure holds the list of (key,value) pairs stored
    * in the database.  Note that only one of the arrays contains valid
    * data for any particular key.
    */
   struct KeyData {
      string                      d_key;		// key name
      int                         d_type;		// type of entry
      int                         d_array_size;		// size of array data
      bool                        d_accessed;		// whether accessed
      bool                        d_from_default;	// from default key
      Pointer<Database> d_database;		// sub-database
      Array<bool>            d_boolean;		// boolean array value
      Array<DatabaseBox>        d_box;		// box array value
      Array<char>            d_char;		// char array value
      Array<dcomplex>        d_complex;		// complex array value
      Array<double>          d_double;		// double array value
      Array<float>           d_float;		// float array value
      Array<int>             d_integer;		// integer array value
      Array<string>          d_string;		// string array value
   };

   /*
    * Private utility routines for managing the database
    */
   bool deleteKeyIfFound(const string& key);
   KeyData* findKeyData(const string& key);
   KeyData* findKeyDataOrExit(const string& key);
   static void indentStream(ostream& os, const int indent);
   void printDatabase(ostream& os, const int indent, const int toprint) const;

   /*
    * Private data members - name and a list of (key,value) pairs
    */
   string d_database_name;
   List<KeyData> d_keyvalues;
};


}
}


#ifndef DEBUG_NO_INLINE
#include "tbox/InputDatabase.I"
#endif
#endif
