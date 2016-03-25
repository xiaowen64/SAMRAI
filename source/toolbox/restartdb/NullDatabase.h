//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/restartdb/NullDatabase.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	A null database that does nothing for all database methods.
//

#ifndef included_tbox_NullDatabase
#define included_tbox_NullDatabase

#include "SAMRAI_config.h"
#include "tbox/Database.h"
#include "tbox/List.h"
#include "tbox/Utilities.h"


namespace SAMRAI {
   namespace tbox {


/**
 * The NullDatabase provides an implementation of the Database 
 * interface with empty methods for the purpose of reducing the
 * the number of guards necessary in methods from other classes that 
 * use databases.
 *
 * See the Database class documentation for a description of the
 * generic database interface.
 *
 */

class NullDatabase : public Database
{
public:
   /**
    * The null database constructor creates an empty database with 
    * the name "null".
    */
   NullDatabase();

   /**
    * The input database destructor deallocates the data in the database.
    */
   virtual ~NullDatabase();

      /**
    * Create a new database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    */
   virtual bool create(const std::string& name);


   /**
    * Open an existing database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    */
   virtual bool open(const std::string& name);


   /**
    * Close the database.
    *
    * Returns true if successful.
    *
    * If the database is currently open then close it.  This should
    * flush all data to the file (if the database is on disk).
    */
   virtual bool close();

   /**
    * Always returns true.
    */
   virtual bool keyExists(const std::string& key);

   /**
    * Return an empty Array<string>.
    */
   virtual Array<std::string> getAllKeys();

   /**
    * Return INVALID.
    */
   virtual enum DataType getArrayType(const std::string& key);

   /**
    * Always returns 0.
    */
   virtual int getArraySize(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isDatabase(const std::string& key);

   /**
    * Returns a pointer to the null database.
    */
   virtual Pointer<Database> putDatabase(const std::string& key);

   /**
    * Returns a pointer to the null database.
    */
   virtual Pointer<Database> getDatabase(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isBool(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putBoolArray(
      const std::string& key, const bool* const data, const int nelements);

   /**
    * Returns an empty Array<bool>.
    */
   virtual Array<bool> getBoolArray(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isDatabaseBox(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBoxArray(
      const std::string& key, const DatabaseBox* const data, const int nelements);

   /**
    * Returns an empty Array<box>.
    */
   virtual Array<DatabaseBox> getDatabaseBoxArray(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isChar(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putCharArray(
      const std::string& key, const char* const data, const int nelements);

   /**
    * Returns an empty Array<char>.
    */
   virtual Array<char> getCharArray(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isComplex(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putComplexArray(
      const std::string& key, const dcomplex* const data, const int nelements);

   /**
    * Returns an empty Array<dcomplex>.
    */
   virtual Array<dcomplex> getComplexArray(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isDouble(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putDoubleArray(
      const std::string& key, const double* const data, const int nelements);

   /**
    * Returns an empty Array<double>.
    */
   virtual Array<double> getDoubleArray(const std::string& key);

   /**
    * Always return true.
    */
   virtual bool isFloat(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putFloatArray(
      const std::string& key, const float* const data, const int nelements);

   /**
    * Returns an empty Array<float>.
    */
   virtual Array<float> getFloatArray(const std::string& key);

   /**
    * Always returns true. 
    */
   virtual bool isInteger(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putIntegerArray(
      const std::string& key, const int* const data, const int nelements);

   /**
    * Returns an empty Array<int>.
    */
   virtual Array<int> getIntegerArray(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool isString(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putStringArray(
      const std::string& key, const std::string* const data, const int nelements);

   /**
    * Returns an empty Array<std::string>.
    */
   virtual Array<std::string> getStringArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual std::string getName();

   /**
    * Does nothing.
    */
   virtual void printClassData(std::ostream& os = pout);

   using Database::putBoolArray;
   using Database::getBoolArray;
   using Database::putDatabaseBoxArray;
   using Database::getDatabaseBoxArray;
   using Database::putCharArray;
   using Database::getCharArray;
   using Database::putComplexArray;
   using Database::getComplexArray;
   using Database::putFloatArray;
   using Database::getFloatArray;
   using Database::putDoubleArray;
   using Database::getDoubleArray;
   using Database::putIntegerArray;
   using Database::getIntegerArray;
   using Database::putStringArray;
   using Database::getStringArray;

private:
   NullDatabase(const NullDatabase&);	// not implemented
   void operator=(const NullDatabase&);		// not implemented

};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/NullDatabase.I"
#endif
#endif
