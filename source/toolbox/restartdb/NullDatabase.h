//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/restartdb/NullDatabase.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	A null database that does nothing for all database methods.
//

#ifndef included_tbox_NullDatabase
#define included_tbox_NullDatabase

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif


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
    * Always returns true.
    */
   virtual bool keyExists(const std::string& key);

   /**
    * Return an empty Array<string>.
    */
   virtual Array<std::string> getAllKeys();

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
   virtual void putBool(const std::string& key, const bool& data);

   /**
    * Does nothing.
    */
   virtual void putBoolArray(const std::string& key, const Array<bool>& data);

   /**
    * Does nothing.
    */
   virtual void putBoolArray(
      const std::string& key, const bool* const data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool getBool(const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool getBoolWithDefault(
      const std::string& key, const bool& defaultvalue); 

   /**
    * Returns an empty Array<bool>.
    */
   virtual Array<bool> getBoolArray(const std::string& key);

   /**
    * Does nothing.  
    */
   virtual void getBoolArray(
      const std::string& key, bool* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isDatabaseBox(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBox(const std::string& key, const DatabaseBox& data);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBoxArray(
      const std::string& key, const Array<DatabaseBox>& data);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBoxArray(
      const std::string& key, const DatabaseBox* const data, const int nelements);

   /**
    * Returns a zero dimension empty box.
    */
   virtual DatabaseBox getDatabaseBox(const std::string& key);

   /**
    * Returns a zero dimension empty box.
    */
   virtual DatabaseBox getDatabaseBoxWithDefault(
      const std::string& key, const DatabaseBox& defaultvalue);

   /**
    * Returns an empty Array<box>.
    */
   virtual Array<DatabaseBox> getDatabaseBoxArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void getDatabaseBoxArray(
      const std::string& key, DatabaseBox* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isChar(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putChar(const std::string& key, const char& data);

   /**
    * Does nothing.
    */
   virtual void putCharArray(
      const std::string& key, const Array<char>& data);

   /**
    * Does nothing.
    */
   virtual void putCharArray(
      const std::string& key, const char* const data, const int nelements);

   /**
    * Always returns 0.
    */
   virtual char getChar(const std::string& key);

   /**
    * Always returns 0.
    */
   virtual char getCharWithDefault(const std::string& key, const char& defaultvalue);

   /**
    * Returns an empty Array<char>.
    */
   virtual Array<char> getCharArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void getCharArray(
      const std::string& key, char* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isComplex(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putComplex(const std::string& key, const dcomplex& data);

   /**
    * Does nothing.
    */
   virtual void putComplexArray(
      const std::string& key, const Array<dcomplex>& data);

   /**
    * Does nothing.
    */
   virtual void putComplexArray(
      const std::string& key, const dcomplex* const data, const int nelements);

   /**
    * Returns a 0.0 + 0.0i
    */
   virtual dcomplex getComplex(const std::string& key);

   /**
    * Returns a 0.0 + 0.0i
    */
   virtual dcomplex getComplexWithDefault(
      const std::string& key, const dcomplex& defaultvalue);


   /**
    * Returns an empty Array<dcomplex>.
    */
   virtual Array<dcomplex> getComplexArray(const std::string& key);

   /** 
    * Does nothing.
    */
   virtual void getComplexArray(
      const std::string& key, dcomplex* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isDouble(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putDouble(const std::string& key, const double& data);

   /**
    * Does nothing.
    */
   virtual void putDoubleArray(
      const std::string& key, const Array<double>& data);

   /**
    * Does nothing.
    */
   virtual void putDoubleArray(
      const std::string& key, const double* const data, const int nelements);

   /**
    * Returns 0.0
    */
   virtual double getDouble(const std::string& key);

   /**
    * Returns 0.0
    */
   virtual double getDoubleWithDefault(
      const std::string& key, const double& defaultvalue);

   /**
    * Returns an empty Array<double>.
    */
   virtual Array<double> getDoubleArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void getDoubleArray(
      const std::string& key, double* data, const int nelements);

   /**
    * Always return true.
    */
   virtual bool isFloat(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putFloat(const std::string& key, const float& data);

   /**
    * Does nothing.
    */
   virtual void putFloatArray(
      const std::string& key, const Array<float>& data);

   /**
    * Does nothing.
    */
   virtual void putFloatArray(
      const std::string& key, const float* const data, const int nelements);

   /**
    * Returns 0.0
    */
   virtual float getFloat(const std::string& key);

   /**
    * Returns 0.0
    */
   virtual float getFloatWithDefault(
      const std::string& key, const float& defaultvalue);

   /**
    * Returns an empty Array<float>.
    */
   virtual Array<float> getFloatArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void getFloatArray(
      const std::string& key, float* data, const int nelements);

   /**
    * Always returns true. 
    */
   virtual bool isInteger(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putInteger(const std::string& key, const int& data);

   /**
    * Does nothing.
    */
   virtual void putIntegerArray(
      const std::string& key, const Array<int>& data);

   /**
    * Does nothing.
    */
   virtual void putIntegerArray(
      const std::string& key, const int* const data, const int nelements);

   /**
    * Returns 0.
    */
   virtual int getInteger(const std::string& key);

   /**
    * Returns 0.
    */
   virtual int getIntegerWithDefault(
      const std::string& key, const int& defaultvalue);

   /**
    * Returns an empty Array<int>.
    */
   virtual Array<int> getIntegerArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void getIntegerArray(
      const std::string& key, int* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isString(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void putString(const std::string& key, const std::string& data);

   /**
    * Does nothing.
    */
   virtual void putStringArray(
      const std::string& key, const Array<std::string>& data);

   /**
    * Does nothing.
    */
   virtual void putStringArray(
      const std::string& key, const std::string* const data, const int nelements);

   /**
    * Returns and empty string.
    */
   virtual std::string getString(const std::string& key);

   /**
    * Returns and empty string.
    */
   virtual std::string getStringWithDefault(
      const std::string& key, const std::string& defaultvalue);

   /**
    * Returns an empty Array<std::string>.
    */
   virtual Array<std::string> getStringArray(const std::string& key);

   /**
    * Does nothing.
    */
   virtual void getStringArray(
      const std::string& key, std::string* data, const int nelements);

   /**
    * Does nothing.
    */
   virtual void printClassData(std::ostream& os = pout);

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
