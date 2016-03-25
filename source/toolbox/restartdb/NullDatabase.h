//
// File:	NullDatabase.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
   virtual bool keyExists(const string& key);

   /**
    * Return an empty Array<string>.
    */
   virtual Array<string> getAllKeys();

   /**
    * Always returns 0.
    */
   virtual int getArraySize(const string& key);

   /**
    * Always returns true.
    */
   virtual bool isDatabase(const string& key);

   /**
    * Returns a pointer to the null database.
    */
   virtual Pointer<Database> putDatabase(const string& key);

   /**
    * Returns a pointer to the null database.
    */
   virtual Pointer<Database> getDatabase(const string& key);

   /**
    * Always returns true.
    */
   virtual bool isBool(const string& key);

   /**
    * Does nothing.
    */
   virtual void putBool(const string& key, const bool& data);

   /**
    * Does nothing.
    */
   virtual void putBoolArray(const string& key, const Array<bool>& data);

   /**
    * Does nothing.
    */
   virtual void putBoolArray(
      const string& key, const bool* const data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool getBool(const string& key);

   /**
    * Always returns true.
    */
   virtual bool getBoolWithDefault(
      const string& key, const bool& defaultvalue); 

   /**
    * Returns an empty Array<bool>.
    */
   virtual Array<bool> getBoolArray(const string& key);

   /**
    * Does nothing.  
    */
   virtual void getBoolArray(
      const string& key, bool* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isDatabaseBox(const string& key);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBox(const string& key, const DatabaseBox& data);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBoxArray(
      const string& key, const Array<DatabaseBox>& data);

   /**
    * Does nothing.
    */
   virtual void putDatabaseBoxArray(
      const string& key, const DatabaseBox* const data, const int nelements);

   /**
    * Returns a zero dimension empty box.
    */
   virtual DatabaseBox getDatabaseBox(const string& key);

   /**
    * Returns a zero dimension empty box.
    */
   virtual DatabaseBox getDatabaseBoxWithDefault(
      const string& key, const DatabaseBox& defaultvalue);

   /**
    * Returns an empty Array<box>.
    */
   virtual Array<DatabaseBox> getDatabaseBoxArray(const string& key);

   /**
    * Does nothing.
    */
   virtual void getDatabaseBoxArray(
      const string& key, DatabaseBox* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isChar(const string& key);

   /**
    * Does nothing.
    */
   virtual void putChar(const string& key, const char& data);

   /**
    * Does nothing.
    */
   virtual void putCharArray(
      const string& key, const Array<char>& data);

   /**
    * Does nothing.
    */
   virtual void putCharArray(
      const string& key, const char* const data, const int nelements);

   /**
    * Always returns 0.
    */
   virtual char getChar(const string& key);

   /**
    * Always returns 0.
    */
   virtual char getCharWithDefault(const string& key, const char& defaultvalue);

   /**
    * Returns an empty Array<char>.
    */
   virtual Array<char> getCharArray(const string& key);

   /**
    * Does nothing.
    */
   virtual void getCharArray(
      const string& key, char* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isComplex(const string& key);

   /**
    * Does nothing.
    */
   virtual void putComplex(const string& key, const dcomplex& data);

   /**
    * Does nothing.
    */
   virtual void putComplexArray(
      const string& key, const Array<dcomplex>& data);

   /**
    * Does nothing.
    */
   virtual void putComplexArray(
      const string& key, const dcomplex* const data, const int nelements);

   /**
    * Returns a 0.0 + 0.0i
    */
   virtual dcomplex getComplex(const string& key);

   /**
    * Returns a 0.0 + 0.0i
    */
   virtual dcomplex getComplexWithDefault(
      const string& key, const dcomplex& defaultvalue);


   /**
    * Returns an empty Array<dcomplex>.
    */
   virtual Array<dcomplex> getComplexArray(const string& key);

   /** 
    * Does nothing.
    */
   virtual void getComplexArray(
      const string& key, dcomplex* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isDouble(const string& key);

   /**
    * Does nothing.
    */
   virtual void putDouble(const string& key, const double& data);

   /**
    * Does nothing.
    */
   virtual void putDoubleArray(
      const string& key, const Array<double>& data);

   /**
    * Does nothing.
    */
   virtual void putDoubleArray(
      const string& key, const double* const data, const int nelements);

   /**
    * Returns 0.0
    */
   virtual double getDouble(const string& key);

   /**
    * Returns 0.0
    */
   virtual double getDoubleWithDefault(
      const string& key, const double& defaultvalue);

   /**
    * Returns an empty Array<double>.
    */
   virtual Array<double> getDoubleArray(const string& key);

   /**
    * Does nothing.
    */
   virtual void getDoubleArray(
      const string& key, double* data, const int nelements);

   /**
    * Always return true.
    */
   virtual bool isFloat(const string& key);

   /**
    * Does nothing.
    */
   virtual void putFloat(const string& key, const float& data);

   /**
    * Does nothing.
    */
   virtual void putFloatArray(
      const string& key, const Array<float>& data);

   /**
    * Does nothing.
    */
   virtual void putFloatArray(
      const string& key, const float* const data, const int nelements);

   /**
    * Returns 0.0
    */
   virtual float getFloat(const string& key);

   /**
    * Returns 0.0
    */
   virtual float getFloatWithDefault(
      const string& key, const float& defaultvalue);

   /**
    * Returns an empty Array<float>.
    */
   virtual Array<float> getFloatArray(const string& key);

   /**
    * Does nothing.
    */
   virtual void getFloatArray(
      const string& key, float* data, const int nelements);

   /**
    * Always returns true. 
    */
   virtual bool isInteger(const string& key);

   /**
    * Does nothing.
    */
   virtual void putInteger(const string& key, const int& data);

   /**
    * Does nothing.
    */
   virtual void putIntegerArray(
      const string& key, const Array<int>& data);

   /**
    * Does nothing.
    */
   virtual void putIntegerArray(
      const string& key, const int* const data, const int nelements);

   /**
    * Returns 0.
    */
   virtual int getInteger(const string& key);

   /**
    * Returns 0.
    */
   virtual int getIntegerWithDefault(
      const string& key, const int& defaultvalue);

   /**
    * Returns an empty Array<int>.
    */
   virtual Array<int> getIntegerArray(const string& key);

   /**
    * Does nothing.
    */
   virtual void getIntegerArray(
      const string& key, int* data, const int nelements);

   /**
    * Always returns true.
    */
   virtual bool isString(const string& key);

   /**
    * Does nothing.
    */
   virtual void putString(const string& key, const string& data);

   /**
    * Does nothing.
    */
   virtual void putStringArray(
      const string& key, const Array<string>& data);

   /**
    * Does nothing.
    */
   virtual void putStringArray(
      const string& key, const string* const data, const int nelements);

   /**
    * Returns and empty string.
    */
   virtual string getString(const string& key);

   /**
    * Returns and empty string.
    */
   virtual string getStringWithDefault(
      const string& key, const string& defaultvalue);

   /**
    * Returns an empty Array<string>.
    */
   virtual Array<string> getStringArray(const string& key);

   /**
    * Does nothing.
    */
   virtual void getStringArray(
      const string& key, string* data, const int nelements);

   /**
    * Does nothing.
    */
   virtual void printClassData(ostream& os = pout);

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
