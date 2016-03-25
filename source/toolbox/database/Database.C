//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/database/Database.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2128 $
// Modified:	$LastChangedDate: 2008-04-11 15:29:55 -0700 (Fri, 11 Apr 2008) $
// Description:	An abstract base class for the SAMRAI database objects
//

#include "tbox/Database.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Database.I"
#endif

#include "tbox/Utilities.h"

namespace SAMRAI {
   namespace tbox {

Database::~Database()
{
}



/*  
 * Boolean
 */ 

/*
 * Create a boolean scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *                                                                      
*/

void Database::putBool(
   const std::string& key, 
   const bool& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   putBoolArray(key, &data, 1);
}

/*
 * Create a boolean array entry in the database with the specified       
 * key 
 */

void Database::putBoolArray(
   const std::string& key, 
   const Array<bool>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putBoolArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putBoolArray() error in database "
		 << getName()
		 << "\n    Attempt to put zero-length array with key = "
		 << key << std::endl);
   }
}

/*                                                                      
 * Get boolean scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a boolean type.
 */

bool Database::getBool(const std::string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   bool ret_val;
   getBoolArray(key, &ret_val, 1);

   return(ret_val);
}

/*
 * Get boolean scalar entry from the database with the specified key  
 * name. An error message is printed and the program exits if the      
 * specified key does not exist in the database or is not associated   
 * with a boolean type.                                                
 */

bool Database::getBoolWithDefault(
   const std::string& key, 
   const bool& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<bool> local_bool = getBoolArray(key);
      bool *locptr = local_bool.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }
}

void Database::getBoolArray(
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
      TBOX_ERROR("Database::getBoolArray() error in database "
		 << getName()
		 << "\n    Incorrect array size = " << nelements
		 << " given for key = " << key
		 << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}

void Database::getScalar(const std::string& key, bool& scalar)
{
   scalar = getBool(key);
}

void Database::putScalar(const std::string& key, const bool scalar)
{
   putBool(key, scalar);
}

void Database::getArray(const std::string& key, Array<bool>& array)
{
   array = getBoolArray(key);
}

void Database::putArray(const std::string& key, const Array<bool> array)
{
   putBoolArray(key, array);
}

/*
*************************************************************************
*                                                                       *
* Create a box entry in the database with the specified                 *
* key name.  A box entry is an array of one.                            *
*                                                                       *
*************************************************************************
*/

void Database::putDatabaseBox(
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

void Database::putDatabaseBoxArray(
   const std::string& key, 
   const Array<DatabaseBox>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putDatabaseBoxArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putDatabaseBoxArray() error in database "
		 << getName()
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

DatabaseBox Database::getDatabaseBox(const std::string& key)
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

DatabaseBox Database::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<DatabaseBox> local_box = getDatabaseBoxArray(key);
      DatabaseBox *locptr = local_box.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}


void Database::getDatabaseBoxArray(
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
      TBOX_ERROR("Database::getDatabaseBoxArray() error in database "
		 << getName()
		 << "\n    Incorrect array size = " << nelements
		 << " given for key = " << key
		 << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*  
 * Char
 */ 

/*
*************************************************************************
*                                                                       *
* Create a char scalar entry in the database with the specified         *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putChar(
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

void Database::putCharArray(
   const std::string& key, 
   const Array<char>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putCharArray(key, data.getPointer(), data.getSize());
   } else { 
      TBOX_ERROR("Database::putCharArray() error in database "
         << getName()
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

char Database::getChar(const std::string& key)
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

char Database::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<char> local_char = getCharArray(key);
      char *locptr = local_char.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}


void Database::getCharArray(
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
      TBOX_ERROR("Database::getCharArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

void Database::getScalar(const std::string& key, char& scalar)
{
   scalar = getChar(key);
}

void Database::putScalar(const std::string& key, const char scalar)
{
   putChar(key, scalar);
}

void Database::getArray(const std::string& key, Array<char>& array)
{
   array = getCharArray(key);
}

void Database::putArray(const std::string& key, const Array<char> array)
{
   putCharArray(key, array);
}

/*  
 * Complex
 */ 

/*
*************************************************************************
*                                                                       *
* Create a complex scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putComplex(
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

void Database::putComplexArray(
   const std::string& key,
   const Array<dcomplex>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putComplexArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putComplexArray() error in database "
         << getName()
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

dcomplex Database::getComplex(const std::string& key)
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

dcomplex Database::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<dcomplex> local_dcomplex = getComplexArray(key);
      dcomplex *locptr = local_dcomplex.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);   
   } else {
      return defaultvalue;
   }

}

void Database::getComplexArray(
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
      TBOX_ERROR("Database::getComplexArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}




void Database::getScalar(const std::string& key, dcomplex& scalar)
{
   scalar = getComplex(key);
}

void Database::putScalar(const std::string& key, const dcomplex scalar)
{
   putComplex(key, scalar);
}

void Database::getArray(const std::string& key, Array<dcomplex>& array)
{
   array = getComplexArray(key);
}

void Database::putArray(const std::string& key, const Array<dcomplex> array)
{
   putComplexArray(key, array);
}

/*  
 * Float
 */ 

/*
*************************************************************************
*                                                                       *
* Create a float scalar entry in the database with the specified        *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putFloat(
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

void Database::putFloatArray(
   const std::string& key,
   const Array<float>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putFloatArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putFloatArray() error in database "
         << getName()
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

float Database::getFloat(const std::string& key)
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

float Database::getFloatWithDefault(
   const std::string& key, 
   const float& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<float> local_float = getFloatArray(key);
      float *locptr = local_float.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}


void Database::getFloatArray(
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
      TBOX_ERROR("Database::getFloatArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

void Database::getScalar(const std::string& key, float& scalar)
{
   scalar = getFloat(key);
}

void Database::putScalar(const std::string& key, const float scalar)
{
   putFloat(key, scalar);
}

void Database::getArray(const std::string& key, Array<float>& array)
{
   array = getFloatArray(key);
}

void Database::putArray(const std::string& key, const Array<float> array)
{
   putFloatArray(key, array);
}

/*  
 * Double
 */ 


/*
*************************************************************************
*                                                                       *
* Create a double scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putDouble(
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

void Database::putDoubleArray(
   const std::string& key,
   const Array<double>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putDoubleArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putDoubleArray() error in database "
         << getName()
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

double Database::getDouble(const std::string& key)
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

double Database::getDoubleWithDefault(
   const std::string& key, 
   const double& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<double> local_double = getDoubleArray(key);
      double *locptr = local_double.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr); 
   } else {
      return defaultvalue;
   }
}

void Database::getDoubleArray(
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
      TBOX_ERROR("Database::getDoubleArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}



void Database::getScalar(const std::string& key, double& scalar)
{
   scalar = getDouble(key);
}

void Database::putScalar(const std::string& key, const double scalar)
{
   putDouble(key, scalar);
}

void Database::getArray(const std::string& key, Array<double>& array)
{
   array = getDoubleArray(key);
}

void Database::putArray(const std::string& key, const Array<double> array)
{
   putDoubleArray(key, array);
}

/*  
 * Integer
 */ 

/*
*************************************************************************
*                                                                       *
* Create a integer scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putInteger(
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

void Database::putIntegerArray(
   const std::string& key, 
   const Array<int>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putIntegerArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putIntegerArray() error in database "
         << getName()
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

int Database::getInteger(const std::string& key)
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

int Database::getIntegerWithDefault(
   const std::string& key, 
   const int& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<int> local_int = getIntegerArray(key);
      int *locptr = local_int.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}

void Database::getIntegerArray(
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
      TBOX_ERROR("Database::getIntegerArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

void Database::getScalar(const std::string& key, int& scalar)
{
   scalar = getInteger(key);
}

void Database::putScalar(const std::string& key, const int scalar)
{
   putInteger(key, scalar);
}

void Database::getArray(const std::string& key, Array<int>& array)
{
   array = getIntegerArray(key);
}

void Database::putArray(const std::string& key, const Array<int> array)
{
   putIntegerArray(key, array);
}

/*
 * String
 */ 

/*
*************************************************************************
*                                                                       *
* Create a string scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void Database::putString(
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

void Database::putStringArray(
   const std::string& key, 
   const Array<std::string>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if ( data.getSize() > 0 ) {
      putStringArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putStringArray() error in database "
         << getName()
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

std::string Database::getString(const std::string& key)
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

std::string Database::getStringWithDefault(
   const std::string& key, 
   const std::string& defaultvalue)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!key.empty());
#endif
   if(keyExists(key)) {
      Array<std::string> local_string = getStringArray(key);
      std::string *locptr = local_string.getPointer();
      return ( locptr == NULL ? defaultvalue : *locptr);
   } else {
      return defaultvalue;
   }

}

void Database::getStringArray(
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
      TBOX_ERROR("Database::getStringArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}

}
}

