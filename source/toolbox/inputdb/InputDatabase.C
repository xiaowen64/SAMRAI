//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/inputdb/InputDatabase.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	An input database structure that stores (key,value) pairs
//

#include "tbox/InputDatabase.h"

#include "tbox/Utilities.h"
#include "tbox/IOStream.h"

#include <stdlib.h>

#ifdef DEBUG_NO_INLINE
#include "tbox/InputDatabase.I"
#endif
#include "tbox/SAMRAI_MPI.h"

#define KEY_DATABASE      (0)
#define KEY_BOOL_ARRAY    (1)
#define KEY_BOX_ARRAY     (2)
#define KEY_CHAR_ARRAY    (3)
#define KEY_COMPLEX_ARRAY (4)
#define KEY_DOUBLE_ARRAY  (5)
#define KEY_FLOAT_ARRAY   (6)
#define KEY_INT_ARRAY     (7)
#define KEY_STRING_ARRAY  (8)

#define PRINT_DEFAULT (1)
#define PRINT_INPUT   (2)
#define PRINT_UNUSED  (4)

#define SSTREAM_BUFFER (4096)

#define INPUT_DB_ERROR(X) do {						\
      pout << "InputDatabase: " << X << std::endl << std::flush;			\
      printClassData(pout);						\
      pout << "Program abort called..." << std::endl << std::flush;               \
      SAMRAI_MPI::abort(); 						\
} while (0)

namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*									*
* The virtual destructor deallocates database data.			*
*									*
*************************************************************************
*/

InputDatabase::~InputDatabase()
{
}

/*
*************************************************************************
*									*
* Return whether the key exists in the database.			*
*									*
*************************************************************************
*/

bool InputDatabase::keyExists(const std::string& key)
{
   return(findKeyData(key) ? true : false);
}

/*
*************************************************************************
*									*
* Return all of the keys in the database.				*
*									*
*************************************************************************
*/

Array<std::string> InputDatabase::getAllKeys()
{
   const int n = d_keyvalues.getNumberItems();
   Array<std::string> keys(n);

   int k = 0;
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {
      keys[k++] = i().d_key;
   }

   return(keys);
}

/*
*************************************************************************
*									*
* Get the size of the array entry associated with the specified key;	*
* return 0 if the key does not exist.					*
*									*
*************************************************************************
*/

int InputDatabase::getArraySize(const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return(keydata ? keydata->d_array_size : 0);
}

/*
*************************************************************************
*									*
* Member functions that manage the database values within the database.	*
*									*
*************************************************************************
*/

bool InputDatabase::isDatabase(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(keydata ? keydata->d_type == KEY_DATABASE : false);
}

Pointer<Database> InputDatabase::putDatabase(const std::string& key)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_DATABASE;
   keydata.d_array_size   = 1;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_database     = new InputDatabase(key);
   d_keyvalues.appendItem(keydata);
   return(keydata.d_database);
}

Pointer<Database> InputDatabase::getDatabase(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if (keydata->d_type != KEY_DATABASE) {
      INPUT_DB_ERROR("Key=" << key << " is not a database...");
   }
   keydata->d_accessed = true;
   return(keydata->d_database);
}

/*
*************************************************************************
*									*
* Member functions that manage boolean values within the database.	*
*									*
*************************************************************************
*/

bool InputDatabase::isBool(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(keydata ? keydata->d_type == KEY_BOOL_ARRAY : false);
}

void InputDatabase::putBool(const std::string& key, const bool& data)
{
   putBoolArray(key, &data, 1);
}

void InputDatabase::putBoolArray(
   const std::string& key, const Array<bool>& data)
{
   this->putBoolArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putBoolArray(
   const std::string& key, const bool* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_BOOL_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_boolean      = Array<bool>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_boolean[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

bool InputDatabase::getBool(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != KEY_BOOL_ARRAY) || (keydata->d_array_size != 1)) {
      INPUT_DB_ERROR("Key=" << key << " is not a boolean scalar...");
   }
   keydata->d_accessed = true;
   return(keydata->d_boolean[0]);
}

bool InputDatabase::getBoolWithDefault(
   const std::string& key, const bool& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getBool(key));
   putBool(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<bool> InputDatabase::getBoolArray(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if (keydata->d_type != KEY_BOOL_ARRAY) {
      INPUT_DB_ERROR("Key=" << key << " is not a boolean...");
   }
   keydata->d_accessed = true;
   return(keydata->d_boolean);
}

void InputDatabase::getBoolArray(
   const std::string& key, bool* data, const int nelements)
{
   Array<bool> tmp = this->getBoolArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage box values within the database.		*
*									*
*************************************************************************
*/

bool InputDatabase::isDatabaseBox(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(keydata ? keydata->d_type == KEY_BOX_ARRAY : false);
}

void InputDatabase::putDatabaseBox(const std::string& key, const DatabaseBox& data)
{
   putDatabaseBoxArray(key, &data, 1);
}

void InputDatabase::putDatabaseBoxArray(
   const std::string& key, const Array<DatabaseBox>& data)
{
   this->putDatabaseBoxArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putDatabaseBoxArray(
   const std::string& key, const DatabaseBox* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_BOX_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_box          = Array<DatabaseBox>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_box[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

DatabaseBox InputDatabase::getDatabaseBox(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != KEY_BOX_ARRAY) || (keydata->d_array_size != 1)) {
      INPUT_DB_ERROR("Key=" << key << " is not a single box...");
   }
   keydata->d_accessed = true;
   return(keydata->d_box[0]);
}

DatabaseBox InputDatabase::getDatabaseBoxWithDefault(
   const std::string& key, const DatabaseBox& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getDatabaseBox(key));
   putDatabaseBox(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<DatabaseBox> InputDatabase::getDatabaseBoxArray(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if (keydata->d_type != KEY_BOX_ARRAY) {
      INPUT_DB_ERROR("Key=" << key << " is not a box...");
   }
   keydata->d_accessed = true;
   return(keydata->d_box);
}

void InputDatabase::getDatabaseBoxArray(
   const std::string& key, DatabaseBox* data, const int nelements)
{
   Array<DatabaseBox> tmp = this->getDatabaseBoxArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage character values within the database.	*
*									*
*************************************************************************
*/

bool InputDatabase::isChar(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(keydata ? keydata->d_type == KEY_CHAR_ARRAY : false);
}

void InputDatabase::putChar(const std::string& key, const char& data)
{
   putCharArray(key, &data, 1);
}

void InputDatabase::putCharArray(
   const std::string& key, const Array<char>& data)
{
   this->putCharArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putCharArray(
   const std::string& key, const char* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_CHAR_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_char         = Array<char>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_char[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

char InputDatabase::getChar(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != KEY_CHAR_ARRAY) || (keydata->d_array_size != 1)) {
      INPUT_DB_ERROR("Key=" << key << " is not a single character...");
   }
   keydata->d_accessed = true;
   return(keydata->d_char[0]);
}

char InputDatabase::getCharWithDefault(
   const std::string& key, const char& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getChar(key));
   putChar(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<char> InputDatabase::getCharArray(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if (keydata->d_type != KEY_CHAR_ARRAY) {
      INPUT_DB_ERROR("Key=" << key << " is not a character...");
   }
   keydata->d_accessed = true;
   return(keydata->d_char);
}

void InputDatabase::getCharArray(
   const std::string& key, char* data, const int nelements)
{
   Array<char> tmp = this->getCharArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage complex values within the database.	*
* Note that complex numbers may be promoted from integers, floats,	*
* and doubles.								*
*									*
*************************************************************************
*/

bool InputDatabase::isComplex(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(!keydata ? false : (keydata->d_type == KEY_COMPLEX_ARRAY
                           || keydata->d_type == KEY_INT_ARRAY
                           || keydata->d_type == KEY_FLOAT_ARRAY
                           || keydata->d_type == KEY_DOUBLE_ARRAY));
}

void InputDatabase::putComplex(const std::string& key, const dcomplex& data)
{
   putComplexArray(key, &data, 1);
}

void InputDatabase::putComplexArray(
   const std::string& key, const Array<dcomplex>& data)
{
   this->putComplexArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putComplexArray(
   const std::string& key, const dcomplex* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_COMPLEX_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_complex      = Array<dcomplex>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_complex[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

dcomplex InputDatabase::getComplex(const std::string& key)
{
   dcomplex value(0.0,0.0);
   KeyData *keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      INPUT_DB_ERROR("Key=" << key << " is not a single complex...");
   }

   switch (keydata->d_type) {
      case KEY_INT_ARRAY:
         value = dcomplex((double) keydata->d_integer[0], 0.0);
         break;
      case KEY_FLOAT_ARRAY:
         value = dcomplex((double) keydata->d_float[0], 0.0);
         break;
      case KEY_DOUBLE_ARRAY:
         value = dcomplex(keydata->d_double[0], 0.0);
         break;
      case KEY_COMPLEX_ARRAY:
         value = keydata->d_complex[0];
         break;
      default:
         INPUT_DB_ERROR("Key=" << key << " is not a single complex...");
   }

   keydata->d_accessed = true;
   return(value);
}

dcomplex InputDatabase::getComplexWithDefault(
   const std::string& key, const dcomplex& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getComplex(key));
   putComplex(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<dcomplex> InputDatabase::getComplexArray(
   const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   Array<dcomplex> array;
   switch (keydata->d_type) {
      case KEY_INT_ARRAY: {
         array = Array<dcomplex>(keydata->d_integer.getSize());
         for (int i = 0; i < keydata->d_integer.getSize(); i++) {
            array[i] = dcomplex((double) keydata->d_integer[i], 0.0);
         }
         break;
      }
      case KEY_FLOAT_ARRAY: {
         array = Array<dcomplex>(keydata->d_float.getSize());
         for (int i = 0; i < keydata->d_float.getSize(); i++) {
            array[i] = dcomplex((double) keydata->d_float[i], 0.0);
         }
         break;
      }
      case KEY_DOUBLE_ARRAY: {
         array = Array<dcomplex>(keydata->d_double.getSize());
         for (int i = 0; i < keydata->d_float.getSize(); i++) {
            array[i] = dcomplex(keydata->d_double[i], 0.0);
         }
         break;
      }
      case KEY_COMPLEX_ARRAY:
         array = keydata->d_complex;
         break;
      default:
         INPUT_DB_ERROR("Key=" << key << " is not a complex...");
   }
   keydata->d_accessed = true;
   return(array);
}

void InputDatabase::getComplexArray(
   const std::string& key, dcomplex* data, const int nelements)
{
   Array<dcomplex> tmp = this->getComplexArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage double values within the database.	*
* Note that doubles may be promoted from integers or floats.		*
*									*
*************************************************************************
*/

bool InputDatabase::isDouble(
   const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(!keydata ? false : (keydata->d_type == KEY_DOUBLE_ARRAY
                           || keydata->d_type == KEY_INT_ARRAY
                           || keydata->d_type == KEY_FLOAT_ARRAY));
}

void InputDatabase::putDouble(
   const std::string& key, const double& data)
{
   putDoubleArray(key, &data, 1);
}

void InputDatabase::putDoubleArray(
   const std::string& key, const Array<double>& data)
{
   this->putDoubleArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putDoubleArray(
   const std::string& key, const double* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_DOUBLE_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_double       = Array<double>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_double[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

double InputDatabase::getDouble(const std::string& key)
{
   double value = 0.0;
   KeyData *keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      INPUT_DB_ERROR("Key=" << key << " is not a single double...");
   }

   switch (keydata->d_type) {
      case KEY_INT_ARRAY:
         value = (double) keydata->d_integer[0];
         break;
      case KEY_FLOAT_ARRAY:
         value = (double) keydata->d_float[0];
         break;
      case KEY_DOUBLE_ARRAY:
         value = keydata->d_double[0];
         break;
      default:
         INPUT_DB_ERROR("Key=" << key << " is not a single double...");
   }

   keydata->d_accessed = true;
   return(value);
}

double InputDatabase::getDoubleWithDefault(
   const std::string& key, const double& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getDouble(key));
   putDouble(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<double> InputDatabase::getDoubleArray(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   Array<double> array;
   switch (keydata->d_type) {
      case KEY_INT_ARRAY: {
         array = Array<double>(keydata->d_integer.getSize());
         for (int i = 0; i < keydata->d_integer.getSize(); i++) {
            array[i] = (double) keydata->d_integer[i];
         }
         break;
      }
      case KEY_FLOAT_ARRAY: {
         array = Array<double>(keydata->d_float.getSize());
         for (int i = 0; i < keydata->d_float.getSize(); i++) {
            array[i] = (double) keydata->d_float[i];
         }
         break;
      }
      case KEY_DOUBLE_ARRAY: {
         array = keydata->d_double;
         break;
      }
      default:
         INPUT_DB_ERROR("Key=" << key << " is not a double...");
   }
   keydata->d_accessed = true;
   return(array);
}

void InputDatabase::getDoubleArray(
   const std::string& key, double* data, const int nelements)
{
   Array<double> tmp = this->getDoubleArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage float values within the database.	*
* Note that floats may be promoted from integers or truncated from	*
* doubles (without a warning).						*
*									*
*************************************************************************
*/

bool InputDatabase::isFloat(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(!keydata ? false : (keydata->d_type == KEY_DOUBLE_ARRAY
                           || keydata->d_type == KEY_INT_ARRAY
                           || keydata->d_type == KEY_FLOAT_ARRAY));
}

void InputDatabase::putFloat(const std::string& key, const float& data)
{
   putFloatArray(key, &data, 1);
}

void InputDatabase::putFloatArray(
   const std::string& key, const Array<float>& data)
{
   this->putFloatArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putFloatArray(
   const std::string& key, const float* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_FLOAT_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_float        = Array<float>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_float[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

float InputDatabase::getFloat(
   const std::string& key)
{
   float value = 0.0;
   KeyData *keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      INPUT_DB_ERROR("Key=" << key << " is not a single float...");
   }

   switch (keydata->d_type) {
      case KEY_INT_ARRAY:
         value = static_cast<float>( keydata->d_integer[0] );
         break;
      case KEY_FLOAT_ARRAY:
         value = keydata->d_float[0];
         break;
      case KEY_DOUBLE_ARRAY:
         value = static_cast<float>( keydata->d_double[0] );
         break;
      default:
         INPUT_DB_ERROR("Key=" << key << " is not a single float...");
   }

   keydata->d_accessed = true;
   return(value);
}

float InputDatabase::getFloatWithDefault(
   const std::string& key, const float& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getFloat(key));
   putFloat(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<float> InputDatabase::getFloatArray(
   const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   Array<float> array;
   switch (keydata->d_type) {
      case KEY_INT_ARRAY: {
         array = Array<float>(keydata->d_integer.getSize());
         for (int i = 0; i < keydata->d_integer.getSize(); i++) {
            array[i] = static_cast<float>( keydata->d_integer[i] );
         }
         break;
      }
      case KEY_FLOAT_ARRAY:
         array = keydata->d_float;
         break;
      case KEY_DOUBLE_ARRAY: {
         array = Array<float>(keydata->d_double.getSize());
         for (int i = 0; i < keydata->d_double.getSize(); i++) {
            array[i] = static_cast<float>( keydata->d_double[i] );
         }
         break;
      }
      default:
         INPUT_DB_ERROR("Key=" << key << " is not a float...");
   }
   keydata->d_accessed = true;
   return(array);
}

void InputDatabase::getFloatArray(
   const std::string& key, float* data, const int nelements)
{
   Array<float> tmp = this->getFloatArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage integer values within the database.	*
*									*
*************************************************************************
*/

bool InputDatabase::isInteger(
   const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(!keydata ? false : keydata->d_type == KEY_INT_ARRAY);
}

void InputDatabase::putInteger(
   const std::string& key, const int& data)
{
   putIntegerArray(key, &data, 1);
}

void InputDatabase::putIntegerArray(
   const std::string& key, const Array<int>& data)
{
   this->putIntegerArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putIntegerArray(
   const std::string& key, const int* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_INT_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_integer      = Array<int>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_integer[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

int InputDatabase::getInteger(
   const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != KEY_INT_ARRAY) || (keydata->d_array_size != 1)) {
      INPUT_DB_ERROR("Key=" << key << " is not an integer scalar...");
   }
   keydata->d_accessed = true;
   return(keydata->d_integer[0]);
}

int InputDatabase::getIntegerWithDefault(
   const std::string& key, const int& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getInteger(key));
   putInteger(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<int> InputDatabase::getIntegerArray(
   const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if (keydata->d_type != KEY_INT_ARRAY) {
      INPUT_DB_ERROR("Key=" << key << " is not an integer...");
   }
   keydata->d_accessed = true;
   return(keydata->d_integer);
}

void InputDatabase::getIntegerArray(
   const std::string& key, int* data, const int nelements)
{
   Array<int> tmp = this->getIntegerArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Member functions that manage string values within the database.	*
*									*
*************************************************************************
*/

bool InputDatabase::isString(const std::string& key)
{
   KeyData *keydata = findKeyData(key);
   return(!keydata ? false : keydata->d_type == KEY_STRING_ARRAY);
}

void InputDatabase::putString(
   const std::string& key, const std::string& data)
{
   putStringArray(key, &data, 1);
}

void InputDatabase::putStringArray(
   const std::string& key, 
   const Array<std::string>& data)
{
   this->putStringArray(key, data.getPointer(), data.getSize());
}

void InputDatabase::putStringArray(
   const std::string& key, const std::string* const data, const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key          = key;
   keydata.d_type         = KEY_STRING_ARRAY;
   keydata.d_array_size   = nelements;
   keydata.d_accessed     = false;
   keydata.d_from_default = false;
   keydata.d_string       = Array<std::string>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_string[i] = data[i];
   }

   d_keyvalues.appendItem(keydata);
}

std::string InputDatabase::getString(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != KEY_STRING_ARRAY) || (keydata->d_array_size != 1)) {
      INPUT_DB_ERROR("Key=" << key << " is not a single string...");
   }
   keydata->d_accessed = true;
   return(keydata->d_string[0]);
}

std::string InputDatabase::getStringWithDefault(
   const std::string& key, const std::string& defaultvalue)
{
   KeyData *keydata = findKeyData(key);
   if (keydata) return(this->getString(key));
   putString(key, defaultvalue);
   d_keyvalues.getLastItem().d_from_default = true;
   return(defaultvalue);
}

Array<std::string> InputDatabase::getStringArray(const std::string& key)
{
   KeyData *keydata = findKeyDataOrExit(key);
   if (keydata->d_type != KEY_STRING_ARRAY) {
      INPUT_DB_ERROR("Key=" << key << " is not a string...");
   }
   keydata->d_accessed = true;
   return(keydata->d_string);
}

void InputDatabase::getStringArray(
   const std::string& key, std::string* data, const int nelements)
{
   Array<std::string> tmp = this->getStringArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      INPUT_DB_ERROR("Incorrect array size=" << nelements << " specified for key="
            << key << " with array size=" << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*									*
* Search the current database for a matching key.  If found, delete	*
* that key and return true.  If the key does not exist, then return	*
* false.								*
*									*
*************************************************************************
*/

bool InputDatabase::deleteKeyIfFound(const std::string& key)
{
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {
      if (i().d_key == key) {
         d_keyvalues.removeItem(i);
         return(true);
      }
   }
   return(false);
}

/*
*************************************************************************
*									*
* Find the key data associated with the specified key and return a	*
* pointer to the record.  If no such key data exists, then return NULL.	*
*									*
*************************************************************************
*/

InputDatabase::KeyData*
InputDatabase::findKeyData(const std::string& key)
{
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {
      if (key == i().d_key) return(&i());
   }
   return(NULL);
}

/*
*************************************************************************
*									*
* Find the key data associated with the specified key and return a	*
* pointer to the record.  If no such key data exists, then exit with	*
* an error message.							*
*									*
*************************************************************************
*/

InputDatabase::KeyData*
InputDatabase::findKeyDataOrExit(const std::string& key)
{
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {
      if (key == i().d_key) return(&i());
   }
   INPUT_DB_ERROR("Key ``" << key << "'' does not exist in the database...");
   return(NULL);
}

/*
*************************************************************************
*									*
* Print the entire input database to the specified output stream.	*
*									*
*************************************************************************
*/

void InputDatabase::printClassData(std::ostream& os)
{
   printDatabase(os, 0, PRINT_DEFAULT | PRINT_INPUT | PRINT_UNUSED);
}

/*
*************************************************************************
*									*
* Print unused input database keys to the specified output stream.	*
*									*
*************************************************************************
*/

void InputDatabase::printUnusedKeys(std::ostream& os) const
{
   printDatabase(os, 0, PRINT_UNUSED);
}

/*
*************************************************************************
*									*
* Print default input database keys to the specified output stream.	*
*									*
*************************************************************************
*/

void InputDatabase::printDefaultKeys(std::ostream& os) const
{
   printDatabase(os, 0, PRINT_DEFAULT);
}

/*
*************************************************************************
*									*
* Indent the output stream by the specified indentation factor.		*
*									*
*************************************************************************
*/

void InputDatabase::indentStream(std::ostream& os, const int indent)
{
   for (int i = 0; i < indent; i++) {
      os << " ";
   }
}

/*
*************************************************************************
*									*
* Print database data to the specified output stream.			*
*									*
*************************************************************************
*/

void InputDatabase::printDatabase(
   std::ostream& os, const int indent, const int toprint) const
{
   /*
    * Get the maximum key width in the output (excluding databases)
    */

   int width = 0;
   for (List<KeyData>::Iterator k(d_keyvalues); k; k++) {
      if ( ( (k().d_from_default) && (toprint & PRINT_DEFAULT))
        || ( (k().d_accessed)     && (toprint & PRINT_INPUT  ))
        || (!(k().d_accessed)     && (toprint & PRINT_UNUSED ))) {
         if (k().d_type != KEY_DATABASE) {
            const int keywidth = k().d_key.length();
            if (keywidth > width) width = keywidth;
         }
      }
   }

   /*
    * Iterate over all non-database keys in the database and output key values
    */

   indentStream(os, indent);
   os << d_database_name << " {\n";
   for (List<KeyData>::Iterator i(d_keyvalues); i; i++) {

      if ( ( (i().d_from_default) && (toprint & PRINT_DEFAULT))
        || ( (i().d_accessed)     && (toprint & PRINT_INPUT  ))
        || (!(i().d_accessed)     && (toprint & PRINT_UNUSED ))) {

#ifndef LACKS_SSTREAM
         std::ostringstream sstream;
#else
         char sstream_buffer[SSTREAM_BUFFER];
         std::ostrstream sstream(sstream_buffer, SSTREAM_BUFFER);
#endif

         switch(i().d_type) {
            case KEY_DATABASE: {
               break;
            }

            case KEY_BOOL_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_boolean.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << (i().d_boolean[j] ? "TRUE" : "FALSE");
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_BOX_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_box.getSize();
               for (int j = 0; j < n; j++) {
                  const int m = i().d_box[j].getDimension();
                  sstream << "[(";
                  for (int k = 0; k < m; k++) {
                     sstream << i().d_box[j].lower(k);
                     if (k < m-1) sstream << ",";
                  }
                  sstream << "),(";
                  for (int l = 0; l < m; l++) {
                     sstream << i().d_box[j].upper(l);
                     if (l < m-1) sstream << ",";
                  }
                  sstream << ")]";
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_CHAR_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_char.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << "'" << i().d_char[j] << "'";
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_COMPLEX_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_complex.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i().d_complex[j];
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_DOUBLE_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_double.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i().d_double[j];
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_FLOAT_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_float.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i().d_float[j];
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_INT_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_integer.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i().d_integer[j];
                  if (j < n-1) sstream << ", ";
               }
               break;
            }

            case KEY_STRING_ARRAY: {
               indentStream(sstream, indent+3);
               sstream << i().d_key;
               indentStream(sstream, width-i().d_key.length());
               sstream << " = ";
               const int n = i().d_string.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << "\"" << i().d_string[j] << "\"";
                  if (j < n-1) sstream << ", ";
               }
               break;
            }
         }

         /*
          * Output whether the key was used or default in column 60
          */

         if (i().d_type != KEY_DATABASE) {
#ifndef LACKS_SSTREAM
            const int tab = 59 - sstream.str().length();
#else
            const int tab = 59 - sstream.pcount();
#endif
            if (tab > 0) indentStream(sstream, tab);
            if (i().d_from_default) {
               sstream << " // from default";
            } else if (i().d_accessed) {
               sstream << " // input used";
            } else {
               sstream << " // input not used";
            }
          
//            sstream << std::endl << ends;
            sstream << std::endl;
            os << sstream.str();
         }
      }
   }

   /*
    * Finally, output all databases in the current key list
    */

   for (List<KeyData>::Iterator j(d_keyvalues); j; j++) {
      if (j().d_type == KEY_DATABASE) {
         Pointer<InputDatabase> db = j().d_database;
         db->printDatabase(os, indent+3, toprint);
      }
   }

   indentStream(os, indent);
   os << "}\n";
}

}
}
