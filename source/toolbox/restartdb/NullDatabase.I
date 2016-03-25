//
// File:	NullDatabase.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An input database structure that stores (key,value) pairs
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
NullDatabase::NullDatabase()
{
}

/*
*************************************************************************
*									*
* The virtual destructor deallocates database data.			*
*									*
*************************************************************************
*/

inline
NullDatabase::~NullDatabase()
{
}

/*
*************************************************************************
*									*
* Always returns true.							*
*									*
*************************************************************************
*/

inline
bool NullDatabase::keyExists(const string& key)
{
   NULL_USE(key);
   return true;
}

/*
*************************************************************************
*									*
* Return an empty Array<string>.					*
*									*
*************************************************************************
*/

inline
Array<string> NullDatabase::getAllKeys()
{ 
   Array<string> keys(0);
   return(keys);
}

/*
*************************************************************************
*									*
* Always returns 0.							* 
*									*
*************************************************************************
*/

inline
int NullDatabase::getArraySize(const string& key)
{
   NULL_USE(key);
   return 0;
}

/*
*************************************************************************
*									*
* Member functions that manage the database values within the database.	*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isDatabase(const string& key)
{
   NULL_USE(key);
   return true;
}

inline
Pointer<Database> NullDatabase::putDatabase(const string& key)
{
   NULL_USE(key);
   return (this);
}

inline
Pointer<Database> NullDatabase::getDatabase(const string& key)
{
   NULL_USE(key);
   Pointer<Database> new_NullDatabase = new NullDatabase();
   return (new_NullDatabase);
}

/*
*************************************************************************
*									*
* Member functions that manage boolean values within the database.	*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isBool(const string& key)
{
   NULL_USE(key);
   return true;
}

inline
void NullDatabase::putBool(const string& key, const bool& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putBoolArray(
   const string& key, const Array<bool>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putBoolArray(
   const string& key, const bool* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
bool NullDatabase::getBool(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
bool NullDatabase::getBoolWithDefault(
   const string& key, 
   const bool& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (true);
}

inline
Array<bool> NullDatabase::getBoolArray(const string& key)
{
   NULL_USE(key);
   Array<bool> empty(0);
   return (empty);
}

inline
void NullDatabase::getBoolArray(
   const string& key, bool* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

/*
*************************************************************************
*									*
* Member functions that manage box values within the database.		*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isDatabaseBox(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putDatabaseBox(const string& key, const DatabaseBox& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDatabaseBoxArray(
   const string& key, const Array<DatabaseBox>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDatabaseBoxArray(
   const string& key, const DatabaseBox* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
DatabaseBox NullDatabase::getDatabaseBox(const string& key)
{
   NULL_USE(key);
   DatabaseBox box;
   return(box);
}

inline
DatabaseBox NullDatabase::getDatabaseBoxWithDefault(
   const string& key, 
   const DatabaseBox& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);

   DatabaseBox box;
   return(box);
}

inline
Array<DatabaseBox> NullDatabase::getDatabaseBoxArray(const string& key)
{
   NULL_USE(key);

   Array<DatabaseBox> empty(0);
   return (empty);
}

inline
void NullDatabase::getDatabaseBoxArray(
   const string& key, DatabaseBox* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

/*
*************************************************************************
*									*
* Member functions that manage character values within the database.	*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isChar(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putChar(const string& key, const char& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putCharArray(
   const string& key, const Array<char>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putCharArray(
   const string& key, const char* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
char NullDatabase::getChar(const string& key)
{
   NULL_USE(key);
   return (0);
}

inline
char NullDatabase::getCharWithDefault(
   const string& key, 
   const char& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (0);
}

inline
Array<char> NullDatabase::getCharArray(const string& key)
{
   NULL_USE(key);

   Array<char> empty(0);
   return (empty);
}

inline
void NullDatabase::getCharArray(
   const string& key, char* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
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

inline
bool NullDatabase::isComplex(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putComplex(const string& key, const dcomplex& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putComplexArray(
   const string& key, const Array<dcomplex>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putComplexArray(
   const string& key, const dcomplex* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
dcomplex NullDatabase::getComplex(const string& key)
{
   NULL_USE(key);

   dcomplex value(0.0,0.0);
   return(value);
}

inline
dcomplex NullDatabase::getComplexWithDefault(
   const string& key, const dcomplex& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);

   dcomplex value(0.0,0.0);
   return(value);
}

inline
Array<dcomplex> NullDatabase::getComplexArray(const string& key)
{
   NULL_USE(key);
  
   Array<dcomplex> empty(0);
   return(empty);
}

inline
void NullDatabase::getComplexArray(
   const string& key, dcomplex* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

/*
*************************************************************************
*									*
* Member functions that manage double values within the database.	*
* Note that doubles may be promoted from integers or floats.		*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isDouble(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putDouble(const string& key, const double& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDoubleArray(
   const string& key, const Array<double>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDoubleArray(
   const string& key, const double* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
double NullDatabase::getDouble(const string& key)
{
   NULL_USE(key);
   return (0.0);
}

inline
double NullDatabase::getDoubleWithDefault(
   const string& key, 
   const double& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (0.0);
}

inline
Array<double> NullDatabase::getDoubleArray(const string& key)
{
   NULL_USE(key);
   Array<double> empty(0);
   return (empty);
}

inline
void NullDatabase::getDoubleArray(
   const string& key, double* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
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

inline
bool NullDatabase::isFloat(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putFloat(const string& key, const float& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putFloatArray(
   const string& key, const Array<float>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putFloatArray(
   const string& key, const float* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
float NullDatabase::getFloat(const string& key)
{
   NULL_USE(key);
   return (0.0);
}

inline
float NullDatabase::getFloatWithDefault(
   const string& key, 
   const float& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (0.0);
}

inline
Array<float> NullDatabase::getFloatArray(const string& key)
{
   NULL_USE(key); 

   Array<float> empty(0);
   return(empty);
}

inline
void NullDatabase::getFloatArray(
   const string& key, float* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

/*
*************************************************************************
*									*
* Member functions that manage integer values within the database.	*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isInteger(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putInteger(const string& key, const int& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putIntegerArray(
   const string& key, const Array<int>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putIntegerArray(
   const string& key, const int* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
int NullDatabase::getInteger(const string& key)
{
   NULL_USE(key);
   return(0);
}

inline
int NullDatabase::getIntegerWithDefault(
   const string& key, 
   const int& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return(0);
}

inline
Array<int> NullDatabase::getIntegerArray(const string& key)
{
   NULL_USE(key);
   
   Array<int> empty(0);
   return (empty);
}

inline
void NullDatabase::getIntegerArray(
   const string& key, 
   int* data, 
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

/*
*************************************************************************
*									*
* Member functions that manage string values within the database.	*
*									*
*************************************************************************
*/

inline
bool NullDatabase::isString(const string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putString(const string& key, const string& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putStringArray(
   const string& key, const Array<string>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putStringArray(
   const string& key, const string* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
string NullDatabase::getString(const string& key)
{
   NULL_USE(key); 
   string empty;
   return (empty);
}

inline
string NullDatabase::getStringWithDefault(
   const string& key, 
   const string& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   string empty;
   return (empty);
}

inline
Array<string> NullDatabase::getStringArray(const string& key)
{
   NULL_USE(key);
   Array<string> empty(0);
   return (empty);
}

inline
void NullDatabase::getStringArray(
   const string& key, string* data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

/*
*************************************************************************
*									*
* Does nothing.								*
*									*
*************************************************************************
*/

inline
void NullDatabase::printClassData(ostream& os)
{
   NULL_USE(os);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
