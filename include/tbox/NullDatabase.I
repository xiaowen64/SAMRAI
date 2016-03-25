//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/restartdb/NullDatabase.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
bool NullDatabase::keyExists(const std::string& key)
{
   NULL_USE(key);
   return true;
}

/*
*************************************************************************
*									*
* Return an empty Array<std::string>.					*
*									*
*************************************************************************
*/

inline
Array<std::string> NullDatabase::getAllKeys()
{ 
   Array<std::string> keys(0);
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
int NullDatabase::getArraySize(const std::string& key)
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
bool NullDatabase::isDatabase(const std::string& key)
{
   NULL_USE(key);
   return true;
}

inline
Pointer<Database> NullDatabase::putDatabase(const std::string& key)
{
   NULL_USE(key);
   return (this);
}

inline
Pointer<Database> NullDatabase::getDatabase(const std::string& key)
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
bool NullDatabase::isBool(const std::string& key)
{
   NULL_USE(key);
   return true;
}

inline
void NullDatabase::putBool(const std::string& key, const bool& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putBoolArray(
   const std::string& key, const Array<bool>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putBoolArray(
   const std::string& key, const bool* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
bool NullDatabase::getBool(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
bool NullDatabase::getBoolWithDefault(
   const std::string& key, 
   const bool& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (true);
}

inline
Array<bool> NullDatabase::getBoolArray(const std::string& key)
{
   NULL_USE(key);
   Array<bool> empty(0);
   return (empty);
}

inline
void NullDatabase::getBoolArray(
   const std::string& key, bool* data, const int nelements)
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
bool NullDatabase::isDatabaseBox(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putDatabaseBox(const std::string& key, const DatabaseBox& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDatabaseBoxArray(
   const std::string& key, const Array<DatabaseBox>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDatabaseBoxArray(
   const std::string& key, const DatabaseBox* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
DatabaseBox NullDatabase::getDatabaseBox(const std::string& key)
{
   NULL_USE(key);
   DatabaseBox box;
   return(box);
}

inline
DatabaseBox NullDatabase::getDatabaseBoxWithDefault(
   const std::string& key, 
   const DatabaseBox& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);

   DatabaseBox box;
   return(box);
}

inline
Array<DatabaseBox> NullDatabase::getDatabaseBoxArray(const std::string& key)
{
   NULL_USE(key);

   Array<DatabaseBox> empty(0);
   return (empty);
}

inline
void NullDatabase::getDatabaseBoxArray(
   const std::string& key, DatabaseBox* data, const int nelements)
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
bool NullDatabase::isChar(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putChar(const std::string& key, const char& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putCharArray(
   const std::string& key, const Array<char>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putCharArray(
   const std::string& key, const char* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
char NullDatabase::getChar(const std::string& key)
{
   NULL_USE(key);
   return (0);
}

inline
char NullDatabase::getCharWithDefault(
   const std::string& key, 
   const char& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (0);
}

inline
Array<char> NullDatabase::getCharArray(const std::string& key)
{
   NULL_USE(key);

   Array<char> empty(0);
   return (empty);
}

inline
void NullDatabase::getCharArray(
   const std::string& key, char* data, const int nelements)
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
bool NullDatabase::isComplex(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putComplex(const std::string& key, const dcomplex& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putComplexArray(
   const std::string& key, const Array<dcomplex>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putComplexArray(
   const std::string& key, const dcomplex* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
dcomplex NullDatabase::getComplex(const std::string& key)
{
   NULL_USE(key);

   dcomplex value(0.0,0.0);
   return(value);
}

inline
dcomplex NullDatabase::getComplexWithDefault(
   const std::string& key, const dcomplex& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);

   dcomplex value(0.0,0.0);
   return(value);
}

inline
Array<dcomplex> NullDatabase::getComplexArray(const std::string& key)
{
   NULL_USE(key);
  
   Array<dcomplex> empty(0);
   return(empty);
}

inline
void NullDatabase::getComplexArray(
   const std::string& key, dcomplex* data, const int nelements)
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
bool NullDatabase::isDouble(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putDouble(const std::string& key, const double& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDoubleArray(
   const std::string& key, const Array<double>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putDoubleArray(
   const std::string& key, const double* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
double NullDatabase::getDouble(const std::string& key)
{
   NULL_USE(key);
   return (0.0);
}

inline
double NullDatabase::getDoubleWithDefault(
   const std::string& key, 
   const double& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (0.0);
}

inline
Array<double> NullDatabase::getDoubleArray(const std::string& key)
{
   NULL_USE(key);
   Array<double> empty(0);
   return (empty);
}

inline
void NullDatabase::getDoubleArray(
   const std::string& key, double* data, const int nelements)
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
bool NullDatabase::isFloat(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putFloat(const std::string& key, const float& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putFloatArray(
   const std::string& key, const Array<float>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putFloatArray(
   const std::string& key, const float* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
float NullDatabase::getFloat(const std::string& key)
{
   NULL_USE(key);
   return (0.0);
}

inline
float NullDatabase::getFloatWithDefault(
   const std::string& key, 
   const float& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return (0.0);
}

inline
Array<float> NullDatabase::getFloatArray(const std::string& key)
{
   NULL_USE(key); 

   Array<float> empty(0);
   return(empty);
}

inline
void NullDatabase::getFloatArray(
   const std::string& key, float* data, const int nelements)
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
bool NullDatabase::isInteger(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putInteger(const std::string& key, const int& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putIntegerArray(
   const std::string& key, const Array<int>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putIntegerArray(
   const std::string& key, const int* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
int NullDatabase::getInteger(const std::string& key)
{
   NULL_USE(key);
   return(0);
}

inline
int NullDatabase::getIntegerWithDefault(
   const std::string& key, 
   const int& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   return(0);
}

inline
Array<int> NullDatabase::getIntegerArray(const std::string& key)
{
   NULL_USE(key);
   
   Array<int> empty(0);
   return (empty);
}

inline
void NullDatabase::getIntegerArray(
   const std::string& key, 
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
bool NullDatabase::isString(const std::string& key)
{
   NULL_USE(key);
   return (true);
}

inline
void NullDatabase::putString(const std::string& key, const std::string& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putStringArray(
   const std::string& key, const Array<std::string>& data)
{
   NULL_USE(key);
   NULL_USE(data);
}

inline
void NullDatabase::putStringArray(
   const std::string& key, const std::string* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
std::string NullDatabase::getString(const std::string& key)
{
   NULL_USE(key); 
   std::string empty;
   return (empty);
}

inline
std::string NullDatabase::getStringWithDefault(
   const std::string& key, 
   const std::string& defaultvalue)
{
   NULL_USE(key);
   NULL_USE(defaultvalue);
   std::string empty;
   return (empty);
}

inline
Array<std::string> NullDatabase::getStringArray(const std::string& key)
{
   NULL_USE(key);
   Array<std::string> empty(0);
   return (empty);
}

inline
void NullDatabase::getStringArray(
   const std::string& key, std::string* data, const int nelements)
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
void NullDatabase::printClassData(std::ostream& os)
{
   NULL_USE(os);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
