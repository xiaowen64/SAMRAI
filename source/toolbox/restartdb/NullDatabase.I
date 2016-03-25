//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/restartdb/NullDatabase.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2620 $
// Modified:	$LastChangedDate: 2008-11-19 14:24:28 -0800 (Wed, 19 Nov 2008) $
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
*                                                                       *
* Create memory data file specified by name.                            *
*                                                                       *
*************************************************************************
*/

inline
bool NullDatabase::create(const std::string& name) 
{
   NULL_USE(name);
   return true;
}


/*
*************************************************************************
*                                                                       *
* Open memory data file specified by name                               *
*                                                                       *
*************************************************************************
*/

inline
bool NullDatabase::open(const std::string& name) 
{
   NULL_USE(name);
   return true;
}

/*
*************************************************************************
*                                                                       *
* Close the open data file.                                             *
*                                                                       *
*************************************************************************
*/

inline
bool NullDatabase::close() 
{
   return true;
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
* Always returns INVALID.	 				        *  
*									*
*************************************************************************
*/

inline
Database::DataType NullDatabase::getArrayType(const std::string& key)
{
   NULL_USE(key);
   return Database::SAMRAI_INVALID;
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
void NullDatabase::putBoolArray(
   const std::string& key, const bool* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<bool> NullDatabase::getBoolArray(const std::string& key)
{
   NULL_USE(key);
   Array<bool> empty(0);
   return (empty);
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
void NullDatabase::putDatabaseBoxArray(
   const std::string& key, const DatabaseBox* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}


inline
Array<DatabaseBox> NullDatabase::getDatabaseBoxArray(const std::string& key)
{
   NULL_USE(key);

   Array<DatabaseBox> empty(0);
   return (empty);
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
void NullDatabase::putCharArray(
   const std::string& key, const char* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<char> NullDatabase::getCharArray(const std::string& key)
{
   NULL_USE(key);

   Array<char> empty(0);
   return (empty);
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
void NullDatabase::putComplexArray(
   const std::string& key, const dcomplex* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<dcomplex> NullDatabase::getComplexArray(const std::string& key)
{
   NULL_USE(key);
  
   Array<dcomplex> empty(0);
   return(empty);
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
void NullDatabase::putDoubleArray(
   const std::string& key, const double* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<double> NullDatabase::getDoubleArray(const std::string& key)
{
   NULL_USE(key);
   Array<double> empty(0);
   return (empty);
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
void NullDatabase::putFloatArray(
   const std::string& key, const float* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<float> NullDatabase::getFloatArray(const std::string& key)
{
   NULL_USE(key); 

   Array<float> empty(0);
   return(empty);
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
void NullDatabase::putIntegerArray(
   const std::string& key, const int* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<int> NullDatabase::getIntegerArray(const std::string& key)
{
   NULL_USE(key);
   
   Array<int> empty(0);
   return (empty);
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
void NullDatabase::putStringArray(
   const std::string& key, const std::string* const data, const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

inline
Array<std::string> NullDatabase::getStringArray(const std::string& key)
{
   NULL_USE(key);
   Array<std::string> empty(0);
   return (empty);
}

inline
std::string NullDatabase::getName(void)
{
   return std::string();
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
