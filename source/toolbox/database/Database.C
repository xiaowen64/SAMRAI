//
// File:	Database.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	An abstract base class for the SAMRAI database objects
//

#include "tbox/Database.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Database.I"
#endif

namespace SAMRAI {
   namespace tbox {

Database::~Database()
{
}

/*  
 * Boolean
 */ 

void Database::getScalar(const string& key, bool& scalar)
{
   scalar = getBool(key);
}

void Database::putScalar(const string& key, const bool scalar)
{
   putBool(key, scalar);
}

void Database::getArray(const string& key, Array<bool>& array)
{
   array = getBoolArray(key);
}

void Database::putArray(const string& key, const Array<bool> array)
{
   putBoolArray(key, array);
}

/*  
 * Char
 */ 

void Database::getScalar(const string& key, char& scalar)
{
   scalar = getChar(key);
}

void Database::putScalar(const string& key, const char scalar)
{
   putChar(key, scalar);
}

void Database::getArray(const string& key, Array<char>& array)
{
   array = getCharArray(key);
}

void Database::putArray(const string& key, const Array<char> array)
{
   putCharArray(key, array);
}

/*  
 * Complex
 */ 


void Database::getScalar(const string& key, dcomplex& scalar)
{
   scalar = getComplex(key);
}

void Database::putScalar(const string& key, const dcomplex scalar)
{
   putComplex(key, scalar);
}

void Database::getArray(const string& key, Array<dcomplex>& array)
{
   array = getComplexArray(key);
}

void Database::putArray(const string& key, const Array<dcomplex> array)
{
   putComplexArray(key, array);
}

/*  
 * Float
 */ 


void Database::getScalar(const string& key, float& scalar)
{
   scalar = getFloat(key);
}

void Database::putScalar(const string& key, const float scalar)
{
   putFloat(key, scalar);
}

void Database::getArray(const string& key, Array<float>& array)
{
   array = getFloatArray(key);
}

void Database::putArray(const string& key, const Array<float> array)
{
   putFloatArray(key, array);
}

/*  
 * Double
 */ 

void Database::getScalar(const string& key, double& scalar)
{
   scalar = getDouble(key);
}

void Database::putScalar(const string& key, const double scalar)
{
   putDouble(key, scalar);
}

void Database::getArray(const string& key, Array<double>& array)
{
   array = getDoubleArray(key);
}

void Database::putArray(const string& key, const Array<double> array)
{
   putDoubleArray(key, array);
}

/*  
 * Integer
 */ 

void Database::getScalar(const string& key, int& scalar)
{
   scalar = getInteger(key);
}

void Database::putScalar(const string& key, const int scalar)
{
   putInteger(key, scalar);
}

void Database::getArray(const string& key, Array<int>& array)
{
   array = getIntegerArray(key);
}

void Database::putArray(const string& key, const Array<int> array)
{
   putIntegerArray(key, array);
}


}
}

