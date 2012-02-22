/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for the SAMRAI database objects
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
void Database::getScalar(
   const std::string& key,
   bool& scalar)
{
   scalar = getBool(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putScalar(
   const std::string& key,
   const bool scalar)
{
   putBool(key, scalar);
}

SAMRAI_INLINE_KEYWORD
void Database::getArray(
   const std::string& key,
   Array<bool>& array)
{
   array = getBoolArray(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putArray(
   const std::string& key,
   const Array<bool> array)
{
   putBoolArray(key, array);
}

SAMRAI_INLINE_KEYWORD
void Database::getScalar(
   const std::string& key,
   char& scalar)
{
   scalar = getChar(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putScalar(
   const std::string& key,
   const char scalar)
{
   putChar(key, scalar);
}

SAMRAI_INLINE_KEYWORD
void Database::getArray(
   const std::string& key,
   Array<char>& array)
{
   array = getCharArray(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putArray(
   const std::string& key,
   const Array<char> array)
{
   putCharArray(key, array);
}

SAMRAI_INLINE_KEYWORD
void Database::getScalar(
   const std::string& key,
   dcomplex& scalar)
{
   scalar = getComplex(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putScalar(
   const std::string& key,
   const dcomplex scalar)
{
   putComplex(key, scalar);
}

SAMRAI_INLINE_KEYWORD
void Database::getArray(
   const std::string& key,
   Array<dcomplex>& array)
{
   array = getComplexArray(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putArray(
   const std::string& key,
   const Array<dcomplex> array)
{
   putComplexArray(key, array);
}

SAMRAI_INLINE_KEYWORD
void Database::getScalar(
   const std::string& key,
   float& scalar)
{
   scalar = getFloat(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putScalar(
   const std::string& key,
   const float scalar)
{
   putFloat(key, scalar);
}

SAMRAI_INLINE_KEYWORD
void Database::getArray(
   const std::string& key,
   Array<float>& array)
{
   array = getFloatArray(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putArray(
   const std::string& key,
   const Array<float> array)
{
   putFloatArray(key, array);
}

SAMRAI_INLINE_KEYWORD
void Database::getScalar(
   const std::string& key,
   double& scalar)
{
   scalar = getDouble(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putScalar(
   const std::string& key,
   const double scalar)
{
   putDouble(key, scalar);
}

SAMRAI_INLINE_KEYWORD
void Database::getArray(
   const std::string& key,
   Array<double>& array)
{
   array = getDoubleArray(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putArray(
   const std::string& key,
   const Array<double> array)
{
   putDoubleArray(key, array);
}

SAMRAI_INLINE_KEYWORD
void Database::getScalar(
   const std::string& key,
   int& scalar)
{
   scalar = getInteger(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putScalar(
   const std::string& key,
   const int scalar)
{
   putInteger(key, scalar);
}

SAMRAI_INLINE_KEYWORD
void Database::getArray(
   const std::string& key,
   Array<int>& array)
{
   array = getIntegerArray(key);
}

SAMRAI_INLINE_KEYWORD
void Database::putArray(
   const std::string& key,
   const Array<int> array)
{
   putIntegerArray(key, array);
}

}
}
