/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2017 Lawrence Livermore National Security, LLC
 * Description:   An memory database structure that stores (key,value) pairs in memory
 *
 ************************************************************************/

#include "SAMRAI/tbox/ConduitDatabase.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"

#include <stdlib.h>

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#define MEMORY_DB_ERROR(X) \
   do {                                         \
      pout << "ConduitDatabase: " << X << std::endl << std::flush;       \
      printClassData(pout);                                             \
      pout << "Program abort called..." << std::endl << std::flush;     \
      SAMRAI_MPI::abort();                                              \
   } while (0)

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 * o */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

const int ConduitDatabase::PRINT_DEFAULT = 1;
const int ConduitDatabase::PRINT_INPUT = 2;
const int ConduitDatabase::PRINT_UNUSED = 4;
const int ConduitDatabase::SSTREAM_BUFFER = 4096;

ConduitDatabase::ConduitDatabase(
   const std::string& name):
   d_database_name(name)
{
   d_node = new conduit::Node();
   d_node->reset();
   d_node->set_dtype(conduit::DataType::object());
}

ConduitDatabase::ConduitDatabase(
   conduit::Node* node):
   d_node(node)
{
   d_database_name.clear();
}

ConduitDatabase::ConduitDatabase(
   const std::string& name,
   conduit::Node* node):
   d_database_name(name),
   d_node(node)
{
}

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

ConduitDatabase::~ConduitDatabase()
{
   if (!(d_node->parent())) {
      delete d_node;
   }
}

/*
 *************************************************************************
 *
 * Create memory data file specified by name.
 *
 *************************************************************************
 */

bool
ConduitDatabase::create(
   const std::string& name)
{
   d_database_name = name;
   d_node->reset();

   return true;
}

/*
 *************************************************************************
 *
 * Open memory data file specified by name
 *
 *************************************************************************
 */

bool
ConduitDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   if (read_write_mode == false) {
      TBOX_ERROR("ConduitDatabase::open: ConduitDatabase only supports\n"
         << "read-write mode.  The read_write_mode flag must be true.");

   }
   d_database_name = name;
   d_node->reset();

   return true;
}

/*
 *************************************************************************
 *
 * Close the open data file.
 *
 *************************************************************************
 */

bool
ConduitDatabase::close()
{
   d_database_name = "";
   d_node->reset();

   return true;
}

/*
 *************************************************************************
 *
 * Return whether the key exists in the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::keyExists(
   const std::string& key)
{
   return (d_node->has_child(key));
}

/*
 *************************************************************************
 *
 * Return all of the keys in the database.
 *
 *************************************************************************
 */

std::vector<std::string>
ConduitDatabase::getAllKeys()
{
   return d_node->child_names();
}

/*
 *************************************************************************
 *
 * Get the type of the array entry associated with the specified key
 *
 *************************************************************************
 */
enum Database::DataType
ConduitDatabase::getArrayType(
   const std::string& key)
{
   findChildNodeOrExit(key);

   if (isBool(key)) {
      return SAMRAI_BOOL;
   } else if (isDatabase(key)) {
      return SAMRAI_DATABASE;
   } else if (isDatabaseBox(key)) {
      return SAMRAI_BOX;
   } else if (isDouble(key)) {
      return SAMRAI_DOUBLE;
   } else if (isFloat(key)) {
      return SAMRAI_FLOAT;
   } else if (isInteger(key)) {
      return SAMRAI_INT;
   } else if (isString(key)) {
      return SAMRAI_STRING;
   } else {
      return SAMRAI_INVALID;
   } 
}

/*
 *************************************************************************
 *
 * Get the size of the array entry associated with the specified key;
 * return 0 if the key does not exist.
 *
 *************************************************************************
 */

size_t
ConduitDatabase::getArraySize(
   const std::string& key)
{
   if (!(*d_node).has_child(key) || isDatabase(key)) {
      return 0;
   } else {
      if (isBool(key)) {
         return (*d_node)[key]["data"].dtype().number_of_elements();
      } else if (isDatabaseBox(key)) {
         return (*d_node)[key]["dimension"].dtype().number_of_elements();
      } if (isString(key)) {
         return (*d_node)[key].number_of_children();
      } else {
         return (*d_node)[key].dtype().number_of_elements();
      }
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage the database values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isDatabase(
   const std::string& key)
{
   bool is_database = false;
   if (d_node->has_child(key)) {
      if ((*d_node)[key].dtype().is_object() && !isBool(key) &&
          !isDatabaseBox(key) && !isString(key)) {
         is_database = true;
      }
   }

   return is_database;
}

std::shared_ptr<Database>
ConduitDatabase::putDatabase(
   const std::string& key)
{
   deleteKeyIfFound(key);

   (*d_node)[key].set_dtype(conduit::DataType::object());
   std::shared_ptr<Database> database(new ConduitDatabase(key, &((*d_node)[key])));
   d_types[key] = SAMRAI_DATABASE;

   return database;
}

std::shared_ptr<Database>
ConduitDatabase::getDatabase(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isDatabase(key)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a database...");
   }
   std::shared_ptr<Database> database(new ConduitDatabase(key, &((*d_node)[key])));
   return database;
}

/*
 *************************************************************************
 *
 * Member functions that manage boolean values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isBool(
   const std::string& key)
{
   bool is_bool = false;
   if (d_node->has_child(key) && (*d_node)[key].has_child("data") &&
       (*d_node)[key].has_child("bool") &&
       (*d_node)[key]["data"].dtype().is_uint8() &&
       (*d_node)[key]["bool"].as_string() == "true") {
      is_bool = true;
   }
   return is_bool;
}

void
ConduitDatabase::putBool(
   const std::string& key,
   const bool& data)
{
   putBoolArray(key, &data, 1);
}

void
ConduitDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   std::vector<conduit::uint8> uint8_vec(nelements, 0);
   for (int i = 0; i < nelements; ++i) {
      if (data[i]) {
         uint8_vec[i] = 1;
      }
   }
   (*d_node)[key]["data"].set(&(uint8_vec[0]), nelements);
   (*d_node)[key]["bool"] = "true";
   d_types[key] = SAMRAI_BOOL;
}

bool
ConduitDatabase::getBool(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isBool(key) ||
       (*d_node)[key]["data"].dtype().number_of_elements() != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a boolean scalar...");
   }
   return static_cast<bool>((*d_node)[key]["data"].as_uint8());
}

bool
ConduitDatabase::getBoolWithDefault(
   const std::string& key,
   const bool& defaultvalue)
{
   if (d_node->has_child(key)) return getBool(key);

   putBool(key, defaultvalue);
   return defaultvalue;
}

std::vector<bool>
ConduitDatabase::getBoolVector(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isBool(key)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a boolean...");
   }
   conduit::uint8_array int_vals = (*d_node)[key]["data"].as_uint8_array();
   unsigned int vec_size = (*d_node)[key].dtype().number_of_elements();
   std::vector<bool> bool_vec(vec_size, false);
   for (unsigned int i = 0; i < vec_size; ++i) {
      if (int_vals[i] != 0) {
         bool_vec[i] = true;
      }
   }
   return bool_vec;
}

void
ConduitDatabase::getBoolArray(
   const std::string& key,
   bool* data,
   const size_t nelements)
{
   std::vector<bool> tmp = getBoolVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage box values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isDatabaseBox(
   const std::string& key)
{
   bool is_box = false;
   if (d_node->has_child(key) && (*d_node)[key].has_child("box") &&
       (*d_node)[key].has_child("dimension") &&
       (*d_node)[key].has_child("lower") &&
       (*d_node)[key].has_child("upper") &&
       (*d_node)[key]["dimension"].dtype().is_uint8() &&
       (*d_node)[key]["lower"].dtype().is_int() &&
       (*d_node)[key]["upper"].dtype().is_int() &&
       (*d_node)[key]["box"].as_string() == "true") {
      is_box = true;
   }
   return is_box;
}

void
ConduitDatabase::putDatabaseBox(
   const std::string& key,
   const DatabaseBox& data)
{
   putDatabaseBoxArray(key, &data, 1);
}

void
ConduitDatabase::putDatabaseBoxVector(
   const std::string& key,
   const std::vector<DatabaseBox>& data)
{
   putDatabaseBoxArray(key, &data[0], data.size());
}

void
ConduitDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   std::vector<conduit::uint8> dim_vec(nelements);
   std::vector<int> lo_vec(nelements * SAMRAI::MAX_DIM_VAL);
   std::vector<int> hi_vec(nelements * SAMRAI::MAX_DIM_VAL);
   for (int i = 0; i < nelements; ++i) {
      dim_vec[i] = data[i].getDimVal();

      for (int d = 0; d < dim_vec[i]; ++d) {
         lo_vec[i*SAMRAI::MAX_DIM_VAL + d] = data[i].lower(d);
         hi_vec[i*SAMRAI::MAX_DIM_VAL + d] = data[i].upper(d);
      }
      for (int d = dim_vec[i]; d < SAMRAI::MAX_DIM_VAL; ++d) {
         lo_vec[i*SAMRAI::MAX_DIM_VAL + d] = 0;
         hi_vec[i*SAMRAI::MAX_DIM_VAL + d] = 0;
      }
   }

   (*d_node)[key]["box"] = "true";
   (*d_node)[key]["dimension"].set(&(dim_vec[0]), nelements);
   (*d_node)[key]["lower"].set(&(lo_vec[0]), nelements*SAMRAI::MAX_DIM_VAL);
   (*d_node)[key]["upper"].set(&(hi_vec[0]), nelements*SAMRAI::MAX_DIM_VAL);
   d_types[key] = SAMRAI_BOX;
}

DatabaseBox
ConduitDatabase::getDatabaseBox(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isDatabaseBox(key) ||
       (*d_node)[key]["dimension"].dtype().number_of_elements() != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single box...");
   }

   tbox::Dimension dim((*d_node)[key]["dimension"].as_uint8());
   int* lower = (*d_node)[key]["lower"].as_int_ptr(); 
   int* upper = (*d_node)[key]["upper"].as_int_ptr(); 

   DatabaseBox db_box(dim, lower, upper);

   return db_box;
}

DatabaseBox
ConduitDatabase::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   if (d_node->has_child(key)) return getDatabaseBox(key);

   putDatabaseBox(key, defaultvalue);
   return defaultvalue;
}

std::vector<DatabaseBox>
ConduitDatabase::getDatabaseBoxVector(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isDatabaseBox(key)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a DatabaseBox...");
   }
   std::vector<DatabaseBox> box_vec;
   size_t vec_size = (*d_node)[key]["dimension"].dtype().number_of_elements();
   conduit::uint8_array dim_vals = (*d_node)[key]["dimension"].as_uint8_array();
   conduit::int_array lo_vals = (*d_node)[key]["lower"].as_int_array();
   conduit::int_array hi_vals = (*d_node)[key]["upper"].as_int_array();
   for (size_t i = 0; i < vec_size; ++i) {
      tbox::Dimension dim(dim_vals[i]);
      int * lower = &(lo_vals[i*SAMRAI::MAX_DIM_VAL]);
      int * upper = &(hi_vals[i*SAMRAI::MAX_DIM_VAL]);
      DatabaseBox db_box(dim, lower, upper);
      box_vec.push_back(db_box); 
   }
   return box_vec;
}

void
ConduitDatabase::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const size_t nelements)
{
   std::vector<DatabaseBox> tmp = getDatabaseBoxVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage character values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isChar(
   const std::string& key)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
   return false;
}

void
ConduitDatabase::putChar(
   const std::string& key,
   const char& data)
{
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");

   putCharArray(key, &data, 1);
}

void
ConduitDatabase::putCharVector(
   const std::string& key,
   const std::vector<char>& data)
{
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
   putCharArray(key, &data[0], data.size());
}

void
ConduitDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
}

char
ConduitDatabase::getChar(
   const std::string& key)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
   return 0;
}

char
ConduitDatabase::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
   return defaultvalue;
}

std::vector<char>
ConduitDatabase::getCharVector(
   const std::string& key)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
   std::vector<char> cvec;
   return cvec;
}

void
ConduitDatabase::getCharArray(
   const std::string& key,
   char* data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
   MEMORY_DB_ERROR("no char implemented in ConduitDatabase");
}

/*
 *************************************************************************
 *
 * Member functions that manage complex values within the database.
 * Note that complex numbers may be promoted from integers, floats,
 * and doubles.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isComplex(
   const std::string& key)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
   return false;
}

void
ConduitDatabase::putComplex(
   const std::string& key,
   const dcomplex& data)
{
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
   putComplexArray(key, &data, 1);
}

void
ConduitDatabase::putComplexVector(
   const std::string& key,
   const std::vector<dcomplex>& data)
{
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
   putComplexArray(key, &data[0], data.size());
}

void
ConduitDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
}

dcomplex
ConduitDatabase::getComplex(
   const std::string& key)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
   dcomplex value(0.0, 0.0);
   return value;
}

dcomplex
ConduitDatabase::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
   return defaultvalue;
}

std::vector<dcomplex>
ConduitDatabase::getComplexVector(
   const std::string& key)
{
   NULL_USE(key);
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
   std::vector<dcomplex> array;
   return array;
}

void
ConduitDatabase::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
   MEMORY_DB_ERROR("no complex implemented in ConduitDatabase");
}

/*
 *************************************************************************
 *
 * Member functions that manage double values within the database.
 * Note that doubles may be promoted from integers or floats.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isDouble(
   const std::string& key)
{
   bool is_double = false;
   if (d_node->has_child(key) && (*d_node)[key].dtype().is_double()) {
      is_double = true;
   }
   return is_double;
}

void
ConduitDatabase::putDouble(
   const std::string& key,
   const double& data)
{
   putDoubleArray(key, &data, 1);
}

void
ConduitDatabase::putDoubleVector(
   const std::string& key,
   const std::vector<double>& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set(data);
   d_types[key] = SAMRAI_DOUBLE;
}

void
ConduitDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const size_t nelements)
{
   std::vector<double> dbl_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      dbl_vec[i] = data[i];
   }
   putDoubleVector(key, dbl_vec);
}

double
ConduitDatabase::getDouble(
   const std::string& key)
{
   double value = 0.0;
   findChildNodeOrExit(key);

   if ((*d_node)[key].dtype().number_of_elements() != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   const conduit::DataType& datatype = (*d_node)[key].dtype();
   if (datatype.is_double()) {
      value = (*d_node)[key].as_double();
   } else if (datatype.is_integer()) {
      value = static_cast<double>((*d_node)[key].as_int());
   } else if (datatype.is_float()) {
      value = static_cast<double>((*d_node)[key].as_float());
   } else {
      MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   return value;
}

double
ConduitDatabase::getDoubleWithDefault(
   const std::string& key,
   const double& defaultvalue)
{
   if (d_node->has_child(key)) return getDouble(key);

   putDouble(key, defaultvalue);
   return defaultvalue;
}

std::vector<double>
ConduitDatabase::getDoubleVector(
   const std::string& key)
{
   findChildNodeOrExit(key);
   std::vector<double> dbl_vec;

   const conduit::DataType& datatype = (*d_node)[key].dtype();
   size_t vec_size = datatype.number_of_elements();
   dbl_vec.resize(vec_size);
   if (datatype.is_double()) {
      conduit::double_array dbl_array = (*d_node)[key].as_double_array();
      for (size_t i = 0; i < vec_size; ++i) {
         dbl_vec[i] = dbl_array[i];
      }
   } else if (datatype.is_integer()) {
      conduit::int_array int_vals = (*d_node)[key].as_int_array();
      for (size_t i = 0; i < vec_size; ++i) {
         dbl_vec[i] = static_cast<double>(int_vals[i]);
      }
   } else if (datatype.is_float()) {
      conduit::float_array float_array = (*d_node)[key].as_float_array();
      for (size_t i = 0; i < vec_size; ++i) {
         dbl_vec[i] = static_cast<double>(float_array[i]);
      }
   } else {
      MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   return dbl_vec;
}

void
ConduitDatabase::getDoubleArray(
   const std::string& key,
   double* data,
   const size_t nelements)
{
   std::vector<double> tmp = getDoubleVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage float values within the database.
 * Note that floats may be promoted from integers or truncated from
 * doubles (without a warning).
 *
 *************************************************************************
 */

bool
ConduitDatabase::isFloat(
   const std::string& key)
{
   bool is_float = false;
   if (d_node->has_child(key) && (*d_node)[key].dtype().is_float()) {
      is_float = true;
   }
   return is_float;
}

void
ConduitDatabase::putFloat(
   const std::string& key,
   const float& data)
{
   putFloatArray(key, &data, 1);
}

void
ConduitDatabase::putFloatVector(
   const std::string& key,
   const std::vector<float>& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set(data);
   d_types[key] = SAMRAI_FLOAT;
}

void
ConduitDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const size_t nelements)
{
   std::vector<float> flt_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      flt_vec[i] = data[i];
   }

   putFloatVector(key, flt_vec);
}

float
ConduitDatabase::getFloat(
   const std::string& key)
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   float value = 0.0;
   findChildNodeOrExit(key);

   if ((*d_node)[key].dtype().number_of_elements() != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   const conduit::DataType& datatype = (*d_node)[key].dtype();
   if (datatype.is_float()) {
      value = (*d_node)[key].as_float();
   } else if (datatype.is_integer()) {
      value = static_cast<float>((*d_node)[key].as_int());
   } else if (datatype.is_double()) {
      value = static_cast<float>((*d_node)[key].as_double());
   } else {
      MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   return value;
}

float
ConduitDatabase::getFloatWithDefault(
   const std::string& key,
   const float& defaultvalue)
{
   if (d_node->has_child(key)) return getFloat(key);

   putFloat(key, defaultvalue);
   return defaultvalue;
}

std::vector<float>
ConduitDatabase::getFloatVector(
   const std::string& key)
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   findChildNodeOrExit(key);
   std::vector<float> flt_vec;

   const conduit::DataType& datatype = (*d_node)[key].dtype();
   size_t vec_size = datatype.number_of_elements();
   flt_vec.resize(vec_size);
   if (datatype.is_float()) {
      conduit::float_array flt_array = (*d_node)[key].as_float_array();
      for (size_t i = 0; i < vec_size; ++i) {
         flt_vec[i] = flt_array[i];
      }
   } else if (datatype.is_integer()) {
      conduit::int_array int_vals = (*d_node)[key].as_int_array();
      for (size_t i = 0; i < vec_size; ++i) {
         flt_vec[i] = static_cast<float>(int_vals[i]);
      }
   } else if (datatype.is_double()) {
      conduit::double_array dbl_array = (*d_node)[key].as_double_array();
      for (size_t i = 0; i < vec_size; ++i) {
         flt_vec[i] = static_cast<float>(dbl_array[i]);
      }
   } else {
      MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   return flt_vec;

}

void
ConduitDatabase::getFloatArray(
   const std::string& key,
   float* data,
   const size_t nelements)
{
   std::vector<float> tmp = getFloatVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage integer values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isInteger(
   const std::string& key)
{
   bool is_int = false;
   if (d_node->has_child(key) && (*d_node)[key].dtype().is_int()) {
      is_int = true;
   }
   return is_int;
}

void
ConduitDatabase::putInteger(
   const std::string& key,
   const int& data)
{
   putIntegerArray(key, &data, 1);
}

void
ConduitDatabase::putIntegerVector(
   const std::string& key,
   const std::vector<int>& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set(data);
   d_types[key] = SAMRAI_INT;
}

void
ConduitDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const size_t nelements)
{
   std::vector<int> int_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      int_vec[i] = data[i];
   }
   putIntegerVector(key, int_vec);
}

int
ConduitDatabase::getInteger(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!(*d_node)[key].dtype().is_int() ||
       (*d_node)[key].dtype().number_of_elements() != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not an integer scalar...");
   }
   return (*d_node)[key].as_int();
}

int
ConduitDatabase::getIntegerWithDefault(
   const std::string& key,
   const int& defaultvalue)
{
   if (d_node->has_child(key)) return getInteger(key);

   putInteger(key, defaultvalue);
   return defaultvalue;
}

std::vector<int>
ConduitDatabase::getIntegerVector(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!(*d_node)[key].dtype().is_int()) {
      MEMORY_DB_ERROR("Key=" << key << " is not an integer...");
   }
   size_t vec_size = (*d_node)[key].dtype().number_of_elements();
   std::vector<int> int_vec(vec_size);
   conduit::int_array int_vals = (*d_node)[key].as_int_array();
   for (size_t i = 0; i < vec_size; ++i) {
      int_vec[i] = int_vals[i];
   }
   return int_vec;
}

void
ConduitDatabase::getIntegerArray(
   const std::string& key,
   int* data,
   const size_t nelements)
{
   std::vector<int> tmp = getIntegerVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage string values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isString(
   const std::string& key)
{
   bool is_string = false;
   if (d_node->has_child(key) &&
       (*d_node)[key].has_child("str0") &&
       (*d_node)[key]["str0"].dtype().is_string()) {
      is_string = true;
   }
   return is_string;
}

void
ConduitDatabase::putString(
   const std::string& key,
   const std::string& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key]["str0"].set_string(data);
   d_types[key] = SAMRAI_STRING;
}

void
ConduitDatabase::putStringVector(
   const std::string& key,
   const std::vector<std::string>& data)
{
   deleteKeyIfFound(key);
   int i = 0;
   for (std::vector<std::string>::const_iterator itr = data.begin();
        itr != data.end(); ++itr) {
      std::stringstream ss;
      ss << i;
      std::string id = "str" + ss.str();
      (*d_node)[key][id].set_string(*itr);
      ++i; 
   }
   d_types[key] = SAMRAI_STRING;
}

void
ConduitDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const size_t nelements)
{
   std::vector<std::string> str_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      str_vec[i] = data[i];
   }
   putStringVector(key, str_vec);
}

std::string
ConduitDatabase::getString(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isString(key) ||
       (*d_node)[key].dtype().number_of_elements() != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single string ...");
   }
   return (*d_node)[key]["str0"].as_string();
}

std::string
ConduitDatabase::getStringWithDefault(
   const std::string& key,
   const std::string& defaultvalue)
{
   if (d_node->has_child(key)) return getString(key);

   putString(key, defaultvalue);
   return defaultvalue;
}

std::vector<std::string>
ConduitDatabase::getStringVector(
   const std::string& key)
{
   findChildNodeOrExit(key);
   if (!isString(key)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a string...");
   }

   size_t nelements = (*d_node)[key].number_of_children();
   std::vector<std::string> str_vec;

   for (size_t i = 0; i < nelements; ++i) {
      std::stringstream ss;
      ss << i;
      std::string id = "str" + ss.str();
      str_vec.push_back((*d_node)[key][id].as_string());
   } 

   return str_vec;
}

void
ConduitDatabase::getStringArray(
   const std::string& key,
   std::string* data,
   const size_t nelements)
{
   std::vector<std::string> tmp = getStringVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

std::string
ConduitDatabase::getName()
{
   return d_database_name;
}

std::string
ConduitDatabase::getName() const
{
   return d_database_name;
}

/*
 *************************************************************************
 *
 * Search the current database for a matching key.  If found, delete
 * that key and return true.  If the key does not exist, then return
 * false.
 *
 *************************************************************************
 */

bool ConduitDatabase::deleteKeyIfFound(
   const std::string& key)
{
   if (d_node->has_child(key)) {
      d_node->remove(key);
      d_types.erase(key);
      return true;
   } else {
      return false;
   }
}

/*
 *************************************************************************
 *
 * Find the child node associated with the specified key, exit with error
 * if not found.
 *
 *************************************************************************
 */

void
ConduitDatabase::findChildNodeOrExit(
   const std::string& key)
{
   if (!d_node->has_child(key)) {
      MEMORY_DB_ERROR("Key ``" << key << "'' does not exist in the database...");
   }
}

void
ConduitDatabase::copyDatabase(const std::shared_ptr<Database>& database)
{
   std::vector<std::string> keys(database->getAllKeys());

   for (std::vector<std::string>::const_iterator k_itr = keys.begin();
        k_itr != keys.end(); ++k_itr) {

      const std::string& key = *k_itr;
      Database::DataType my_type = database->getArrayType(key);

      size_t size =  database->getArraySize(key);
      if (my_type == SAMRAI_DATABASE) {
         std::shared_ptr<Database> child_db = database->getDatabase(key);
         std::shared_ptr<ConduitDatabase> new_db =
            SAMRAI_SHARED_PTR_CAST<ConduitDatabase, Database>(
               putDatabase(key));
         new_db->copyDatabase(child_db); 
      } else if (my_type == SAMRAI_BOOL) {
         if (size == 1) {
            putBool(key, database->getBool(key));
         } else if (size % 2 == 0) {
            bool barray[size];
            database->getBoolArray(key, barray, size);
            putBoolArray(key, barray, size);
         } else {
            std::vector<bool> bvec(database->getBoolVector(key));
            putBoolVector(key, bvec);
         }
      } else if (my_type == SAMRAI_CHAR) {
         if (size == 1) {
            putChar(key, database->getChar(key));
         } else if (size % 2 == 0) {
            char barray[size];
            database->getCharArray(key, barray, size);
            putCharArray(key, barray, size);
         } else {
            std::vector<char> bvec(database->getCharVector(key));
            putCharVector(key, bvec);
         }
      } else if (my_type == SAMRAI_INT) {
         if (size == 1) {
            putInteger(key, database->getInteger(key));
         } else if (size % 2 == 0) {
            int barray[size];
            database->getIntegerArray(key, barray, size);
            putIntegerArray(key, barray, size);
         } else {
            std::vector<int> bvec(database->getIntegerVector(key));
            putIntegerVector(key, bvec);
         }
      } else if (my_type == SAMRAI_COMPLEX) {
         if (size == 1) {
            putComplex(key, database->getComplex(key));
         } else if (size % 2 == 0) {
            dcomplex barray[size];
            database->getComplexArray(key, barray, size);
            putComplexArray(key, barray, size);
         } else {
            std::vector<dcomplex> bvec(database->getComplexVector(key));
            putComplexVector(key, bvec);
         }
      } else if (my_type == SAMRAI_DOUBLE) {
         if (size == 1) {
            putDouble(key, database->getDouble(key));
         } else if (size % 2 == 0) {
            double barray[size];
            database->getDoubleArray(key, barray, size);
            putDoubleArray(key, barray, size);
         } else {
            std::vector<double> bvec(database->getDoubleVector(key));
            putDoubleVector(key, bvec);
         }
      } else if (my_type == SAMRAI_FLOAT) {
         if (size == 1) {
            putFloat(key, database->getFloat(key));
         } else if (size % 2 == 0) {
            float barray[size];
            database->getFloatArray(key, barray, size);
            putFloatArray(key, barray, size);
         } else {
            std::vector<float> bvec(database->getFloatVector(key));
            putFloatVector(key, bvec);
         }
      } else if (my_type == SAMRAI_STRING) {
         if (size == 1) {
            putString(key, database->getString(key));
         } else if (size % 2 == 0) {
            std::string barray[size];
            database->getStringArray(key, barray, size);
            putStringArray(key, barray, size);
         } else {
            std::vector<std::string> bvec(database->getStringVector(key));
            putStringVector(key, bvec);
         }
      } else if (my_type == SAMRAI_BOX) {
         if (size == 1) {
            putDatabaseBox(key, database->getDatabaseBox(key));
         } else if (size % 2 == 0) {
            DatabaseBox barray[size];
            database->getDatabaseBoxArray(key, barray, size);
            putDatabaseBoxArray(key, barray, size);
         } else {
            std::vector<DatabaseBox> bvec(database->getDatabaseBoxVector(key));
            putDatabaseBoxVector(key, bvec);
         }
      }

/*
                   SAMRAI_DATABASE,
                   SAMRAI_BOOL,
                   SAMRAI_CHAR,
                   SAMRAI_INT,
                   SAMRAI_COMPLEX,
                   SAMRAI_DOUBLE,
                   SAMRAI_FLOAT,
                   SAMRAI_STRING,
                   SAMRAI_BOX };
*/

   }
}


/*
 *************************************************************************
 *
 * Print the entire database to the specified output stream.
 *
 *************************************************************************
 */

void
ConduitDatabase::printClassData(
   std::ostream& os)
{
   printDatabase(os, 0, PRINT_DEFAULT | PRINT_INPUT | PRINT_UNUSED);
}

/*
 *************************************************************************
 *
 * Print database data to the specified output stream.
 *
 *************************************************************************
 */

void
ConduitDatabase::printDatabase(
   std::ostream& os,
   const int indent,
   const int toprint) const
{
// need to implement print
#if 0
   /*
    * Get the maximum key width in the output (excluding databases)
    */

   int width = 0;
   for (std::list<KeyData>::const_iterator k = d_keyvalues.begin();
        k != d_keyvalues.end(); ++k) {
      if (((k->d_from_default) && (toprint & PRINT_DEFAULT))
          || ((k->d_accessed) && (toprint & PRINT_INPUT))
          || (!(k->d_accessed) && (toprint & PRINT_UNUSED))) {
         if (k->d_type != Database::SAMRAI_DATABASE) {
            const int keywidth = static_cast<int>(k->d_key.length());
            if (keywidth > width) {
               width = keywidth;
            }
         }
      }
   }

   /*
    * Iterate over all non-database keys in the database and output key values
    */

   indentStream(os, indent);
   os << d_database_name << " {\n";
   for (std::list<KeyData>::const_iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); ++i) {

      if (((i->d_from_default) && (toprint & PRINT_DEFAULT))
          || ((i->d_accessed) && (toprint & PRINT_INPUT))
          || (!(i->d_accessed) && (toprint & PRINT_UNUSED))) {

#ifndef LACKS_SSTREAM
         std::ostringstream sstream;
#else
         char sstream_buffer[SSTREAM_BUFFER];
         std::ostrstream sstream(sstream_buffer, SSTREAM_BUFFER);
#endif

         switch (i->d_type) {

            case Database::SAMRAI_INVALID: {
               break;
            }

            case Database::SAMRAI_DATABASE: {
               break;
            }

            case Database::SAMRAI_BOOL: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<bool>::size_type n = i->d_boolean.size();
               for (std::vector<bool>::size_type j = 0; j < n; ++j) {
                  sstream << (i->d_boolean[j] ? "TRUE" : "FALSE");
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_BOX: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<DatabaseBox>::size_type n = i->d_box.size();
               for (std::vector<DatabaseBox>::size_type j = 0; j < n; ++j) {
                  const int m = i->d_box[j].getDimVal();
                  sstream << "[(";
                  for (int k = 0; k < m; ++k) {
                     sstream << i->d_box[j].lower(k);
                     if (k < m - 1) {
                        sstream << ",";
                     }
                  }
                  sstream << "),(";
                  for (int l = 0; l < m; ++l) {
                     sstream << i->d_box[j].upper(l);
                     if (l < m - 1) {
                        sstream << ",";
                     }
                  }
                  sstream << ")]";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_CHAR: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<char>::size_type n = i->d_char.size();
               for (std::vector<char>::size_type j = 0; j < n; ++j) {
                  sstream << "'" << i->d_char[j] << "'";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_COMPLEX: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<dcomplex>::size_type n = i->d_complex.size();
               for (std::vector<dcomplex>::size_type j = 0; j < n; ++j) {
                  sstream << i->d_complex[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_DOUBLE: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<double>::size_type n = i->d_double.size();
               for (std::vector<double>::size_type j = 0; j < n; ++j) {
                  sstream << i->d_double[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_FLOAT: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<float>::size_type n = i->d_float.size();
               for (std::vector<float>::size_type j = 0; j < n; ++j) {
                  sstream << i->d_float[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_INT: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<int>::size_type n = i->d_integer.size();
               for (std::vector<int>::size_type j = 0; j < n; ++j) {
                  sstream << i->d_integer[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_STRING: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - static_cast<int>(i->d_key.length()));
               sstream << " = ";
               const std::vector<std::string>::size_type n = i->d_string.size();
               for (std::vector<std::string>::size_type j = 0; j < n; ++j) {
                  sstream << "\"" << i->d_string[j] << "\"";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }
         }

         /*
          * Output whether the key was used or default in column 60
          */

         if (i->d_type != Database::SAMRAI_DATABASE) {
#ifndef LACKS_SSTREAM
            const int tab = static_cast<int>(59 - sstream.str().length());
#else
            const int tab = static_cast<int>(59 - sstream.pcount());
#endif
            if (tab > 0) {
               indentStream(sstream, tab);
            }
            if (i->d_from_default) {
               sstream << " // from default";
            } else if (i->d_accessed) {
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

   for (std::list<KeyData>::const_iterator j = d_keyvalues.begin();
        j != d_keyvalues.end(); ++j) {
      if (j->d_type == Database::SAMRAI_DATABASE) {
         std::shared_ptr<ConduitDatabase> db(
            SAMRAI_SHARED_PTR_CAST<ConduitDatabase, Database>(j->d_database));
         db->printDatabase(os, indent + 3, toprint);
      }
   }

   indentStream(os, indent);
   os << "}\n";
#endif
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
