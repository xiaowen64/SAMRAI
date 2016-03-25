//
// File:        main-HDF5.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Tests HDF database in SAMRAI
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/DatabaseBox.h"
#include "tbox/Complex.h"
#include "tbox/HDFDatabase.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include <string>

using namespace std;
using namespace SAMRAI;


// Number of (non-abortive) failures.
int number_of_failures = 0;

float db_float_val = 3.14159;
int db_int_val = 4;

dcomplex arraydb_compArray0 = dcomplex(1,2);
dcomplex arraydb_compArray1 = dcomplex(2,3);
dcomplex arraydb_compArray2 = dcomplex(3,4);
bool arraydb_boolArray0 = true;
bool arraydb_boolArray1 = false;
bool arraydb_boolArray2 = false;
int arraydb_intArray0 = 0;
int arraydb_intArray1 = 1;
int arraydb_intArray2 = 2;
int arraydb_intArray3 = 3;
int arraydb_intArray4 = 4;
string arraydb_stringArray0 = "This is 1 test.";
string arraydb_stringArray1 = "This is 2nd test.";
string arraydb_stringArray2 = "This is a long 3rd test.";
float arraydb_floatArray0 = 0*1.2;
float arraydb_floatArray1 = 1*1.2;
float arraydb_floatArray2 = 2*1.2;
float arraydb_floatArray3 = 3*1.2;
float arraydb_floatArray4 = 4*1.2;
double arraydb_doubleArray0 = 0*1.111111;
double arraydb_doubleArray1 = 1*1.111111;
double arraydb_doubleArray2 = 2*1.111111;
double arraydb_doubleArray3 = 3*1.111111;
double arraydb_doubleArray4 = 4*1.111111;
double arraydb_doubleArray5 = 5*1.111111;
char arraydb_charArray0 = 'a';
char arraydb_charArray1 = 'b';
tbox::DatabaseBox arraydb_boxArray0; 
tbox::DatabaseBox arraydb_boxArray1; 
tbox::DatabaseBox arraydb_boxArray2; 

float scalardb_float1 = 1.111;
float scalardb_float2 = 2.222;
float scalardb_float3 = 3.333;
double scalardb_full_thisDouble = 123.456;
dcomplex scalardb_full_thisComplex = dcomplex(2.3, 4.5);
int scalardb_full_thisInt = 89;
float scalardb_full_thisFloat = 9.9;
bool scalardb_full_thisBool = true;
string scalardb_full_thisString = "This is a test.";
char scalardb_full_thisChar = 'q';
int ilo[3] = {0,0,0};
int ihi[3] = {1,1,1};
tbox::DatabaseBox scalardb_full_thisBox(3, ilo, ihi);


void writeHDFTestData(tbox::Pointer<tbox::Database> db);
void readHDFTestData(tbox::Pointer<tbox::Database> db);
void testDatabaseContents(tbox::Pointer<tbox::Database> db,
                          const string& tag);

class RestartTester : public tbox::Serializable 
{
public:   

   RestartTester() 
   {
      tbox::RestartManager::getManager()->registerRestartItem("RestartTester",
                                                             this);
   }

   virtual ~RestartTester() {}

   void putToDatabase(tbox::Pointer<tbox::Database> db) 
   {
      writeHDFTestData(db);
   }

   void getFromDatabase() 
   {
      tbox::Pointer<tbox::Database> root_db =
         tbox::RestartManager::getManager()->getRootDatabase();

      tbox::Pointer<tbox::Database> db;
      if (root_db->isDatabase("RestartTester")) {
         db = root_db->getDatabase("RestartTester");
      } 

      readHDFTestData(db);
   }

};

int main(int argc, char *argv[]) 
{
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   tbox::PIO::logAllNodes("HDF5test.log");
// tbox::PIO::logOnlyNodeZero("HDF5test.log");

#ifdef HAVE_HDF5

   tbox::pout << "\n--- HDF5 database tests BEGIN ---" << endl;

   tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

   RestartTester hdf_tester;

   tbox::pout << "\n--- HDF5 write database tests BEGIN ---" << endl;

   arraydb_boxArray0.setDimension(3);
   arraydb_boxArray1.setDimension(2);
   arraydb_boxArray2.setDimension(1);

   restart_manager->writeRestartFile("test_dir", 0);

   tbox::pout << "\n--- HDF5 write database tests END ---" << endl;

   tbox::pout << "\n--- HDF5 read database tests BEGIN ---" <<  endl;

   restart_manager->openRestartFile("test_dir", 0, tbox::MPI::getNodes());

   hdf_tester.getFromDatabase();

   restart_manager->closeRestartFile();

   tbox::pout << "\n--- HDF5 read database tests END ---" << endl;

   tbox::pout << "\n--- HDF5 database tests END ---" << endl;

#endif

   if (number_of_failures == 0) {
      tbox::pout << "\nPASSED:  HDF5" << endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();

   return number_of_failures;

}

/*
 * Write database and test contents.
 */
void writeHDFTestData(tbox::Pointer<tbox::Database> db)
{
   if (db.isNull()) {
      tbox::perr << "FAILED: - Test #0-write: database to write to is null" << endl;
      tbox::MPI::abort();
   }

   /*
    * Build database hierarchy and test.
    */

   tbox::Pointer<tbox::Database> arraydb = db->putDatabase("Array Entries");

   tbox::Pointer<tbox::Database> scalardb = db->putDatabase("Scalar Entries");
   tbox::Pointer<tbox::Database> scalardb_empty = scalardb->putDatabase("Empty");
   tbox::Pointer<tbox::Database> scalardb_full = scalardb->putDatabase("Full");

   if (arraydb.isNull()) {
      tbox::perr << "FAILED: - Test #1a-write: `arraydb' is null" << endl;
      tbox::MPI::abort();
   }
   if (scalardb.isNull()) {
      tbox::perr << "FAILED: - Test #1b-write: `scalardb' is null" << endl;
      tbox::MPI::abort();
   }
   if (scalardb_empty.isNull()) {
      tbox::perr << "FAILED: - Test #1c-write: `scalardb_empty' is null" << endl;
      tbox::MPI::abort();
   }
   if (scalardb_full.isNull()) {
      tbox::perr << "FAILED: - Test #1d-write: `scalardb_full' is null" << endl;
      tbox::MPI::abort();
   }

   /*
    * Set array values and write to database hierarchy.
    */

   tbox::Array<dcomplex> arraydb_compArray(3);
   arraydb_compArray[0] = arraydb_compArray0;
   arraydb_compArray[1] = arraydb_compArray1;
   arraydb_compArray[2] = arraydb_compArray2;

   tbox::Array<bool> arraydb_boolArray(3);
   arraydb_boolArray[0] = arraydb_boolArray0;
   arraydb_boolArray[1] = arraydb_boolArray1;
   arraydb_boolArray[2] = arraydb_boolArray2;

   tbox::Array<int> arraydb_intArray(5);
   arraydb_intArray[0] = arraydb_intArray0;
   arraydb_intArray[1] = arraydb_intArray1;
   arraydb_intArray[2] = arraydb_intArray2;
   arraydb_intArray[3] = arraydb_intArray3;
   arraydb_intArray[4] = arraydb_intArray4;

   tbox::Array<string> arraydb_stringArray(3);
   arraydb_stringArray[0] = arraydb_stringArray0;
   arraydb_stringArray[1] = arraydb_stringArray1;
   arraydb_stringArray[2] = arraydb_stringArray2;

   tbox::Array<float> arraydb_floatArray(5);
   arraydb_floatArray[0] = arraydb_floatArray0;
   arraydb_floatArray[1] = arraydb_floatArray1;
   arraydb_floatArray[2] = arraydb_floatArray2;
   arraydb_floatArray[3] = arraydb_floatArray3;
   arraydb_floatArray[4] = arraydb_floatArray4;

   tbox::Array<double> arraydb_doubleArray(6);
   arraydb_doubleArray[0] = arraydb_doubleArray0;
   arraydb_doubleArray[1] = arraydb_doubleArray1;
   arraydb_doubleArray[2] = arraydb_doubleArray2;
   arraydb_doubleArray[3] = arraydb_doubleArray3;
   arraydb_doubleArray[4] = arraydb_doubleArray4;
   arraydb_doubleArray[5] = arraydb_doubleArray5;

   tbox::Array<char> arraydb_charArray(2);
   arraydb_charArray[0] = arraydb_charArray0;
   arraydb_charArray[1] = arraydb_charArray1;

   tbox::Array<tbox::DatabaseBox> arraydb_boxArray(3);
   arraydb_boxArray[0] = arraydb_boxArray0;
   arraydb_boxArray[1] = arraydb_boxArray1;
   arraydb_boxArray[2] = arraydb_boxArray2;

   db->putFloat("float_val", db_float_val);
   db->putInteger("int_val", db_int_val); 

   arraydb->putComplexArray("ComplexArray", arraydb_compArray);
   arraydb->putDatabaseBoxArray("BoxArray", arraydb_boxArray);
   arraydb->putBoolArray("BoolArray", arraydb_boolArray);
   arraydb->putIntegerArray("IntArray", arraydb_intArray);
   arraydb->putStringArray("StringArray", arraydb_stringArray);
   arraydb->putFloatArray("FloatArray", arraydb_floatArray);
   arraydb->putDoubleArray("DoubleArray", arraydb_doubleArray);
   arraydb->putCharArray("CharArray", arraydb_charArray);

   scalardb->putFloat("float1", scalardb_float1);
   scalardb->putFloat("float2", scalardb_float2);
   scalardb->putFloat("float3", scalardb_float3);

   scalardb_full->putDouble("thisDouble", scalardb_full_thisDouble);
   scalardb_full->putComplex("thisComplex", scalardb_full_thisComplex);
   scalardb_full->putInteger("thisInt", scalardb_full_thisInt);
   scalardb_full->putFloat("thisFloat", scalardb_full_thisFloat);
   scalardb_full->putBool("thisBool", scalardb_full_thisBool);
   scalardb_full->putString("thisString", scalardb_full_thisString);
   scalardb_full->putChar("thisChar", scalardb_full_thisChar);
   scalardb_full->putDatabaseBox("thisBox", scalardb_full_thisBox);


   testDatabaseContents(db, "write");
}

/*
 * Read database and test contents.
 */
void readHDFTestData(tbox::Pointer<tbox::Database> db)
{
   testDatabaseContents(db, "read");
}

/*
 * Test contents of database.
 */
void testDatabaseContents(tbox::Pointer<tbox::Database> db,
                          const string& tag)
{

   if (db.isNull()) {
      tbox::perr << "FAILED: - Test #0-" << tag 
           << ": database to read from is null" << endl;
      ++number_of_failures;
   }

   tbox::Pointer<tbox::Database> arraydb = db->getDatabase("Array Entries");

   tbox::Pointer<tbox::Database> scalardb = db->getDatabase("Scalar Entries");
   tbox::Pointer<tbox::Database> scalardb_empty = scalardb->getDatabase("Empty");
   tbox::Pointer<tbox::Database> scalardb_full = scalardb->getDatabase("Full");

   if (arraydb.isNull()) {
      tbox::perr << "FAILED: - Test #1a-" << tag
           << ": `arraydb' is null" << endl;
      ++number_of_failures;
   }
   if (scalardb.isNull()) {
      tbox::perr << "FAILED: - Test #1b-" << tag 
           << ": `scalardb' is null" << endl;
      ++number_of_failures;
   }
   if (scalardb_empty.isNull()) {
      tbox::perr << "FAILED: - Test #1c-" << tag 
           << ": `scalardb_empty' is null" << endl;
      ++number_of_failures;
   }
   if (scalardb_full.isNull()) {
      tbox::perr << "FAILED: - Test #1d-" << tag 
           << ": `scalardb_full' is null" << endl;
      ++number_of_failures;
   }

   tbox::Array<string> dbkeys = db->getAllKeys();
   tbox::Array<string> arraydbkeys = arraydb->getAllKeys();
   tbox::Array<string> scalardbkeys = scalardb->getAllKeys();
   tbox::Array<string> scalardb_emptykeys = scalardb_empty->getAllKeys();
   tbox::Array<string> scalardb_fullkeys = scalardb_full->getAllKeys();
 
   int i, nkeys; 

   if (dbkeys.getSize() != 4) {
      tbox::perr << "FAILED: - Test #2a-" << tag 
           << ": # `db' keys wrong" << endl;
      ++number_of_failures;
   }
   nkeys =arraydbkeys.getSize();
   if (arraydbkeys.getSize() != 8) {
      tbox::perr << "FAILED: - Test #2b-" << tag 
           << ": # `arraydb' keys wrong"
	   << "\n\tFound " << nkeys << " keys:"
	   ;
      ++number_of_failures;
      for ( i=0; i<nkeys; ++i ) {
	 tbox::pout << "\n\t\t" << i << ": '" << arraydbkeys[i] << "'";
      }
      tbox::pout << endl;
   }
   if (scalardbkeys.getSize() != 5) {
      tbox::perr << "FAILED: - Test #2c-" << tag 
           << ": # `scalardb' keys wrong" << endl;
      ++number_of_failures;
   }
   if (scalardb_emptykeys.getSize() != 0) {
      tbox::perr << "FAILED: - Test #2d-" << tag 
           << ": # `scalardb_empty' keys wrong" << endl;
      ++number_of_failures;
   }
   if (scalardb_fullkeys.getSize() != 8) {
      tbox::perr << "FAILED: - Test #2e-" << tag 
           << ": `scalardb_full' is null" << endl;
      ++number_of_failures;
   }

  
   if (!db->isDatabase("Array Entries")) { 
      tbox::perr << "FAILED: - #3a-" << tag 
           << ": `Array Entries' not a database" << endl;
      ++number_of_failures;
   } 
   if (!db->isDatabase("Scalar Entries")) { 
      tbox::perr << "FAILED: - #3b-" << tag 
           << ": `Scalar Entries' not a database" << endl;
      ++number_of_failures;
   }
   float tdb_float_val = db->getFloat("float_val");
   if (tdb_float_val != db_float_val) {
      tbox::perr << "FAILED: - Test #3c-" << tag 
           << ": `RestartTester' database"
           << "\n   Returned `float_val' = " << tdb_float_val 
           << "  , Expected = " << db_float_val << endl;
      ++number_of_failures;
   }
   int tdb_int_val = db->getInteger("int_val");   
   if (tdb_int_val != db_int_val) {
      tbox::perr << "FAILED: - Test #3d-" << tag 
           << ": `RestartTester' database"
           << "\n   Returned `int_val' = " << tdb_int_val 
           << "  , Expected = " << db_int_val << endl;
      ++number_of_failures;
   }

   
   /*
    * Set array values to test database.
    */

   tbox::Array<dcomplex> arraydb_compArray(3);
   arraydb_compArray[0] = arraydb_compArray0;
   arraydb_compArray[1] = arraydb_compArray1;
   arraydb_compArray[2] = arraydb_compArray2;

   tbox::Array<bool> arraydb_boolArray(3);
   arraydb_boolArray[0] = arraydb_boolArray0;
   arraydb_boolArray[1] = arraydb_boolArray1;
   arraydb_boolArray[2] = arraydb_boolArray2;

   tbox::Array<int> arraydb_intArray(5);
   arraydb_intArray[0] = arraydb_intArray0;
   arraydb_intArray[1] = arraydb_intArray1;
   arraydb_intArray[2] = arraydb_intArray2;
   arraydb_intArray[3] = arraydb_intArray3;
   arraydb_intArray[4] = arraydb_intArray4;

   tbox::Array<string> arraydb_stringArray(3);
   arraydb_stringArray[0] = arraydb_stringArray0;
   arraydb_stringArray[1] = arraydb_stringArray1;
   arraydb_stringArray[2] = arraydb_stringArray2;

   tbox::Array<float> arraydb_floatArray(5);
   arraydb_floatArray[0] = arraydb_floatArray0;
   arraydb_floatArray[1] = arraydb_floatArray1;
   arraydb_floatArray[2] = arraydb_floatArray2;
   arraydb_floatArray[3] = arraydb_floatArray3;
   arraydb_floatArray[4] = arraydb_floatArray4;

   tbox::Array<double> arraydb_doubleArray(6);
   arraydb_doubleArray[0] = arraydb_doubleArray0;
   arraydb_doubleArray[1] = arraydb_doubleArray1;
   arraydb_doubleArray[2] = arraydb_doubleArray2;
   arraydb_doubleArray[3] = arraydb_doubleArray3;
   arraydb_doubleArray[4] = arraydb_doubleArray4;
   arraydb_doubleArray[5] = arraydb_doubleArray5;

   tbox::Array<char> arraydb_charArray(2);
   arraydb_charArray[0] = arraydb_charArray0;
   arraydb_charArray[1] = arraydb_charArray1;

   tbox::Array<tbox::DatabaseBox> arraydb_boxArray(3);
   arraydb_boxArray[0] = arraydb_boxArray0;
   arraydb_boxArray[1] = arraydb_boxArray1;
   arraydb_boxArray[2] = arraydb_boxArray2;

   int tsize = 0;

   tbox::Array<dcomplex> tarraydb_compArray = 
      arraydb->getComplexArray("ComplexArray");
   tsize = tarraydb_compArray.getSize();
   if (tsize != arraydb_compArray.getSize()) {
      tbox::perr << "FAILED: - Test #4a-" << tag 
           << ": `Array Entries' database"
           << "\n   Returned `ComplexArray' size = " << tsize 
           << "  , Expected = " << arraydb_compArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_compArray[i] != arraydb_compArray[i]) {
         tbox::perr << "FAILED: - Test #4b-" << tag 
              << ": `Array Entries' database"
              << "\n   `ComplexArray' entry " << i << " incorrect" << endl; 
	 ++number_of_failures;
      }
   }
   tbox::Array<bool> tarraydb_boolArray = 
      arraydb->getBoolArray("BoolArray");
   tsize = tarraydb_boolArray.getSize();
   if (tsize != arraydb_boolArray.getSize()) {
      tbox::perr << "FAILED: - Test #4c-" << tag 
           << ": `Array Entries' database" 
           << "\n   Returned `BoolArray' size = " << tsize 
           << "  , Expected = " << arraydb_boolArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_boolArray[i] != arraydb_boolArray[i]) {
         tbox::perr << "FAILED: - Test #4d-" << tag 
              << ": `Array Entries' database"
              << "\n   `BoolArray' entry " << i << " incorrect"
              << "\n   " << tarraydb_boolArray[i] << " should be " << arraydb_boolArray[i] << endl;
	 ++number_of_failures;
      }
   }
   tbox::Array<int> tarraydb_intArray = 
      arraydb->getIntegerArray("IntArray");
   tsize = tarraydb_intArray.getSize();
   if (tsize != arraydb_intArray.getSize()) {
      tbox::perr << "FAILED: - Test #4e-" << tag 
           << ": `Array Entries' database" 
           << "\n   Returned `IntArray' size = " << tsize 
           << "  , Expected = " << arraydb_intArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_intArray[i] != arraydb_intArray[i]) {
         tbox::perr << "FAILED: - Test #4f-" << tag 
              << ": `Array Entries' database"
              << "\n   `IntArray' entry " << i << " incorrect" << endl;
	 ++number_of_failures;
      }
   }
   tbox::Array<string> tarraydb_stringArray = 
      arraydb->getStringArray("StringArray");
   tsize = tarraydb_stringArray.getSize();
   if (tsize != arraydb_stringArray.getSize()) {
      tbox::perr << "FAILED: - Test #4g-" << tag 
           << ": `Array Entries' database" 
           << "\n   Returned `StringArray' size = " << tsize 
           << "  , Expected = " << arraydb_stringArray[i] << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_stringArray[i] != arraydb_stringArray[i]) {
         tbox::perr << "FAILED: - Test #4h-" << tag 
              << ": `Array Entries' database"
              << "\n   `StringArray' entry " << i << " incorrect" << endl;
	 ++number_of_failures;
      }
   }
   tbox::Array<float> tarraydb_floatArray = 
      arraydb->getFloatArray("FloatArray");
   tsize = tarraydb_floatArray.getSize();
   if (tsize != arraydb_floatArray.getSize()) {
      tbox::perr << "FAILED: - Test #4i-" << tag 
           << ": `Array Entries' database" 
           << "\n   Returned `FloatArray' size = " << tsize 
           << "  , Expected = " << arraydb_floatArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_floatArray[i] != arraydb_floatArray[i]) {
         tbox::perr << "FAILED: - Test #4j-" << tag 
              << ": `Array Entries' database"
              << "\n   `FloatArray' entry " << i << " incorrect" << endl;
	 ++number_of_failures;
      }
   }
   tbox::Array<double> tarraydb_doubleArray = 
      arraydb->getDoubleArray("DoubleArray");
   tsize = tarraydb_doubleArray.getSize();
   if (tsize != arraydb_doubleArray.getSize()) {
      tbox::perr << "FAILED: - Test #4k-" << tag 
           << ": `Array Entries' database" 
           << "\n   Returned `DoubleArray' size = " << tsize 
           << "  , Expected = " << arraydb_doubleArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_doubleArray[i] != arraydb_doubleArray[i]) {
         tbox::perr << "FAILED: - Test #4l-" << tag 
              << ": `Array Entries' database"
              << "\n   `DoubleArray' entry " << i << " incorrect" << endl;
	 ++number_of_failures;
      }
   }
   tbox::Array<char> tarraydb_charArray =
      arraydb->getCharArray("CharArray");
   tsize = tarraydb_charArray.getSize();
   if (tsize != arraydb_charArray.getSize()) {
      tbox::perr << "FAILED: - Test #4m-" << tag 
           << ": `Array Entries' database"
           << "\n   Returned `CharArray' size = " << tsize
           << "  , Expected = " << arraydb_charArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_charArray[i] != arraydb_charArray[i]) {
         tbox::perr << "FAILED: - Test #4l-" << tag 
              << ": `Array Entries' database"
              << "\n   `CharArray' entry " << i << " incorrect" << endl;
	 ++number_of_failures;
      }
   }
   tbox::Array<tbox::DatabaseBox> tarraydb_boxArray = 
      arraydb->getDatabaseBoxArray("BoxArray");
   tsize = tarraydb_boxArray.getSize();
   if (tsize != arraydb_boxArray.getSize()) {
      tbox::perr << "FAILED: - Test #4o-" << tag 
           << ": `Array Entries' database" 
           << "\n   Returned `BoxArray' size = " << tsize 
           << "  , Expected = " << arraydb_boxArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (!(tarraydb_boxArray[i] == arraydb_boxArray[i])) {
         tbox::perr << "FAILED: - Test #4p-" << tag 
              << ": `Array Entries' database"
              << "\n   `BoxArray' entry " << i << " incorrect" << endl;
	 ++number_of_failures;
      }
   }

 
   if (!scalardb->isDatabase("Empty")) {
      tbox::perr << "FAILED: - #5a-" << tag 
           << ": `Empty' not a database" << endl;
      ++number_of_failures;
   }
   if (!scalardb->isDatabase("Full")) {
      tbox::perr << "FAILED: - #5b-" << tag 
           << ": `Full' not a database" << endl;
      ++number_of_failures;
   } 
   double tscalardb_float1 = scalardb->getFloat("float1");
   if (tscalardb_float1 != scalardb_float1) {
      tbox::perr << "FAILED: - Test #5c-" << tag 
           << ": `Scalar Entries' database"
           << "\n   Returned `float1' = " << tscalardb_float1
           << "  , Expected = " << scalardb_float1 << endl;
      ++number_of_failures;
   }
   double tscalardb_float2 = scalardb->getFloat("float2");
   if (tscalardb_float2 != scalardb_float2) {
      tbox::perr << "FAILED: - Test #5d-" << tag 
           << ": `Scalar Entries' database"
           << "\n   Returned `float2' = " << tscalardb_float2
           << "  , Expected = " << scalardb_float2 << endl;
      ++number_of_failures;
   }
   double tscalardb_float3 = scalardb->getFloat("float3");
   if (tscalardb_float3 != scalardb_float3) {
      tbox::perr << "FAILED: - Test #5e-" << tag 
           << ": `Scalar Entries' database"
           << "\n   Returned `float3' = " << tscalardb_float3
           << "  , Expected = " << scalardb_float3 << endl;
      ++number_of_failures;
   }


   double tscalardb_full_thisDouble = scalardb_full->getDouble("thisDouble");
   if (tscalardb_full_thisDouble != scalardb_full_thisDouble) {
      tbox::perr << "FAILED: - Test #6a-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisDouble' = " << tscalardb_full_thisDouble
           << "  , Expected = " << scalardb_full_thisDouble << endl;
      ++number_of_failures;
   }
   dcomplex tscalardb_full_thisComplex = 
      scalardb_full->getComplex("thisComplex");
   if (tscalardb_full_thisComplex != scalardb_full_thisComplex) {
      tbox::perr << "FAILED: - Test #6b-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisComplex' = " << tscalardb_full_thisComplex
           << "  , Expected = " << scalardb_full_thisComplex << endl;
      ++number_of_failures;
   }
   int tscalardb_full_thisInt = scalardb_full->getInteger("thisInt");
   if (tscalardb_full_thisInt != scalardb_full_thisInt) {
      tbox::perr << "FAILED: - Test #6c-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisInt' = " << tscalardb_full_thisInt
           << "  , Expected = " << scalardb_full_thisInt << endl;
      ++number_of_failures;
   }
   float tscalardb_full_thisFloat = scalardb_full->getFloat("thisFloat");
   if (tscalardb_full_thisFloat != scalardb_full_thisFloat) {
      tbox::perr << "FAILED: - Test #6d-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisFloat' = " << tscalardb_full_thisFloat
           << "  , Expected = " << scalardb_full_thisFloat << endl;
      ++number_of_failures;
   }
   bool tscalardb_full_thisBool = scalardb_full->getBool("thisBool");
   if (tscalardb_full_thisBool != scalardb_full_thisBool) {
      tbox::perr << "FAILED: - Test #6e-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisBool' = " << tscalardb_full_thisBool
           << "  , Expected = " << scalardb_full_thisBool << endl;
      ++number_of_failures;
   }
   string tscalardb_full_thisString = scalardb_full->getString("thisString");
   if (tscalardb_full_thisString != scalardb_full_thisString) {
      tbox::perr << "FAILED: - Test #6f-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisString' = " << tscalardb_full_thisString
           << "  , Expected = " << scalardb_full_thisString << endl;
      ++number_of_failures;
   }
   char tscalardb_full_thisChar = scalardb_full->getChar("thisChar");
   if (tscalardb_full_thisChar != scalardb_full_thisChar) {
      tbox::perr << "FAILED: - Test #6g-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisChar' = " << tscalardb_full_thisChar
           << "  , Expected = " << scalardb_full_thisChar << endl;
      ++number_of_failures;
   }
   tbox::DatabaseBox tscalardb_full_thisBox = scalardb_full->getDatabaseBox("thisBox");
   if (!(tscalardb_full_thisBox == scalardb_full_thisBox)) {
      tbox::perr << "FAILED: - Test #6h-" << tag 
           << ": `Full' database"
           << "\n   Returned `thisBox' does not match Expected value" << endl;
      ++number_of_failures;
   }

}
