//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/branches/smith84/source/test/restartdb/main-HDF5.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2105 $
// Modified:    $LastChangedDate: 2008-04-01 12:36:17 -0700 (Tue, 01 Apr 2008) $
// Description: Global values for the restart tests
//


#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"
#include "tbox/Array.h"
#include "tbox/DatabaseBox.h"
#include "tbox/Complex.h"
#include "tbox/SAMRAI_MPI.h"
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

dcomplex arraydb_dcomplexArray0 = dcomplex(1,2);
dcomplex arraydb_dcomplexArray1 = dcomplex(2,3);
dcomplex arraydb_dcomplexArray2 = dcomplex(3,4);
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
