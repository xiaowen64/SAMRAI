//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/branches/smith84/source/test/restartdb/main-HDF5.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2104 $
// Modified:    $LastChangedDate: 2008-04-01 10:02:41 -0700 (Tue, 01 Apr 2008) $
// Description: Some simple generic database test functions 
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
extern int number_of_failures;

/**
 * Write database and test contents.
 */
void setupTestData(void);

/**
 * Write database and test contents.
 */
void writeTestData(tbox::Pointer<tbox::Database> db);

/**
 * Read database and test contents.
 */
void readTestData(tbox::Pointer<tbox::Database> db);

/**
 * Test contents of database.
 */
void testDatabaseContents(tbox::Pointer<tbox::Database> db,
                          const string& tag);
