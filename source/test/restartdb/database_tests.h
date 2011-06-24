/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Some simple generic database test functions 
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/RestartManager.h"
#include <string>

using namespace std;
using namespace SAMRAI;

// Number of (non-abortive) failures.
extern int number_of_failures;

/**
 * Write database and test contents.
 */
void
setupTestData(
   void);

/**
 * Write database and test contents.
 */
void
writeTestData(
   tbox::Pointer<tbox::Database> db);

/**
 * Read database and test contents.
 */
void
readTestData(
   tbox::Pointer<tbox::Database> db);

/**
 * Test contents of database.
 */
void
testDatabaseContents(
   tbox::Pointer<tbox::Database> db,
   const string& tag);
