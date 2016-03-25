//
// File:        main-var_db.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Main program to test variable database operations
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include "tbox/PIO.h"
#include "tbox/MPI.h"
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif

#include "tbox/SAMRAIManager.h"

#include "CellVariable.h"
#include "FaceVariable.h"
#include "IntVector.h"
#include "OuterfaceVariable.h"
#include "tbox/Pointer.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include <string>
using namespace std;

using namespace SAMRAI;

int main( int argc, char *argv[] ) {
   tbox::MPI::init(&argc, &argv); 
   tbox::SAMRAIManager::startup();
// tbox::PIO::logOnlyNodeZero("main-var_db.log");
   tbox::PIO::logAllNodes("vdbtest.log");

   hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();

   tbox::Pointer<hier::VariableContext> current_context = var_db->getContext("CURRENT");

   hier::IntVector<NDIM> nghosts(4); 
   hier::IntVector<NDIM> fluxghosts(1); 
   hier::IntVector<NDIM> zero_ghosts(0);
  
   /* State variable */
   tbox::Pointer< pdat::CellVariable<NDIM,double> > uval = 
                                   new pdat::CellVariable<NDIM,double>("uval",1);

   /* Flux variable */
   tbox::Pointer< pdat::FaceVariable<NDIM,double> > flux =
                                   new pdat::FaceVariable<NDIM,double>("flux",1);

   /* Register uval using ready made context 
    * and by creating a new one on the fly 
    */
   const int uval_current_id = 
      var_db->registerVariableAndContext(uval,
                                         current_context,
                                         zero_ghosts);

   /* Register flux using ready created context and by using 
    * context already in the database but one we don't have pointer to.
    */
// use of void eliminates compiler warning
// const int flux_current_id = 
   (void)
      var_db->registerVariableAndContext(flux,
                                         current_context,
                                         zero_ghosts);

   const int flux_scratch_id = 
      var_db->registerVariableAndContext(flux,
                                         var_db->getContext("SCRATCH"),
                                         nghosts);

   tbox::Pointer< pdat::OuterfaceVariable<NDIM,double> > fluxsum = 
      new pdat::OuterfaceVariable<NDIM,double>("fluxsum", 1);

// use of void eliminates compiler warning
// const int fluxsum_current_id = 
   (void)
      var_db->registerVariableAndContext(fluxsum,
                                         current_context,
                                         zero_ghosts);

   /*
    * Now, run the variable database through the ringer... 
    */

   tbox::plog << "\n\nPrintout #1 of hier::Variable tbox::Database (after initial registration)..." << endl;
   var_db->printClassData(tbox::plog);

   // Test #1: Check Context functions...

   // Test #1a: hier::VariableDatabase::checkContextExists()
   tbox::plog << "Test #1a: hier::VariableDatabase::checkContextExists()..." << endl;
   if (!var_db->checkContextExists("SCRATCH")) {
      tbox::perr << "FAILED: - Test #1a: hier::VariableDatabase::checkContextExists()\n"
           << "SCRATCH context added to var_db, but not found" << endl;
   }

   // Test #1b: hier::VariableDatabase::checkContextExists()
   tbox::plog << "Test #1b: hier::VariableDatabase::checkContextExists()..." << endl;
   if (!var_db->checkContextExists("CURRENT")) {
      tbox::perr << "FAILED: - Test #1b: hier::VariableDatabase::checkContextExists()\n"
           << "CURRENT context added to var_db, but not found" << endl;
   }

   // Test #1c: hier::VariableDatabase::checkContextExists()
   tbox::plog << "Test #1c: hier::VariableDatabase::checkContextExists()..." << endl;
   if (var_db->checkContextExists("dummy")) {
      tbox::perr << "FAILED: - Test #1c: hier::VariableDatabase::checkContextExists()\n"
           << "dummy context not added to var_db, but found" << endl;
   }

   // Adding dummy context to variable databse...
   tbox::Pointer<hier::VariableContext> dummy_ctxt = var_db->getContext("dummy");

   // Test #1d: hier::VariableDatabase::checkContextExists()
   tbox::plog << "Test #1d: hier::VariableDatabase::checkContextExists()..." << endl;
   if (!var_db->checkContextExists("dummy")) {
      tbox::perr << "FAILED: - Test #1d: hier::VariableDatabase::checkContextExists()\n"
           << "dummy context added to var_db, but not found" << endl;
   }


   // Test #2,3: Check hier::Variable functions...

   // Test #2a: hier::VariableDatabase::checkVariableExists()
   tbox::plog << "Test #2a: hier::VariableDatabase::checkVariableExists()..." << endl;
   if (!var_db->checkVariableExists("uval")) {
      tbox::perr << "FAILED: - Test #2a: hier::VariableDatabase::checkVariableExists()\n"
           << "uval variable added to var_db, but not found" << endl;
   }

   // Test #2b: hier::VariableDatabase::checkVariableExists()
   tbox::plog << "Test #2b: hier::VariableDatabase::checkVariableExists()..." << endl;
   if (!var_db->checkVariableExists("flux")) {
      tbox::perr << "FAILED: - Test #2b: hier::VariableDatabase::checkVariableExists()\n"
           << "flux variable added to var_db, but not found" << endl;
   }
   
   // Test #2c: hier::VariableDatabase::checkVariableExists()
   tbox::plog << "Test #2c: hier::VariableDatabase::checkVariableExists()..." << endl;
   if (!var_db->checkVariableExists("fluxsum")) {
      tbox::perr << "FAILED: - Test #2c: hier::VariableDatabase::checkVariableExists()\n"
           << "fluxsum variable added to var_db, but not found" << endl;
   }

   // Test #2d: hier::VariableDatabase::checkVariableExists()
   tbox::plog << "Test #2d: hier::VariableDatabase::checkVariableExists()..." << endl;
   if (var_db->checkVariableExists("dummy")) {
      tbox::perr << "FAILED: - Test #2d: hier::VariableDatabase::checkVariableExists()\n"
           << "dummy variable not added to var_db, but found" << endl;
   }

   // Test #3a: hier::VariableDatabase::getVariable()
   tbox::plog << "Test #3a: hier::VariableDatabase::getVariable()..." << endl;
   tbox::Pointer<hier::Variable<NDIM> > tvar_uval = var_db->getVariable("uval");
   if (tvar_uval.isNull()) {
      tbox::perr << "FAILED: - Test #3a: hier::VariableDatabase::getVariable()\n"
           << "uval variable added to var_db, but returning NULL" << endl;
   }

   // Test #3b: hier::VariableDatabase::getVariable()
   tbox::plog << "Test #3b: hier::VariableDatabase::getVariable()..." << endl;
   tbox::Pointer<hier::Variable<NDIM> > tvar_flux = var_db->getVariable("flux");
   if (tvar_flux.isNull()) {
      tbox::perr << "FAILED: - Test #3b: hier::VariableDatabase::getVariable()\n"
           << "flux variable added to var_db, but returning NULL" << endl;
   }

   // Test #3c: hier::VariableDatabase::getVariable()
   tbox::plog << "Test #3c: hier::VariableDatabase::getVariable()..." << endl;
   tbox::Pointer<hier::Variable<NDIM> > tvar_fluxsum = var_db->getVariable("fluxsum");
   if (tvar_fluxsum.isNull()) {
      tbox::perr << "FAILED: - Test #3c: hier::VariableDatabase::getVariable()\n"
           << "fluxsum variable added to var_db, but returning NULL" << endl;
   }

   // Test #3d: hier::VariableDatabase::getVariable()
   tbox::plog << "Test #3d: hier::VariableDatabase::getVariable()..." << endl;
   //   tbox::perr << "Attempt to get variable named dummy..." << endl;
   tbox::Pointer<hier::Variable<NDIM> > tvar_dummy = var_db->getVariable("dummy");
   if (!tvar_dummy.isNull()) {
      tbox::perr << "FAILED: - Test #3d: hier::VariableDatabase::getVariable()\n"
           << "dummy variable not added to var_db, but not returning NULL" 
           << endl;
   }

   // Test #4: Check instance identifier assignments

   // Test #4a: hier::Variable::getInstanceIdentifier()
   tbox::plog << "Test #4a: hier::Variable::getInstanceIdentifier()..." << endl;
   int uval_id = tvar_uval->getInstanceIdentifier();
   if (uval_id != 0) {
      tbox::perr << "FAILED: - Test #4a: hier::Variable::getInstanceIdentifier()\n"
           << "uval should have id = 0, but has id = " << uval_id <<endl;
   }

   // Test #4b: hier::Variable::getInstanceIdentifier()
   tbox::plog << "Test #4b: hier::Variable::getInstanceIdentifier()..." << endl;
   int flux_id = tvar_flux->getInstanceIdentifier();
   if (flux_id != 1) {
      tbox::perr << "FAILED: - Test #4b: hier::Variable::getInstanceIdentifier()\n"
           << "flux should have id = 1, but has id = " << flux_id <<endl;
   }

   // Test #4c: hier::Variable::getInstanceIdentifier()
   tbox::plog << "Test #4c: hier::Variable::getInstanceIdentifier()..." << endl;
   int fluxsum_id = tvar_fluxsum->getInstanceIdentifier();
   if (fluxsum_id != 2) {
      tbox::perr << "FAILED: - Test #4c: hier::Variable::getInstanceIdentifier()\n"
           << "fluxsum should have id = 2, but has id = " << fluxsum_id <<endl;
   }

   // Test #5: Attempt to register (uval,CURRENT) again
   tbox::plog << "Test #5: Attempt to register (uval,CURRENT) again..." << endl;
   tbox::Pointer<hier::VariableContext> tctxt_current = var_db->getContext("CURRENT");
   hier::IntVector<NDIM> tzero_ghosts(0);
   int ti = var_db->registerVariableAndContext(
      tvar_uval, tctxt_current, tzero_ghosts);
   if (ti != 0) {
      tbox::perr << "FAILED: - Test #5: Re-registering a variable and context\n"
           << "Original id = 0 should be returned, but "
           << ti << "returned" << endl;
   }

   // Test #6a: hier::VariableDatabase::mapVariableAndContextToIndex()
   tbox::plog << "Test #6a: hier::VariableDatabase::mapVariableAndContextToIndex()..." << endl;
   ti = var_db->mapVariableAndContextToIndex(tvar_uval, tctxt_current);
   if (ti != 0) {
     tbox::perr << "FAILED: - Test #6a: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
          << "(uval,CURRENT) should be mapped to 0, but is mapped to " 
          << ti << endl;
   }

   // Test #6b: hier::VariableDatabase::mapVariableAndContextToIndex()
   tbox::plog << "Test #6b: hier::VariableDatabase::mapVariableAndContextToIndex()..." << endl;
   tvar_uval = var_db->getVariable("uval");
   tbox::Pointer<hier::VariableContext> tctxt_scratch = var_db->getContext("SCRATCH");
   ti = var_db->mapVariableAndContextToIndex(tvar_uval, tctxt_scratch); 
   if (ti != -1) {
     tbox::perr << "FAILED: - Test #6b: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
          << "(uval,SCRATCH) should be mapped to -1, but is mapped to " 
          << ti << endl;
   }

   // Test #6c: hier::VariableDatabase::mapVariableAndContextToIndex()
   tbox::plog << "Test #6c: hier::VariableDatabase::mapVariableAndContextToIndex()..." << endl;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > dummy_var =
   new pdat::CellVariable<NDIM,double>("dummy",3);
   tctxt_scratch = var_db->getContext("SCRATCH");
   ti = var_db->mapVariableAndContextToIndex(dummy_var, tctxt_scratch);
   if (ti != -1) {
     tbox::perr << "FAILED: - Test #6c: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
          << "(dummy,SCRATCH) should be mapped to -1, but is mapped to " 
          << ti << endl;
   }

   // Test #6d: hier::VariableDatabase::mapVariableAndContextToIndex()
   tbox::plog << "Test #6d: hier::VariableDatabase::mapVariableAndContextToIndex()..." << endl;
   tvar_uval = var_db->getVariable("uval");
   tbox::Pointer<hier::VariableContext> tctxt_random = new hier::VariableContext("RANDOM");
   ti = var_db->mapVariableAndContextToIndex(tvar_uval, tctxt_random); 
   if (ti != -1) {
     tbox::perr << "FAILED: - Test #6d: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
          << "(uval,RANDOM) should be mapped to -1, but is mapped to " 
          << ti << endl;
   }

   // Test #7a: hier::VariableDatabase::mapIndexToVariableAndContext()
   tbox::plog << "Test #7a: hier::VariableDatabase::mapIndexToVariableAndContext()..." << endl;
   int search_id = 2;
   tbox::Pointer<hier::Variable<NDIM> > search_var;
   tbox::Pointer<hier::VariableContext> search_ctxt;
   string flux_variable("flux");
   string scratch_variable("SCRATCH");

   // searching for index = 2
   if (!var_db->mapIndexToVariableAndContext(
           search_id, search_var, search_ctxt)) {
      tbox::perr << "FAILED: - Test #7a.1: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
           << "Problem finding a (variable,context) pair for index = 2" 
           << endl;
      if (search_var->getName() != flux_variable) {
         tbox::perr << "FAILED: - Test #7a.2: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
              << "Returned var name should be \"flux\" but is "
              << search_var->getName() << endl;
      }
      if (search_ctxt->getName() != scratch_variable) {
         tbox::perr << "FAILED: - Test #7a.3: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
              << "Returned context name should be \"SCRATCH\" but is "
              << search_ctxt->getName() << endl;
      }
/*      tbox::plog << "Var name = " << tvar->getName() << " = flux?, " 
           << "Context name = " << tctxt->getName() << " = CURRENT?" << endl;
*/
/*   } else {
     tbox::plog << "Houston, we have a problem looking for index "
          << search_id << endl;
*/
   }

   // Test #7b: hier::VariableDatabase::mapIndexToVariableAndContext()
   tbox::plog << "Test #7b: hier::VariableDatabase::mapIndexToVariableAndContext()" << endl; 
   search_id = 20;
   
   // searching for index = 20
   if (var_db->mapIndexToVariableAndContext(
          search_id, search_var, search_ctxt)) {
      tbox::perr << "FAILED: - Test #7b.1: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
           << "Something maps to index = 20 when nothing should.\n"
           << "Variable name: " << search_var->getName() << "\n"
           << "Context name: " << search_ctxt->getName() 
           << endl;
   } else {

      if (!search_var.isNull()) {
         tbox::perr << "FAILED: - Test #7b.2: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
              << "search_var should be NULL" << endl;
      }
      if (!search_ctxt.isNull()) {
         tbox::perr << "FAILED: - Test #7b.3: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
              << "search_ctxt should be NULL" << endl;
      }

   }

   // Test #7c: hier::VariableDatabase::mapIndexToVariable()
   tbox::plog << "Test #7c: hier::VariableDatabase::mapIndexToVariable()..." << endl; 
   search_id = 2;

   // searching for index = 2
   if (!var_db->mapIndexToVariable(search_id, search_var)) {
      tbox::perr << "FAILED: - Test #7c.1: hier::VariableDatabase::mapIndexToVariable()\n"
           << "Problem finding a (variable) for index = 2"
           << endl;
      if (search_var->getName() != flux_variable) {
         tbox::perr << "FAILED: - Test #7c.2: hier::VariableDatabase::mapIndexToVariable()\n"
              << "Returned var name should be \"flux\" but is "
              << search_var->getName() << endl;
      }
   }

   // Test #7d: hier::VariableDatabase::mapIndexToVariable()
   tbox::plog << "Test #7d: hier::VariableDatabase::mapIndexToVariable()..." << endl;
   search_id = 20;

   // searching for index = 20
   if (var_db->mapIndexToVariable(search_id, search_var)) {
      tbox::perr << "FAILED: - Test #7d.1: hier::VariableDatabase::mapIndexToVariable()\n"
           << "Something maps to index = 20 when nothing should.\n"
           << "Variable name: " << search_var->getName() << "\n"
           << endl;
   } else {

      if (!search_var.isNull()) {
         tbox::perr << "FAILED: - Test #7d.2: hier::VariableDatabase::mapIndexToVariable()\n"
              << "search_var should be NULL" << endl;
      }

   }

   // Testing new registration functions and ability of variable 
   // database to add and delete descriptor indices...

   tbox::plog << "\nPrintout #2 of hier::Variable tbox::Database (after tests 1-7." 
        << "\n Only difference with previous printout should be addition"
        << "\n of \"dummy\" context).....\n";
   var_db->printClassData(tbox::plog);

   // Test #8a: Checking mapping in variable database
   tbox::plog << "Test #8a: Checking mapping in variable database..." << endl;
   if (!var_db->checkVariablePatchDataIndex(uval, uval_current_id)) {
      tbox::perr << "FAILED: - Test #8a: "
           << "VariableDatabase:checkVariablePatchDataIndex()\n"
           << "uval should be map to current id in patch descriptor" << endl;
   }

   // Test #8b: Checking restration of existing patch data index...
   tbox::plog << "Test #8b: Checking restration of existing patch data index..." << endl;
   int test_id = var_db->registerPatchDataIndex(uval, uval_current_id);
  
   if (test_id != uval_current_id) {
      tbox::perr << "FAILED: - Test #8b: "
           << "VariableDatabase:registerPatchDataIndex()\n"
           << "re-registering current uval id should return same id" << endl;
   }
   
   // Test #8c: Testing registration of new cloned factory to variable...
   tbox::plog << "Test #8c: Testing registration of new cloned factory to variable..." << endl;
   int new_id = 
      var_db->registerClonedPatchDataIndex(uval, uval_current_id);

   if ( (new_id < 0) || (new_id == uval_current_id) ) {
      tbox::perr << "FAILED: - Test #8c: "
           << "VariableDatabase:registerClonedPatchDataIndex()\n"
           << "cloning current uval id return invalid id" << endl;
   }

   // Test #8d: Testing registration of new cloned factory to variable...
   tbox::plog << "Test #8d: Testing registration of new cloned factory to variable..." << endl;
   tbox::Pointer<hier::Variable<NDIM> > tvar;
   if ( !var_db->mapIndexToVariable(new_id, tvar)
       || (tvar->getName() != "uval") ) {
      tbox::perr << "FAILED: - Test #8d: "
           << "VariableDatabase:mapIndexToVariable()\n"
           << "descriptor id = " << new_id 
           << " should now map to uval" << endl;
   }

   tbox::plog << "\nPrintout #3 of hier::Variable tbox::Database. (after tests 8a-8d"
        << "\n Descriptor index " << new_id 
        << " should map to newly-created uval descriptor index)" << endl;
   var_db->printClassData(tbox::plog);

   // Test #8e: Testing removal of new variable-descriptor id mapping
   var_db->removePatchDataIndex(new_id);

   tbox::plog << "\nPrintout #4 of hier::Variable tbox::Database. (after test 8d"
        << "\n Descriptor index " << new_id 
        << " should no longer be mapped to a variable)..." << endl;
   var_db->printClassData(tbox::plog);

   tvar.setNull();
   if (var_db->mapIndexToVariable(new_id, tvar)) {
      tbox::perr << "FAILED: - Test #8e: "
           << "VariableDatabase:removePatchDataIndex()\n"
           << "descriptor id = " << new_id
           << " should no longer map to uval variable" << endl;
   }

   // Test #8f-g: Testing whether inconsistent mapping is allowed...
   tbox::plog << "Test #8f-g: Testing whether inconsistent mapping is allowed..." << endl;
   if (!var_db->checkVariablePatchDataIndex(flux, flux_scratch_id)) {
      tbox::perr << "FAILED: - Test #8f: "
           << "VariableDatabase:checkVariablePatchDataIndex()\n"
           << "flux should be mapped to scratch flux id" << endl;
   }
   if (var_db->checkVariablePatchDataIndex(uval, flux_scratch_id)) {
      tbox::perr << "FAILED: - Test #8g: "
           << "VariableDatabase:checkVariablePatchDataIndex()\n"
           << "uval should not map to scratch flux id" << endl;
   }

   // Test #8h-i: Testing removal of existing descriptor id...
   tbox::plog << "Test #8h-i: Testing removal of existing descriptor id..." << endl;
   var_db->removePatchDataIndex(flux_scratch_id);

   int tindex = 
      var_db->mapVariableAndContextToIndex(flux,
                                           var_db->getContext("SCRATCH"));
   if (tindex > -1) {
      tbox::perr << "FAILED: - Test #8h: "
           << "VariableDatabase:removePatchDataIndex()\n"
           << "flux-SCRATCH mapping should no longer be in database" << endl;
   }
 
   tvar.setNull();
   if (var_db->mapIndexToVariable(flux_scratch_id, tvar)) {
      tbox::perr << "FAILED: - Test #8i: "
           << "VariableDatabase:mapIndexToVariable()\n"
           << "descriptor id = " << flux_scratch_id
           << " should no longer map to flux variable" << endl;
   }

   tbox::plog << "\nPrintout #5 of hier::Variable tbox::Database. (after tests 8e-i\n"
       << "Should be identical to first printout except for the "
       << "\naddition of \n\"dummy\" context and removal of " 
       << "flux scratch variable-context mapping" << endl;
   var_db->printClassData(tbox::plog);

   /*
    * Tests that end in program abort...
    */
#if 0
   // Abort Test #1
   tbox::Pointer< pdat::CellVariable<NDIM,double> > dummy =
   new pdat::CellVariable<NDIM,double>("uval",2);

   tbox::plog << "Attempt to add a different variable named uval." 
        << "This should bomb!!" << endl;
   var_db->addVariable(dummy);

   // Abort Test #2
   tbox::plog << "Attempt to register uval, CURRENT again w/ wrong ghosts." 
        << "This should bomb!!" << endl;
   tctxt = var_db->getContext("CURRENT");
   tvar = var_db->getVariable("uval");
   g = hier::IntVector<NDIM>(2);
   ti = var_db->registerVariableAndContext(tvar, tctxt, g);
   tbox::plog << "uval, CURRENT at index = " << ti << endl;

   // Abort Test #3
   tbox::plog << "Attempt to register uval with fake CURRENT context."
        << "This should bomb!!" << endl;
   tctxt = new hier::VariableContext("CURRENT");
   tvar = var_db->getVariable("uval");
   g = hier::IntVector<NDIM>(0);
   ti = var_db->registerVariableAndContext(tvar, tctxt, g);
   tbox::plog << "uval, fake CURRENT at index = " << ti << endl;

   // Abort Test #4
   tbox::plog << "Attempt to map uval to descriptor id for flux."
        << "This should bomb!!" << endl;
   var_db->addVariablePatchDataIndex(uval, flux_scratch_id);
#endif

   tbox::pout << "\nPASSED:  vdbtest" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();

   return(0); 
}
