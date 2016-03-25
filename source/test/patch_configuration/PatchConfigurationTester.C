//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/patch_configuration/PatchConfigurationTester.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Class to test patch configuration utility
//

#include "PatchConfigurationTester.h"


/*
 * Header files for SAMRAI library classes
 */
#include "BoxArray.h"
#include "BoundaryLookupTable.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "GridGeometry.h"
#include "Index.h"
#include "LoadBalancer.h"
#include "ProcessorMapping.h"
#include "tbox/Utilities.h"

/*************************************************************************
 *
 * Constructor and Destructor for PatchConfigurationTester class.
 *
 ************************************************************************/

PatchConfigurationTester::PatchConfigurationTester(
   const string& object_name,
   tbox::Pointer<tbox::InputDatabase> input_db,
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());
   TBOX_ASSERT(!hierarchy.isNull());
#endif

   d_object_name   = object_name;
   d_hierarchy  = hierarchy; 

   d_num_levels = 2;
   d_npatches_on_coarsest = hier::IntVector<NDIM>(1);

   if (input_db->keyExists("test_to_run")) {
      d_test_to_run = input_db->getString("test_to_run");
      if ( (d_test_to_run != "HIERARCHY_AT_ONCE") &&
           (d_test_to_run != "LEVELS_COARSE_TO_FINE") &&
           (d_test_to_run != "LEVELS_FINE_TO_COARSE") ) {
         TBOX_ERROR("Input error in " << d_object_name << ": "
                    << "'test_to_run' input string must be either "
                    << "'HIERARCHY_AT_ONCE',"
                    << "'LEVELS_COARSE_TO_FINE',"
                    << " or 'LEVELS_FINE_TO_COARSE'." << endl);
      }
   } else {
      TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nNo 'test_to_run' item found!" << endl);
   }

   if (input_db->keyExists("GriddingParameters")) {
      getGriddingParametersFromInput(
         input_db->getDatabase("GriddingParameters")); 
   } else {
      TBOX_ERROR("Input error in " << d_object_name << ":"
                 << "\nNo 'GriddingParameters' input found!" << endl);
   }

}

PatchConfigurationTester::~PatchConfigurationTester()
{
   if (d_test_to_run == "HIERARCHY_AT_ONCE") {
      d_patch_config_util->clear();
   }

   if (d_test_to_run == "LEVELS_COARSE_TO_FINE") {
      for (int ln = 0; ln < d_hierarchy->getNumberLevels(); ln++) {
         tbox::Pointer< hier::PatchLevel<NDIM> > level =
            d_hierarchy->getPatchLevel(ln);
         d_patch_config_util->clear(level);
      }
   }

   if (d_test_to_run == "LEVELS_FINE_TO_COARSE") {
      for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
         tbox::Pointer< hier::PatchLevel<NDIM> > level =
            d_hierarchy->getPatchLevel(ln);
         d_patch_config_util->clear(level);
      }
   }
}

/*************************************************************************
 *
 * Build patch hierarchy.
 *
 ************************************************************************/

void PatchConfigurationTester::buildPatchHierarchy()
{
   tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom = 
      d_hierarchy->getGridGeometry();
   const hier::BoxArray<NDIM> domain_boxes = 
      grid_geom->getPhysicalDomain();

   tbox::Pointer<tbox::InputDatabase> lb_input_db = 
      new tbox::InputDatabase("LoadBalancer");

   tbox::Array<int> num_proc(NDIM);
   for (int idim = 0; idim < NDIM; idim++) {
      num_proc[idim] = d_npatches_on_coarsest(idim);
   }

   lb_input_db->putIntegerArray("processor_layout", num_proc);

   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = 
      new mesh::LoadBalancer<NDIM>(lb_input_db);

   hier::BoxList<NDIM> in0_boxes(domain_boxes);
   hier::BoxArray<NDIM> level0_boxes;
   hier::ProcessorMapping level0_mapping;
   int level_num = 0;
   hier::IntVector<NDIM> cut_factor(hier::IntVector<NDIM>(1));
   hier::IntVector<NDIM> bad_interval(hier::IntVector<NDIM>(1));
 
   hier::Box<NDIM> domain_box = in0_boxes.getBoundingBox();
   hier::IntVector<NDIM> box_size = 
      domain_box.upper() - domain_box.lower() + hier::IntVector<NDIM>(1);
   hier::IntVector<NDIM> largest_patch_size = 
      box_size/d_npatches_on_coarsest;
   hier::IntVector<NDIM> smallest_patch_size = 
      largest_patch_size - hier::IntVector<NDIM>(1);

   tbox::plog << "domain_box = " << domain_box << endl;
   tbox::plog << "d_npatches_on_coarsest = " << d_npatches_on_coarsest << endl;
   tbox::plog << "processor_layout = ";
   for (int idim = 0; idim < NDIM; idim++) {
      tbox::plog << num_proc[idim] << " , ";
   }
   tbox::plog << endl;
   tbox::plog << "cut_factor = " << cut_factor << endl;
   tbox::plog << "bad_interval = " << bad_interval << endl;
   tbox::plog << "box_size = " << box_size << endl;
   tbox::plog << "largest_patch_size = " << largest_patch_size << endl;
   tbox::plog << "smallest_patch_size = " << smallest_patch_size << endl;
 
   load_balancer->loadBalanceBoxes(level0_boxes, level0_mapping,
                                   in0_boxes, d_hierarchy, level_num,
                                   grid_geom->getPhysicalDomain(),
                                   hier::IntVector<NDIM>(1),
                                   smallest_patch_size,
                                   largest_patch_size,
                                   cut_factor, bad_interval);

#if 0 
   tbox::plog << "\n\nPatch-Processor Mapping - Level 0" << endl;
   tbox::Array<int> mapping = level0_mapping.getProcessorMapping();
   for (int i = 0; i < mapping.getSize(); i++) {
      tbox::plog << "  Patch<NDIM>: " << i << " : " << level0_boxes.getBox(i) 
           << " : " << mapping[i] << endl;
   }
#endif

   hier::IntVector<NDIM> ratio_to_level_zero(1);
 
   d_hierarchy->makeNewPatchLevel(level_num, ratio_to_level_zero,
                                  level0_boxes, level0_mapping);
 
   for (int ln = 1; ln < d_num_levels; ln++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln-1);
      hier::BoxArray<NDIM> coarser_level_boxes = coarser_level->getBoxes();
 
      ratio_to_level_zero = coarser_level->getRatio() * d_ratio_to_coarser[ln];
 
      hier::BoxList<NDIM> new_level_domain(coarser_level_boxes);
      new_level_domain.refine(d_ratio_to_coarser[ln]);
 
      hier::BoxArray<NDIM> new_level_boxes;
      hier::ProcessorMapping new_level_mapping;

      hier::BoxArray<NDIM> physical_domain;
      grid_geom->computePhysicalDomain(physical_domain,
                                       ratio_to_level_zero);

      cut_factor = d_ratio_to_coarser[ln];
 
      load_balancer->loadBalanceBoxes(new_level_boxes, new_level_mapping,
                                      new_level_domain, d_hierarchy, ln,
                                      physical_domain,
                                      ratio_to_level_zero,
                                      smallest_patch_size,
                                      largest_patch_size,
                                      cut_factor, bad_interval);

#if 0 
      tbox::plog << "\n\nPatch-Processor Mapping - Level " << ln << endl;
      tbox::Array<int> mapping = new_level_mapping.getProcessorMapping();
      for (int i = 0; i < mapping.getSize(); i++) {
         tbox::plog << "  Patch<NDIM>: " << i << " : " << new_level_boxes.getBox(i) 
              << " : " << mapping[i] << endl;
      }
#endif
 
      d_hierarchy->makeNewPatchLevel(ln, ratio_to_level_zero,
                                     new_level_boxes, new_level_mapping);
   }

}

/*************************************************************************
 *
 * Setup information about patch configuration.
 *
 ************************************************************************/

void PatchConfigurationTester::setupPatchConfiguration()
{

   d_patch_config_util = 
      new hier::PatchConfigurationUtilities<NDIM>("Tester Utils",
                                                  d_hierarchy);

   if (d_test_to_run == "HIERARCHY_AT_ONCE") {
      d_patch_config_util->initialize();
   } 

   if (d_test_to_run == "LEVELS_COARSE_TO_FINE") {
      for (int ln = 0; ln < d_hierarchy->getNumberLevels(); ln++) {
         tbox::Pointer< hier::PatchLevel<NDIM> > level = 
            d_hierarchy->getPatchLevel(ln);
         d_patch_config_util->initialize(level);
      }
   }
   
   if (d_test_to_run == "LEVELS_FINE_TO_COARSE") {
      for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
         tbox::Pointer< hier::PatchLevel<NDIM> > level = 
            d_hierarchy->getPatchLevel(ln);
         d_patch_config_util->initialize(level);
      }
   }
   
}

/*************************************************************************
 *
 * Check patch configuration data.
 *
 ************************************************************************/

int PatchConfigurationTester::checkPatchConfiguration(ostream& os) const
{
   int fail_count = 0;

   bool all_levels_correct = true;

   for (int ln = 0; ln < d_hierarchy->getNumberLevels(); ++ln) {

      tbox::pout << "\n\n***************************************************"
                 << "\nChecking patch configuration data for level: ln = " << ln << endl;

      if ( d_patch_config_util->levelIsSet( hier::LevelNumber(ln) ) ) {

         bool level_correct = true;

         tbox::Pointer< hier::PatchLevel<NDIM> > level =
            d_hierarchy->getPatchLevel(ln);

         const hier::ProcessorMapping& mapping = level->getProcessorMapping();
         for (int ip = 0; ip < mapping.getSizeOfMappingArray(); ++ip) {

            bool patch_correct = true;

            if ( mapping.isMappingLocal(ip) ) {

               if ( d_patch_config_util->patchIsSet(hier::PatchNumber(ip),
                                                    hier::LevelNumber(ln)) ) {

                 bool neighbors_good = checkPatchNeighbors(ip, ln); 
                 if (!neighbors_good) {
                    fail_count++;
                    patch_correct = false;
                    tbox::perr << "\tPatch neighbor check FAILED for patch " << ip
                               << " on patch level " << ln << endl;
                 }

                 if ( ln < d_hierarchy->getFinestLevelNumber() ) {
                    bool finer_level_overlap_good = checkFinerLevelPatchOverlap(ip, ln);
                    if (!finer_level_overlap_good) {
                       fail_count++;
                       patch_correct = false;
                       tbox::perr << "\tFiner patch overlap check FAILED for patch " << ip
                                  << " on patch level " << ln << endl;
                    }       
                 }

                 if ( ln > 0 ) {
                    bool coarser_level_overlap_good = checkCoarserLevelPatchOverlap(ip, ln);
                    if (!coarser_level_overlap_good) {
                       fail_count++;
                       patch_correct = false;
                       tbox::perr << "\tCoarser patch overlap check FAILED for patch " << ip
                                  << " on patch level " << ln << endl;
                    }
                 }

               } else {
                  fail_count++;
                  patch_correct = false;
                  tbox::perr << "\tPatch check FAILED: data for patch " << ip
                             << " on patch level " << ln 
                             << " is not set in PatchConfigurationUtilities object." 
                             << endl;
               }
                
            } else {

               if ( d_patch_config_util->patchIsSet(hier::PatchNumber(ip),
                                                    hier::LevelNumber(ln)) ) {
                  fail_count++;
                  patch_correct = false;
                  tbox::perr << "\tPatch check FAILED: data for non-local patch " 
                             << ip << " on patch level " << ln 
                             << " is set in PatchConfigurationUtilities object." 
                             << endl;
               }

            }

            level_correct &= patch_correct;

         }

         if (level_correct) {
            tbox::pout << "\n\tAll data correct on level " << ln << endl;
         }

         all_levels_correct &= level_correct;

      } else {
         all_levels_correct = false; 
         tbox::perr << "\tLevel check FAILED: data for patch level " << ln 
                    << " is not set in PatchConfigurationUtilities object."
                    << endl;
      }

   }

   if (all_levels_correct) {
      tbox::pout << "\nAll data correct on all levels." << endl;
   } 
#if 0
   d_patch_config_util->printClassData(os);
#endif

   return(fail_count);
}

/*************************************************************************
 *
 * Private member functions used to check test results.
 *
 ************************************************************************/

bool PatchConfigurationTester::checkPatchNeighbors(
   int ip, int ln) const
{
   bool data_good = true; 

   tbox::Pointer< hier::PatchLevel<NDIM> > level =
      d_hierarchy->getPatchLevel(ln);
   const int npatches = level->getNumberOfPatches();
   const hier::BoxArray<NDIM>& level_boxes = level->getBoxes();

   hier::BoundaryLookupTable<NDIM>* blut =
      hier::BoundaryLookupTable<NDIM>::getLookupTable();
  
   const hier::Box<NDIM>& pbox = level->getPatch(ip)->getBox();

   for (int codim = NDIM; codim > 0; --codim) {

      const int blut_loc_idx = codim - 1; 
      const int num_locations =
         blut->getMaxLocationIndices()[blut_loc_idx];

      const tbox::Array< hier::PatchConfigurationUtilities<NDIM>::NeighborPatchInfo >&
         util_neighbors = d_patch_config_util->
            getCodimensionNeighborPatchInfo( codim,
                                             hier::PatchNumber(ip),
                                             hier::LevelNumber(ln) );

      tbox::Array< hier::PatchConfigurationUtilities<NDIM>::NeighborPatchInfo >
         my_neighbors(num_locations);

      int num_my_neighbors = 0;

      tbox::List< hier::IntVector<NDIM> >::Iterator
         is( level->getShiftsForPatch(ip) );
      bool zero_shift = true;

      while (is || zero_shift) {

         hier::IntVector<NDIM> nshift(0);
         if (!zero_shift) {
            nshift = -is();
         }

         for (int loc = 0; loc < num_locations; ++loc) {

            const tbox::Array<int>& dirs = blut->getDirections(loc, codim);

            hier::Box<NDIM> bregion( pbox );
            for (int i = 0; i < dirs.size(); ++i) {
               if ( blut->isUpper(loc, codim, i) ) {
                  bregion.lower(dirs[i]) = pbox.upper(dirs[i]) + 1;
                  bregion.upper(dirs[i]) = pbox.upper(dirs[i]) + 1;
               } else {
                  bregion.lower(dirs[i]) = pbox.lower(dirs[i]) - 1;
                  bregion.upper(dirs[i]) = pbox.lower(dirs[i]) - 1;
               }
            }

            for (int ni = 0; ni < npatches; ++ni) {

               if ( !zero_shift || (ni != ip) ) {

                  const hier::Box<NDIM> shifted_nbox(
                     hier::Box<NDIM>::shift(level_boxes.getBox(ni), nshift) );

                  const hier::Box<NDIM> intersection = shifted_nbox * bregion;                       

                  if ( !intersection.empty() ) {

                     my_neighbors[num_my_neighbors].d_neighbor_patch_number = ni;
                     my_neighbors[num_my_neighbors].d_neighbor_type = codim;
                     my_neighbors[num_my_neighbors].d_location_index = loc;
                     my_neighbors[num_my_neighbors].d_neighbor_shift = nshift;

                     num_my_neighbors++;

                  }  // add neighbor to my array of neighbors

               }  // patch may only be a neighbor of itself when shift is non-zero
         
            }  // iterate over all level patches to check for neighbors

         }  // iterate over locations for codimension

         if (!zero_shift) {
            is++;
         } else {
            zero_shift = false;
         }

      }  // iterate over valid shifts for codimension

      my_neighbors.resizeArray(num_my_neighbors);

      data_good &= comparePatchNeighbors(ip, ln, codim, util_neighbors, my_neighbors); 

   }  // iterate over potential neighbor codimensions

   return(data_good);
}

bool PatchConfigurationTester::checkFinerLevelPatchOverlap(
   int ip, int ln) const
{
   bool data_good = true;

   if ( ln < d_hierarchy->getFinestLevelNumber() ) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level =
         d_hierarchy->getPatchLevel(ln); 
      tbox::Pointer< hier::PatchLevel<NDIM> > finer_level =
         d_hierarchy->getPatchLevel(ln+1);

      hier::BoxArray<NDIM> cf_boxes = finer_level->getBoxes();
      cf_boxes.coarsen( finer_level->getRatioToCoarserLevel() );

      tbox::Array<int> fpatch_ids(cf_boxes.size());

      const hier::Box<NDIM>& pbox = level->getPatch(ip)->getBox();

      int pcount = 0; 
      for (int ifp = 0; ifp < cf_boxes.size(); ++ifp) {
         if ( pbox.intersects(cf_boxes.getBox(ifp)) ) {
            fpatch_ids[pcount] = ifp;
            pcount++;
         }
      }
      fpatch_ids.resizeArray(pcount);

      const tbox::Array<int>& fpatch_ids_util = 
         d_patch_config_util->getFinerLevelOverlapPatchIndices(
            hier::PatchNumber(ip), hier::LevelNumber(ln) ); 
    
      if ( fpatch_ids.size() != fpatch_ids_util.size() ) {
         data_good = false; 
         tbox::perr << "\tFiner patch overlap check FAILED: patch "
                    << ip << " on level " << ln
                    << "\n number of fine patches does not match."
                    << "\n    PatchConfigUtils gives " 
                    << fpatch_ids_util.size() << " fine overlap patches."
                    << "\n    Test computes " 
                    << fpatch_ids.size() << " fine overlap patches."
                    << endl; 
      }
 
      if (data_good) {
         bool ids_match = true;
         int ifp = 0;
         while ( ids_match && (ifp < fpatch_ids.size()) ) { 
            ids_match = (fpatch_ids[ifp] == fpatch_ids_util[ifp]);
            ifp++; 
         }
         if (!ids_match) {
            data_good = false;
            tbox::perr << "\tFiner patch overlap check FAILED: patch "
                       << ip << " on level " << ln
                       << "\n fine patch overlap ids incorrect. Check log file."
                       << endl;
            tbox::plog << "Fine overlap patch ids for patch " 
                       << ip << " on level " << ln << "..." << endl;
            tbox::plog << "   Test computes fine patch ids = ...." << endl;
            for (ifp = 0; ifp < fpatch_ids.size() - 1; ++ifp) {
               tbox::plog << fpatch_ids[ifp] << " , ";
            }
            if ( ifp < fpatch_ids.size() ) {
               tbox::plog << fpatch_ids[ifp] << endl;
            }
            tbox::plog << "   PatchConfigUtils gives fine patch ids = ...." << endl;
            for (ifp = 0; ifp < fpatch_ids_util.size() - 1; ++ifp) {
               tbox::plog << fpatch_ids_util[ifp] << " , ";
            }
            if ( ifp < fpatch_ids_util.size() ) {
               tbox::plog << fpatch_ids_util[ifp] << endl;
            }
         }
      }
       
   }

   return(data_good);
}

bool PatchConfigurationTester::checkCoarserLevelPatchOverlap(
   int ip, int ln) const
{
   bool data_good = true;

   if ( ln > 0 ) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level =
         d_hierarchy->getPatchLevel(ln); 
      tbox::Pointer< hier::PatchLevel<NDIM> > coarser_level =
         d_hierarchy->getPatchLevel(ln-1);

      hier::BoxArray<NDIM> fc_boxes = coarser_level->getBoxes();
      fc_boxes.refine( level->getRatioToCoarserLevel() );

      tbox::Array<int> cpatch_ids(fc_boxes.size());

      const hier::Box<NDIM>& pbox = level->getPatch(ip)->getBox();

      int pcount = 0; 
      for (int icp = 0; icp < fc_boxes.size(); ++icp) {
         if ( pbox.intersects(fc_boxes.getBox(icp)) ) {
            cpatch_ids[pcount] = icp;
            pcount++;
         }
      }
      cpatch_ids.resizeArray(pcount);

      const tbox::Array<int>& cpatch_ids_util = 
         d_patch_config_util->getCoarserLevelOverlapPatchIndices(
            hier::PatchNumber(ip), hier::LevelNumber(ln) ); 
    
      if ( cpatch_ids.size() != cpatch_ids_util.size() ) {
         data_good = false; 
         tbox::perr << "\tCoarser patch overlap check FAILED: patch "
                    << ip << " on level " << ln
                    << "\n number of coarse patches does not match."
                    << "\n    PatchConfigUtils gives " 
                    << cpatch_ids_util.size() << " coarse overlap patches."
                    << "\n    Test computes " 
                    << cpatch_ids.size() << " coarse overlap patches."
                    << endl;
      }
 
      if (data_good) {
         bool ids_match = true;
         int icp = 0;
         while ( ids_match && (icp < cpatch_ids.size()) ) { 
            ids_match = (cpatch_ids[icp] == cpatch_ids_util[icp]);
            icp++; 
         }
         if (!ids_match) {
            data_good = false;
            tbox::perr << "\tFiner patch overlap check FAILED: patch "
                       << ip << " on level " << ln
                       << "\n fine patch overlap ids incorrect. Check log file."
                       << endl;
            tbox::plog << "Fine overlap patch ids for patch " 
                       << ip << " on level " << ln << "..." << endl;
            tbox::plog << "   Test computes fine patch ids = ...." << endl;
            for (icp = 0; icp < cpatch_ids.size() - 1; ++icp) {
               tbox::plog << cpatch_ids[icp] << " , ";
            }
            if ( icp < cpatch_ids.size() ) {
               tbox::plog << cpatch_ids[icp] << endl;
            }
            tbox::plog << "   PatchConfigUtils gives fine patch ids = ...." << endl;
            for (icp = 0; icp < cpatch_ids_util.size() - 1; ++icp) {
               tbox::plog << cpatch_ids_util[icp] << " , ";
            }
            if ( icp < cpatch_ids_util.size() ) {
               tbox::plog << cpatch_ids_util[icp] << endl;
            }
         }
      }
       
   }

   return(data_good);
}

bool PatchConfigurationTester::comparePatchNeighbors(
   int ip, 
   int ln,
   int codim,
   const tbox::Array< hier::PatchConfigurationUtilities<NDIM>::NeighborPatchInfo >&
      util_neighbors,
   const tbox::Array< hier::PatchConfigurationUtilities<NDIM>::NeighborPatchInfo >&
      test_neighbors) const
{
   bool data_good = true;
 
   if (util_neighbors.size() != test_neighbors.size()) {
      data_good = false;
      tbox::perr << "\tPatch neighbor check FAILED: patch "
                 << ip << " on level " << ln
                 << "\n number of neighbors of codim " << codim << "does not match."
                 << "\n    PatchConfigUtils gives " 
                 << util_neighbors.size() << " neighbors."
                 << "\n    Test computes " 
                 << test_neighbors.size() << " neighbors."
                 << endl;
   }

   if (data_good) {
      bool neighbors_match = true;
      tbox::Array<bool> good_neighbors(util_neighbors.size());
      for (int ni = 0; ni < util_neighbors.size(); ++ni) {
         good_neighbors[ni] = 
            (util_neighbors[ni].d_neighbor_patch_number ==
             test_neighbors[ni].d_neighbor_patch_number);
         good_neighbors[ni] &= 
            (util_neighbors[ni].d_neighbor_type ==
             test_neighbors[ni].d_neighbor_type);
         good_neighbors[ni] &= 
            (util_neighbors[ni].d_location_index ==
             test_neighbors[ni].d_location_index);
         good_neighbors[ni] &= 
            (util_neighbors[ni].d_neighbor_shift ==
             test_neighbors[ni].d_neighbor_shift);
         neighbors_match &= good_neighbors[ni];
      }
      if (!neighbors_match) {
         data_good = false;
         tbox::perr << "\tPatch neighbor check FAILED: patch "
                    << ip << " on level " << ln 
                    << " for codim " << codim << endl;
         for (int ni = 0; ni < util_neighbors.size(); ++ni) {
            if (!good_neighbors[ni]) { 
               tbox::plog << "  Bad patch neighbor info at neighbor index " << ni << endl;
               tbox::plog << "   Neighbor patch number: Util -> " 
                          << util_neighbors[ni].d_neighbor_patch_number
                          << " , Test -> "
                          << test_neighbors[ni].d_neighbor_patch_number << endl;
               tbox::plog << "   Neighbor type: Util -> " 
                          << util_neighbors[ni].d_neighbor_type
                          << " , Test -> "
                          << test_neighbors[ni].d_neighbor_type << endl;
               tbox::plog << "   Neighbor location: Util -> " 
                          << util_neighbors[ni].d_location_index
                          << " , Test -> "
                          << test_neighbors[ni].d_location_index << endl;
               tbox::plog << "   Neighbor shift: Util -> " 
                          << util_neighbors[ni].d_neighbor_shift
                          << " , Test -> "
                          << test_neighbors[ni].d_neighbor_shift << endl;
             }
          }
      }
   }

   return(data_good);
}

/*************************************************************************
 *
 * Print patch hierarchy information to output stream.
 *
 ************************************************************************/

void PatchConfigurationTester::printHierarchyData(
   ostream& os) const
{
   for (int ln = 0; ln < d_hierarchy->getNumberLevels(); ln++) {
      tbox::Pointer< hier::PatchLevel<NDIM> > level = 
         d_hierarchy->getPatchLevel(ln);

      os << "\n Patches on level: " << level->getLevelNumber() << endl;

      const hier::BoxArray<NDIM>& level_boxes = level->getBoxes();
      for (int ib = 0; ib < level_boxes.getNumberOfBoxes(); ib++) {
         os << "   Patch[" << ib << "] = " 
            << level_boxes.getBox(ib) << endl;
      }
   }
}

/*************************************************************************
 *
 * Read gridding parameter data from input file.
 *
 ************************************************************************/

void PatchConfigurationTester::getGriddingParametersFromInput(
   tbox::Pointer<tbox::InputDatabase> gridding_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!gridding_db.isNull());
#endif
 
   int ln;
 
   d_num_levels = 
      gridding_db->getIntegerWithDefault("num_levels", d_num_levels);
   if (d_num_levels < 1) {
      TBOX_ERROR("Input error for " << d_object_name << ":  "
                 << "Key data 'num_levels' found in input is < 1" << endl);
   }
   if ( (d_test_to_run == "COARSEN_TEST") && 
        d_num_levels > 2) {
      d_num_levels = 2;
      TBOX_WARNING("Setting 'num_levels' to 2 for COARSEN_TEST" << endl);
   }

   d_ratio_to_coarser.resizeArray(d_num_levels);
   d_ratio_to_coarser[0] = hier::IntVector<NDIM>(1);

   tbox::Pointer<tbox::Database> ratio_to_coarser_db;
   if ( d_num_levels > 1) {

      if (!gridding_db->keyExists("ratio_to_coarser")) {
         TBOX_ERROR("Input error for " << d_object_name << ":  "
                    << "Key data `ratio_to_coarser' not found in input.");
      } else {
         ratio_to_coarser_db = gridding_db->getDatabase("ratio_to_coarser");
         for (ln = 1; ln < d_num_levels; ln++) {
            char level_name[32];
            sprintf(level_name, "level_%d", ln);
 
            if (!ratio_to_coarser_db->keyExists(level_name)) {
                TBOX_ERROR("Input error for " << d_object_name << ":  "
                           <<"Key data `" << level_name
                           << "' not found in ratio_to_coarser input" << endl);
            }
 
            int* temp_ratio_to_coarser = d_ratio_to_coarser[ln];
            ratio_to_coarser_db->getIntegerArray(level_name,
                                                 temp_ratio_to_coarser, NDIM);
            if (d_ratio_to_coarser[ln].max() < 1) {
               TBOX_ERROR("Input error for " << d_object_name << ":  "
                          <<"Key data `" << level_name
                          << "' in ratio_to_coarser input has entry < 1" << endl);
            }
         }
      }

   }

   int* temp_npatches_on_coarsest = d_npatches_on_coarsest;
   gridding_db->getIntegerArray("npatches_on_coarsest",
                                 temp_npatches_on_coarsest, NDIM);
   if (d_npatches_on_coarsest.min() < 1) {
      TBOX_ERROR("Input error for " << d_object_name << ":  "
                 << "npatches_on_coarsest input has entries < 1" << endl);
   }

}
