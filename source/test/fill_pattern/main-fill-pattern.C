//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/test/fill_pattern/main-fill_pattern.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2901 $
// Modified:    $LastChangedDate: 2009-02-12 15:21:53 -0800 (Thu, 12 Feb 2009) $// Description: Main program for testing SAMRAI VariableFillPatterns
//


#include "SAMRAI_config.h"

#include "SkeletonGridGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "NodeData.h"
#include "NodeVariable.h"
#include "FirstLayerCellFillPattern.h"
#include "FirstLayerCellNoCornersFillPattern.h"
#include "SecondLayerNodeFillPattern.h"
#include "SecondLayerNodeNoCornersFillPattern.h"
#include "RefineAlgorithm.h"
#include "BoxArray.h"
#include "PatchHierarchy.h"
#include "tbox/SAMRAIManager.h"

#include <cstring>
#include <vector>

using namespace std;
using namespace SAMRAI;



/*
 *
 * Test program to test VariableFillPattern implementations
 *
 */

void txt2boxes(const char* txt, hier::BoxArray<NDIM>& boxes)
{
   // algorithm:
   // find width
   // find height
   // find x locations in i,j
   // foreach x1 in x:
   //   find x' where i>i1, j>j1
   //     foreach x2 in x':
   //       if 4 corners && interior cells blank, (x1,x2) is a box
   // translate coordinates into cell-centered, lower-left origin based

   int width = -1;
   for (unsigned int idx = 0; idx < strlen(txt)-1; idx++) {
      if ( ('x' == txt[idx]  || '.' == txt[idx]) &&
           ('.' == txt[idx+1] || '|' == txt[idx+1]) ) {
         width = idx+1;
         break;
      }
   }
   if (-1 == width) {
      cout << "error in box txt" << endl;
      exit(1);
   }

   int height = strlen(txt)/width;

   // Find cell height
   int cell_height = (height-1)/2;
   int cell_max = cell_height-1;

   // make vector of x locations
//   vector<pair<int,int> > ix;
   tbox::List<hier::IntVector<NDIM> > ix;
   for (unsigned int idx = 0; idx < strlen(txt); idx++) {
      if ('x' == txt[idx]) {
         int j = idx / width;
         int i = idx - j*width;
         hier::IntVector<NDIM> pt;
         pt(0) = i;
         pt(1) = j; 
//         ix.push_back(pair<int,int>(i,j));
         ix.appendItem(pt);
      }
   }

   // foreach x1 in x
//   vector< pair<int,int> >::iterator it;
//   for (it = ix.begin(); it != ix.end(); it++) {
   for (tbox::List<hier::IntVector<NDIM> >::Iterator it(ix); it; it++) {

//      vector< pair<int,int> >::iterator it2;

      // We need to gather all potential boxes rooted here, and then
      // only take the smallest one.

//      vector< hier::Box<NDIM> > boxes_here;
//      boxes_here.clear();
      hier::BoxList<NDIM> boxes_here;

//      for (it2 = ix.begin(); it2 != ix.end(); it2++) {
      for (tbox::List<hier::IntVector<NDIM> >::Iterator it2(ix); it2; it2++) {

//         if (it2->first > it->first &&
//             it2->second > it->second) {
         if (it2()(0) > it()(0) &&
             it2()(1) > it()(1)) {

            bool isbox = true;

            // If the two other corners exist, and...
//            int i1 = it->first;
//            int j1 = it2->second;
            int i1 = it()(0);
            int j1 = it2()(1);
            int idx1 = j1*width +i1;
            if (txt[idx1] != 'x') isbox = false;

//            int i2 = it2->first;
//            int j2 = it->second;
            int i2 = it2()(0);
            int j2 = it()(1);
            int idx2 = j2*width +i2;
            if (txt[idx2] != 'x') isbox = false;

            // ...interior cells contain no corners
//            for (int i = it->first+1; i < it2->first; i++) {
//               for (int j = it->second+1; j < it2->second; j++) {
            for (int i = it()(0)+1; i < it2()(0); i++) {
               for (int j = it()(1)+1; j < it2()(1); j++) {
                  int idx = j*width +i;
                  if ('x' == txt[idx]) isbox = false;
                  if ('-' == txt[idx]) isbox = false;
                  if ('|' == txt[idx]) isbox = false;
               }
            }

            if (isbox) {

               // Translate indices into node centered coords
//               int i0 = it->first/4;
//               int i1 = it2->first/4;
//               int j0 = it->second/2;
//               int j1 = it2->second/2;
               int i0 = it()(0)/4;
               i1 = it2()(0)/4;
               int j0 = it()(1)/2;
               j1 = it2()(1)/2;

               i1--;
               j1--;

               // Flip coordinates vertically.
               j0 = cell_max-j0;
               j1 = cell_max-j1;

               // Lower left uses j1, upper right j0
               int tmp = j1;
               j1 = j0;
               j0 = tmp;

               hier::Box<NDIM> abox(hier::Index<NDIM>(i0,j0),
                                    hier::Index<NDIM>(i1,j1));
//               boxes_here.push_back(abox);
               boxes_here.appendItem(abox);
            }
         }
      }

      // Find smallest box at this 'x'
      if (boxes_here.size()) {

//         hier::Box<NDIM> smallest_box(boxes_here[0]);
         hier::Box<NDIM> smallest_box(boxes_here.getFirstItem());

//         for (vector< hier::Box<NDIM> >::iterator it = boxes_here.begin();
//              it != boxes_here.end(); it++) {
         for (hier::BoxList<NDIM>::Iterator itb(boxes_here); itb; itb++) {
//            if ((*it).numberCells() < smallest_box.numberCells()) {
//               smallest_box = *it;
            if (itb().numberCells() < smallest_box.numberCells()) {
               smallest_box = itb();
            }
         }

         boxes.resizeBoxArray(boxes.getNumberOfBoxes()+1);
         int bidx = boxes.getNumberOfBoxes()-1;
         boxes[bidx] = smallest_box;
      }

   }

   // Shift all boxes into SAMRAI coordinates
   for (int idx = 0; idx < boxes.getNumberOfBoxes(); idx++) {
      boxes[idx].shift(-hier::IntVector<NDIM>(2));
   }
}

int txt_width(const char* txt)
{
   int width = -1;
   for (unsigned int idx = 0; idx < strlen(txt)-1; idx++) {
      if ( ('x' == txt[idx]  || '.' == txt[idx]) &&
           ('.' == txt[idx+1] || '|' == txt[idx+1]) ) {
         width = idx+1;
         break;
      }
   }
   if (-1 == width) {
      cout << "error in box txt" << endl;
      exit(1);
   }
   return width;
}

bool txt_next_val(const char* txt,
                  int& idx,
                  const hier::PatchData<NDIM>& data,
                  int* datapt,
                  bool is_node)
{
   // Find text size
   int txt_w = txt_width(txt);
   int txt_h = strlen(txt)/txt_w;

   // Find grid size
   int grid_height = (txt_h-1)/2;
   int grid_width = (txt_w-1)/4;
   int grid_max = grid_height-1;

   int cnt_max = 10000; // limit infinite loop possibility
   int cnt = 0;
   do {

      //
      // Translate domain local idx into grid idx
      //

      //const hier::Box& ghost_box = data.getGhostBox();
      hier::Box<NDIM> ghost_box;
      if (is_node) {
         ghost_box = pdat::NodeGeometry<NDIM>::toNodeBox(data.getGhostBox());
      } else {
         ghost_box = data.getGhostBox();
      }
      // Translate domain idx to domain coordinates
      int domain_i = idx % ghost_box.numberCells(0);
      int domain_j = idx / ghost_box.numberCells(0);

      hier::Box<NDIM> shifted_box(hier::Box<NDIM>::shift(ghost_box,
                                      hier::IntVector<NDIM>(2)));
      // Translate domain coordinates into grid zone coordintes
      int di = shifted_box.lower()(0);
      int dj = shifted_box.lower()(1);

      int grid_i = domain_i +di;
      int grid_j = domain_j +dj;

      // If we outside the grid, there cannot be a value here
      if (grid_i < 0 || grid_j < 0) {
         idx++;
         continue;
      }
      if (grid_i > grid_width || grid_j > grid_height) {
         idx++;
         continue;
      }
      // Translate grid coords to text coordinates.  Text coordinates
      // have j increasing downwards.
      int txt_zone_i = grid_i*4+2;
      int txt_zone_j = (grid_max-grid_j)*2+1;

      int txt_node_i = grid_i*4;
      int txt_node_j = (grid_max-grid_j)*2+2;

      // Translate text coordinates to txt idx
      unsigned int txt_zone_idx = txt_zone_i+txt_zone_j*txt_w;
      int txt_node_idx = txt_node_i+txt_node_j*txt_w;

      // If we're past the end of the txt, return false
      if (txt_zone_idx > strlen(txt)) {
         return false;
      }

      // Check for non-zero zone data
      if (' ' != txt[txt_zone_idx]) {

         istringstream valstr(&txt[txt_zone_idx]);
         valstr >> *datapt;
         return true;
      }

      // Check for numeric node data
      if ('0' == txt[txt_node_idx] ||
          '1' == txt[txt_node_idx] ||
          '2' == txt[txt_node_idx] ||
          '3' == txt[txt_node_idx] ||
          '4' == txt[txt_node_idx] ||
          '5' == txt[txt_node_idx] ||
          '6' == txt[txt_node_idx] ||
          '7' == txt[txt_node_idx] ||
          '8' == txt[txt_node_idx] ||
          '9' == txt[txt_node_idx]) {

         istringstream valstr(&txt[txt_node_idx]);
         valstr >> *datapt;
         return true;
      }

      idx++; // advance to next domain idx

   } while ( cnt++ < cnt_max );

   cout << "Data reading loop exceeded maximum iterations"
        << __LINE__ << " in "
        << __FILE__ << endl;

   exit(1);
}


void txt2data(const char* txt,
              const hier::PatchData<NDIM>& data,
              int* datptr,
              bool zero_out,
              bool is_node)
{
   if (zero_out) memset(datptr, 0., data.getGhostBox().size()*sizeof(int));

   int idx = 0;
   int datapt;

   while ( txt_next_val(txt, idx, data, &datapt, is_node) ) {
      datptr[idx++] = datapt;
   }
}


/*
  Acceptance test cases.  First iteration at an executable
  specification of the desired change.
*/

bool SingleLevelTestCase(const char* levelboxes_txt,
                         const char* initialdata_txt[],
                         const char* finaldata_txt[],
                         tbox::Pointer<hier::Variable<NDIM> > variable, 
                         tbox::Pointer<xfer::VariableFillPattern<NDIM> > fill_pattern)
{
   const std::string& pattern_name = fill_pattern->getPatternName();

   hier::BoxArray<NDIM> level_boxes;
   txt2boxes(levelboxes_txt, level_boxes);

   hier::Box<NDIM> domain_box;
   for (int i = 0; i < level_boxes.getNumberOfBoxes(); i++) {
      domain_box += level_boxes[i];
   }

   hier::BoxArray<NDIM> physical_domain(domain_box);

   tbox::Pointer<geom::SkeletonGridGeometry<NDIM> > geom =
      new geom::SkeletonGridGeometry<NDIM>("GridGeometry", level_boxes);

   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = new
      hier::PatchHierarchy<NDIM>("hier", geom);

//   hier::MappedBoxHierarchy &lh = hierarchy->getMappedBoxHierarchy();


//   lh.setMappedBoxLevelParameters(0,
//                                  hier::IntVector<NDIM>(1),
//                                  hier::IntVector<NDIM>(1),
//                                  hier::IntVector<NDIM>(1));

//   hier::MappedBoxLevel mblevel( NULL, hier::IntVector<NDIM>(1) );


   const int num_nodes = tbox::SAMRAI_MPI::getNodes();
   const int num_boxes = level_boxes.size();
//   hier::GlobalId::LocalIndex local_id = 0;
   tbox::Array<int> mapping(num_boxes);
   for ( int i=0; i<num_boxes; ++i ) {

      if (i < num_boxes / num_nodes) {
         mapping[i] = 0;
      } else {
         mapping[i] = 1;
      }

   }

   hier::ProcessorMapping proc_mapping(mapping);

   int level_no = 0;
   hierarchy->makeNewPatchLevel(level_no, hier::IntVector<NDIM>(1),
                                level_boxes, proc_mapping);

   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(0);

   // There is one variable-context pair with a gcw of 2
   
   xfer::RefineAlgorithm<NDIM> refine_alg;

   tbox::Pointer<hier::VariableContext> context =
      hier::VariableDatabase<NDIM>::getDatabase()->getContext("CONTEXT");

   hier::IntVector<NDIM> ghost_cell_width(2);

   int data_id =
      hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
         variable, context, ghost_cell_width);

   refine_alg.registerRefine(data_id, data_id, data_id, NULL, fill_pattern);

   level->allocatePatchData(data_id);

   TBOX_ASSERT(tbox::SAMRAI_MPI::getNodes() <= 2);

   if (pattern_name == "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN" ||
       pattern_name == "FIRST_LAYER_CELL_FILL_PATTERN") {
      // Loop over each patch and initialize data
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
         tbox::Pointer<hier::Patch<NDIM> > patch(level->getPatch(p()));  
         tbox::Pointer<pdat::CellData<NDIM,int> > cdata  =
            patch->getPatchData(data_id);

         int data_txt_id = p();
//         if (tbox::SAMRAI_MPI::getRank() == 1) {
//            data_txt_id += (num_boxes / num_nodes);
//         }

         txt2data(initialdata_txt[data_txt_id], *cdata,
                  cdata->getPointer(), false, false);
      }
   } else if (pattern_name == "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN" ||
              pattern_name == "SECOND_LAYER_NODE_FILL_PATTERN") {
      // Loop over each patch and initialize data
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
         tbox::Pointer<hier::Patch<NDIM> > patch(level->getPatch(p()));
         tbox::Pointer<pdat::NodeData<NDIM,int> > ndata  =
            patch->getPatchData(data_id);

         int data_txt_id = p();
//         if (tbox::SAMRAI_MPI::getRank() == 1) {
//            data_txt_id += (num_boxes / num_nodes);
//         }

         txt2data(initialdata_txt[data_txt_id], *ndata,
                  ndata->getPointer(), false, true);
      }
   }
/*
   hier::OverlapConnectorUtil connect_util(dim);

   hier::Connector peer_cnect;
   peer_cnect.initialize(*(level->getMappedBoxLevel()),
                         *(level->getMappedBoxLevel()),
                         hier::IntVector<NDIM>(2));
  
   connect_util.findEdges(peer_cnect);
*/
   // Create and run comm schedule
   refine_alg.createSchedule(level)->fillData(0.0, false);

   // Check for expected data
   bool failed = false;


   if (pattern_name == "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN" ||
       pattern_name == "FIRST_LAYER_CELL_FILL_PATTERN") {
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
         tbox::Pointer<hier::Patch<NDIM> > patch(level->getPatch(p()));
         tbox::Pointer<pdat::CellData<NDIM,int> > cdata  =
            patch->getPatchData(data_id);

         pdat::CellData<NDIM,int> expected(cdata->getBox(),
                                              cdata->getDepth(),
                                              ghost_cell_width);

         int data_txt_id = p();
//         if (tbox::SAMRAI_MPI::getRank() == 1) {
//            data_txt_id += (num_boxes / num_nodes);
//         }

         txt2data(finaldata_txt[data_txt_id],
                  expected, expected.getPointer(), false, false);

         for (pdat::CellData<NDIM,int>::Iterator ci(cdata->getGhostBox());
              ci; ci++) {
            if ((*cdata)(ci()) != expected(ci())) {
               failed = true;
            }
         }

      }
   } else if (pattern_name == "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN" ||
              pattern_name == "SECOND_LAYER_NODE_FILL_PATTERN") {
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
         tbox::Pointer<hier::Patch<NDIM> > patch(level->getPatch(p()));
         tbox::Pointer<pdat::NodeData<NDIM,int> > ndata  =
            patch->getPatchData(data_id);

         pdat::NodeData<NDIM,int> expected(ndata->getBox(),
                                              ndata->getDepth(),
                                              ghost_cell_width);

         int data_txt_id = p();
//         if (tbox::SAMRAI_MPI::getRank() == 1) {
//            data_txt_id += (num_boxes / num_nodes);
//         }

         txt2data(finaldata_txt[data_txt_id],
                  expected, expected.getPointer(), false, true);

         for (pdat::NodeData<NDIM,int>::Iterator ni(ndata->getGhostBox());
              ni; ni++) {
            if ((*ndata)(ni()) != expected(ni())) {
               failed = true;
            }
         }

      }
   }

   if (failed) {
      tbox::perr << "FAILED: - Test of "<< pattern_name << endl;
   }

   return failed; 
}


/*
  This tests FirstLayerCellNoCornersFillPattern ..
*/

bool Test_FirstLayerCellNoCornersVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 4
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 2
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 0 data before comm
   
   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 2 . 2 . 0 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // expected patch 3 data after comm
   
   const char* final3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 3 . 3 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // expected patch 0 data after comm
   
   const char* final0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 0 . 0 . 2 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 1 . 1 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;


   const char* initial_txt[4] = {initial0_txt, initial1_txt,
                                 initial2_txt, initial3_txt};
   const char* final_txt[4] = {final0_txt, final1_txt,
                               final2_txt, final3_txt};

   tbox::Pointer<pdat::CellVariable<NDIM,int> > var =
      new pdat::CellVariable<NDIM,int>("1cellnocorners");

   tbox::Pointer<pdat::FirstLayerCellNoCornersFillPattern<NDIM> > fill_pattern =
      new pdat::FirstLayerCellNoCornersFillPattern<NDIM>();

   return SingleLevelTestCase(levelboxes_txt,
                              initial_txt,
                              final_txt,
                              var,
                              fill_pattern);
}

/*
  This tests FirstLayerCellFillPattern
*/

bool Test_FirstLayerCellVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 4
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 2
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;


   const char* final0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 2 . 2 . 3 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;


   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 2 . 3 . 3 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;


   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 0 . 0 . 1 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;


   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 0 . 1 . 1 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
      ;


   const char* initial_txt[4] = {initial0_txt, initial1_txt,
                                 initial2_txt, initial3_txt};
   const char* final_txt[4] = {final0_txt, final1_txt,
                               final2_txt, final3_txt};

   tbox::Pointer<pdat::CellVariable<NDIM,int> > var =
      new pdat::CellVariable<NDIM,int>("1cell");

   tbox::Pointer<pdat::FirstLayerCellFillPattern<NDIM> > fill_pattern =
      new pdat::FirstLayerCellFillPattern<NDIM>();

   return SingleLevelTestCase(levelboxes_txt,
                              initial_txt,
                              final_txt,
                              var,
                              fill_pattern);
}



bool Test_SecondLayerNodeNoCornersVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;


   const char* final0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 2 . 2 . 2 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 1 . 3 . 3 . 3 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;


   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 0 . 0 . 0 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;


   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 3 . 1 . 1 . 1 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   const char* initial_txt[4] = {initial0_txt, initial1_txt,
                                 initial2_txt, initial3_txt};
   const char* final_txt[4] = {final0_txt, final1_txt,
                               final2_txt, final3_txt};

   tbox::Pointer<pdat::NodeVariable<NDIM,int> > var =
      new pdat::NodeVariable<NDIM,int>("secondnodenocorners");

   tbox::Pointer<pdat::SecondLayerNodeNoCornersFillPattern<NDIM> >
   fill_pattern =
      new pdat::SecondLayerNodeNoCornersFillPattern<NDIM>();

   return SingleLevelTestCase(levelboxes_txt,
                              initial_txt,
                              final_txt,
                              var,
                              fill_pattern);
}

bool Test_SecondLayerNodeVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;


   const char* final0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 2 . 2 . 2 . 3 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 2 . 3 . 3 . 3 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;


   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 0 . 0 . 0 . 1 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;


   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 0 . 1 . 1 . 1 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
      ;

   const char* initial_txt[4] = {initial0_txt, initial1_txt,
                                 initial2_txt, initial3_txt};
   const char* final_txt[4] = {final0_txt, final1_txt,
                               final2_txt, final3_txt};

   tbox::Pointer<pdat::NodeVariable<NDIM,int> > var =
      new pdat::NodeVariable<NDIM,int>("secondnode");

   tbox::Pointer<pdat::SecondLayerNodeFillPattern<NDIM> >
   fill_pattern =
      new pdat::SecondLayerNodeFillPattern<NDIM>();

   return SingleLevelTestCase(levelboxes_txt,
                              initial_txt,
                              final_txt,
                              var,
                              fill_pattern);
}

int main(int argc, char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc,&argv);
   tbox::SAMRAIManager::startup();

   int failures = 0;
   
   failures += Test_FirstLayerCellNoCornersVariableFillPattern();
   failures += Test_FirstLayerCellVariableFillPattern();
   failures += Test_SecondLayerNodeNoCornersVariableFillPattern();
   failures += Test_SecondLayerNodeVariableFillPattern();

   if (failures == 0) { 
      tbox::pout << "\nPASSED:  fill_pattern" << endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
   
   return(failures);
}
