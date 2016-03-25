//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/dataaccess/main.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Main program for testing SAMRAI data access
//

#include "SAMRAI_config.h"

#include "Box.h"
#include "CellData.h"
#include "EdgeData.h"
#include "FaceData.h"
#include "NodeData.h"
#include "SideData.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include "tbox/SAMRAIManager.h"

using namespace std;
using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   int error_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("iteratortest.log");

      /*
       * Regular pointer tests.
       */

      hier::Index<NDIM> box_lower(0);
      hier::Index<NDIM> box_upper;

      for (int d = 0; d < NDIM; d++) {
	 box_upper(d) = (d+4)*3;
      }

      hier::Box<NDIM> box(box_lower, box_upper);

      pdat::CellData<NDIM,double> cell_data(box, 1, hier::IntVector<NDIM>(0));
      pdat::FaceData<NDIM,double> face_data(box, 1, hier::IntVector<NDIM>(0));
      pdat::NodeData<NDIM,double> node_data(box, 1, hier::IntVector<NDIM>(0));
      pdat::EdgeData<NDIM,double> edge_data(box, 1, hier::IntVector<NDIM>(0));
      pdat::SideData<NDIM,double> side_data(box, 1, hier::IntVector<NDIM>(0));

      /*
       * The tests of the iterators first fill the patch data by directly
       * accessing the pointer to the double data.  The iterators should access
       * the data in exactly the same order that it was filled.
       */

      /*
       * The tests of the index classes loop over the box with a CellIterator.
       * The CellIterator gives each cell-centered index within the box.
       * For each cell-centered index, the tests loop over every possible
       * set of parameters for the datatype-specific index class (EdgeIndex,
       * SideIndex, etc.).  Then the test checks that accessing the data using
       * the index class retrieves data from the correct point in the data array.
       */ 


      /*
       * Test CellIterator
       */

      double* cell_ptr = cell_data.getPointer();

      int cell_data_size = box.size();

      for (int i = 0; i < cell_data_size; i++) {
	 cell_ptr[i] = (double) i;
      }

      int j = 0;
      for (pdat::CellIterator<NDIM> ci(box); ci; ci++) {
	 if (!tbox::MathUtilities<double>::equalEps(cell_data(ci()), cell_ptr[j])) {
	    tbox::perr << "FAILED: - CellIterator test" << endl;
	    ++error_count;
	 }
	 j++;
      }

      /*
       * Test FaceIterator
       */

      double* face_ptr[NDIM];

      for (int axis = 0; axis < NDIM; axis++) {
 
	 face_ptr[axis] = face_data.getPointer(axis);

	 hier::Box<NDIM> face_box = pdat::FaceGeometry<NDIM>::toFaceBox(box,
									axis);

	 int face_data_size = face_box.size();

	 for (int i = 0; i < face_data_size; i++) {
	    face_ptr[axis][i] = (double) i;
	 }

	 j = 0;
	 for (pdat::FaceIterator<NDIM> fi(box, axis); fi; fi++) {
	    if (!tbox::MathUtilities<double>::equalEps(face_data(fi()), face_ptr[axis][j])) {
	       tbox::perr << "FAILED: - FaceIterator test" << endl;
	       ++error_count;
	    }
	    j++;
	 } 
      }

      /*
       * Test FaceIndex
       */

      for (pdat::CellIterator<NDIM> ifc(box); ifc; ifc++) {

	 for (int a = 0; a < NDIM; a++) {

	    hier::Box<NDIM> face_box = pdat::FaceGeometry<NDIM>::toFaceBox(box,
									   a);
	    hier::Index<NDIM> flo = face_box.lower();
	    hier::Index<NDIM> fhi = face_box.upper();
	    for (int f = 0; f < 2; f++) {

	       pdat::FaceIndex<NDIM> findx(ifc(), a, f);

	       int offset = 0;
	       for (int i = NDIM-1; i > 0; i--) {
		  offset = (fhi(i-1)-flo(i-1)+1)*
		     (findx(i)-flo(i) + offset);
	       }
	       offset += findx(0)-flo(0);

	       if (!tbox::MathUtilities<double>::equalEps(face_data(findx), face_ptr[a][offset])) {
		  tbox::perr << "FAILED: - FaceIndex test" << endl;
		  ++error_count;
	       }
	    }
	 }
      }

      /*
       * Test NodeIterator
       */

      double* node_ptr = node_data.getPointer();

      hier::Box<NDIM> node_box = pdat::NodeGeometry<NDIM>::toNodeBox(box);

      hier::Index<NDIM> nlo = node_box.lower();
      hier::Index<NDIM> nhi = node_box.upper();

      int node_data_size = node_box.size();

      for (int i = 0; i < node_data_size; i++) {
	 node_ptr[i] = (double) i;
      }

      j = 0;
      for (pdat::NodeIterator<NDIM> ni(box); ni; ni++) {
	 if (!tbox::MathUtilities<double>::equalEps(node_data(ni()), node_ptr[j])) {
	    tbox::perr << "FAILED: - NodeIterator test" << endl;
	    ++error_count;
	 }
	 j++;
      }

      /*
       * Test NodeIndex
       */

      for (pdat::CellIterator<NDIM> inc(box); inc; inc++) {

	 hier::IntVector<NDIM> corner(0);

	 bool all_corners_complete = false;

	 while (!all_corners_complete) {

	    pdat::NodeIndex<NDIM> nindx(inc(), corner);

	    int offset = 0;
	    for (int i = NDIM-1; i > 0; i--) {
	       offset = (nhi(i-1)-nlo(i-1)+1)*
		  (nindx(i)-nlo(i) + offset);
	    }
	    offset += nindx(0)-nlo(0);

	    if (!tbox::MathUtilities<double>::equalEps(node_data(nindx), node_ptr[offset])) {
	       tbox::perr << "FAILED: - NodeIndex test" << endl;
	       ++error_count;
	    }

	    int u;
	    for (u = 0; u < NDIM; u++) {
	       if (corner(u) == 1) {
		  corner(u) = 0;
	       } else {
		  corner(u)++;
		  break;
	       }
	    }
	    if (u == NDIM) {
	       all_corners_complete = true;
	    }
	 }
      }

      /*
       * Test EdgeIterator
       */

      double* edge_ptr[NDIM];

      for (int axis = 0; axis < NDIM; axis++) {
                                                                                
	 edge_ptr[axis] = edge_data.getPointer(axis);
                                                                                
	 hier::Box<NDIM> edge_box = pdat::EdgeGeometry<NDIM>::toEdgeBox(box,
									axis);
	 int edge_data_size = edge_box.size();
                                                                                
	 for (int i = 0; i < edge_data_size; i++) {
	    edge_ptr[axis][i] = (double) i;
	 }

	 j = 0;
	 for (pdat::EdgeIterator<NDIM> ei(box, axis); ei; ei++) {
	    if (!tbox::MathUtilities<double>::equalEps(edge_data(ei()), edge_ptr[axis][j])) {
	       tbox::perr << "FAILED: - EdgeIterator test" << endl;
	       ++error_count;
	    }
	    j++;
	 }
      }

      /*
       * Test EdgeIndex
       */

      for (pdat::CellIterator<NDIM> iec(box); iec; iec++) {
	 for (int a = 0; a < NDIM; a++) {

	    hier::Box<NDIM> edge_box = pdat::EdgeGeometry<NDIM>::toEdgeBox(box,
									   a);
	    hier::Index<NDIM> elo = edge_box.lower();
	    hier::Index<NDIM> ehi = edge_box.upper();

	    for (int f = 0; f < (1 << (NDIM-1)); f++) {
	       pdat::EdgeIndex<NDIM> eindx(iec(), a, f);

	       int offset = 0;
	       for (int i = NDIM-1; i > 0; i--) {
		  offset = (ehi(i-1)-elo(i-1)+1)*
		     (eindx(i)-elo(i) + offset);
	       }
	       offset += eindx(0)-elo(0);

	       if (!tbox::MathUtilities<double>::equalEps(edge_data(eindx), edge_ptr[a][offset])) {
		  tbox::perr << "FAILED: - EdgeIndex test" << endl;
		  ++error_count;
	       }
	    }
	 }
      }

      /*
       * Test SideIterator
       */

      double* side_ptr[NDIM];

      for (int axis = 0; axis < NDIM; axis++) {
                                                                                
	 side_ptr[axis] = side_data.getPointer(axis);
                                                                                
	 hier::Box<NDIM> side_box = pdat::SideGeometry<NDIM>::toSideBox(box,
									axis);
	 int side_data_size = side_box.size();
                                                                                
	 for (int i = 0; i < side_data_size; i++) {
	    side_ptr[axis][i] = (double) i;
	 }

	 j = 0;
	 for (pdat::SideIterator<NDIM> si(box, axis); si; si++) {
	    if (!tbox::MathUtilities<double>::equalEps(side_data(si()), side_ptr[axis][j])) {
	       tbox::perr << "FAILED: - SideIterator test" << endl;
	       ++error_count;
	    }
	    j++;
	 }
      }

      /*
       * Test SideIndex
       */

      for (pdat::CellIterator<NDIM> isc(box); isc; isc++) {
	 for (int a = 0; a < NDIM; a++) {

	    hier::Box<NDIM> side_box = pdat::SideGeometry<NDIM>::toSideBox(box,
									   a);
	    hier::Index<NDIM> slo = side_box.lower();
	    hier::Index<NDIM> shi = side_box.upper();

	    for (int f = 0; f < 2; f++) {
	       pdat::SideIndex<NDIM> sindx(isc(), a, f);

	       int offset = 0;
	       for (int i = NDIM-1; i > 0; i--) {
		  offset = (shi(i-1)-slo(i-1)+1)*
		     (sindx(i)-slo(i) + offset);
	       }
	       offset += sindx(0)-slo(0);

	       if (!tbox::MathUtilities<double>::equalEps(side_data(sindx), side_ptr[a][offset])) {
		  tbox::perr << "FAILED: - SideIndex test" << endl;
		  ++error_count;
	       }
	    }
	 }
      }

      if ( error_count == 0 ) {
	 tbox::pout << "\nPASSED:  dataaccess" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();

   return( error_count );
}
