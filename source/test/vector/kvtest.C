//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/vector/kvtest.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1800 $
// Modified:    $LastChangedDate: 2007-12-18 14:27:05 -0800 (Tue, 18 Dec 2007) $
// Description: Main program to test SAMRAI-KINSOL vector interface.
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include "tbox/SAMRAI_MPI.h"
#include <string>
using namespace std;
#include "tbox/PIO.h"

#include "tbox/SAMRAIManager.h"

#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellData.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "FaceVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "NodeData.h"
#include "NodeIndex.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "ProcessorMapping.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableContext.h"
#include "VariableDatabase.h"

#if HAVE_SUNDIALS

#include "Sundials_SAMRAIVector.h"
#include "SundialsAbstractVector.h"
#include "sundials/sundials_nvector.h"
#include "solv_NVector.h"

#endif




#ifndef NULL
#define NULL (0)
#endif

#define NCELL_VARS    2
#define NFACE_VARS    2
#define NNODE_VARS    4

using namespace SAMRAI;

int main( int argc, char *argv[] ) {

   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {


      tbox::PIO::logOnlyNodeZero("kvtest.log");

#if HAVE_SUNDIALS

      int ln, iv;

      // Make a dummy hierarchy domain
      double lo[3] = {0.0, 0.0, 0.0};
      double hi[3] = {1.0, 0.5, 0.5};

      hier::Box<3> coarse0(hier::Index<3>(0,0,0), hier::Index<3>(4,2,2));
      hier::Box<3> coarse1(hier::Index<3>(5,0,0), hier::Index<3>(9,2,2));
      hier::Box<3> coarse2(hier::Index<3>(0,0,3), hier::Index<3>(4,2,4));
      hier::Box<3> coarse3(hier::Index<3>(5,0,3), hier::Index<3>(9,2,4));
      hier::Box<3> coarse4(hier::Index<3>(0,3,0), hier::Index<3>(4,4,2));
      hier::Box<3> coarse5(hier::Index<3>(5,3,0), hier::Index<3>(9,4,2));
      hier::Box<3> coarse6(hier::Index<3>(0,3,3), hier::Index<3>(4,4,4));
      hier::Box<3> coarse7(hier::Index<3>(5,3,3), hier::Index<3>(9,4,4));
      hier::Box<3> fine0(hier::Index<3>(4,4,4), hier::Index<3>(7,5,5));
      hier::Box<3> fine1(hier::Index<3>(4,4,6), hier::Index<3>(7,5,7));
      hier::Box<3> fine2(hier::Index<3>(4,6,4), hier::Index<3>(7,7,5));
      hier::Box<3> fine3(hier::Index<3>(4,6,6), hier::Index<3>(7,7,7));
      hier::Box<3> fine4(hier::Index<3>(8,4,4), hier::Index<3>(13,5,5));
      hier::Box<3> fine5(hier::Index<3>(8,4,6), hier::Index<3>(13,5,7));
      hier::Box<3> fine6(hier::Index<3>(8,6,4), hier::Index<3>(13,7,5));
      hier::Box<3> fine7(hier::Index<3>(8,6,6), hier::Index<3>(13,7,7));
      hier::IntVector<3> ratio(2);

      hier::BoxArray<3> coarse_domain(8);
      hier::BoxArray<3> fine_boxes(8);
      coarse_domain(0) = coarse0;
      coarse_domain(1) = coarse1;
      coarse_domain(2) = coarse2;
      coarse_domain(3) = coarse3;
      coarse_domain(4) = coarse4;
      coarse_domain(5) = coarse5;
      coarse_domain(6) = coarse6;
      coarse_domain(7) = coarse7;
      fine_boxes(0) = fine0;
      fine_boxes(1) = fine1;
      fine_boxes(2) = fine2;
      fine_boxes(3) = fine3;
      fine_boxes(4) = fine4;
      fine_boxes(5) = fine5;
      fine_boxes(6) = fine6;
      fine_boxes(7) = fine7;

      hier::BoxList<3> coarse_domain_list(coarse_domain);
      hier::BoxList<3> fine_level_list(fine_boxes);
      coarse_domain_list.coalesceBoxes();
      fine_level_list.coalesceBoxes();

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(coarse_domain_list.getNumberItems() == 1);
      TBOX_ASSERT(fine_level_list.getNumberItems() == 1);
#endif

      hier::Box<3> coarse_domain_box = coarse_domain_list.getFirstItem();
      hier::Box<3> fine_level_box = fine_level_list.getFirstItem();

      tbox::Pointer<geom::CartesianGridGeometry<3> > geometry =
	 new geom::CartesianGridGeometry<3>("CartesianGeometry", lo, hi, coarse_domain);

      tbox::Pointer<hier::PatchHierarchy<3> > hierarchy = new hier::PatchHierarchy<3>("PatchHierarchy",
										      geometry);

      // Note: For these simple tests we allow at most 2 processors.
      const int nproc = tbox::SAMRAI_MPI::getNodes();
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(nproc < 3);
#endif
      const int n_coarse_boxes = coarse_domain.getNumberOfBoxes();
      const int n_fine_boxes = fine_boxes.getNumberOfBoxes();
      hier::ProcessorMapping mapping0(n_coarse_boxes);
      hier::ProcessorMapping mapping1(n_fine_boxes);

      int ib;
      for (ib = 0; ib < n_coarse_boxes; ib++) {
	 mapping0.setProcessorAssignment(ib, ib%nproc);
      }

      for (ib = 0; ib < n_fine_boxes; ib++) {
	 mapping1.setProcessorAssignment(ib, ib%nproc);
      }

      hierarchy->makeNewPatchLevel(0, hier::IntVector<3>(1), coarse_domain, mapping0);
      hierarchy->makeNewPatchLevel(1, ratio, fine_boxes, mapping1);

      // Create instance of hier::Variable<NDIM> database
      hier::VariableDatabase<3>* variable_db = hier::VariableDatabase<3>::getDatabase();
      tbox::Pointer< hier::VariableContext > dummy = variable_db->getContext("dummy");
      const hier::IntVector<3> no_ghosts(0);

      // Make some dummy variables and data on the hierarchy
      tbox::Pointer< pdat::CellVariable<3,double> > cvar[NCELL_VARS];
      int cvindx[NCELL_VARS];
      cvar[0] = new pdat::CellVariable<3,double>("cvar0", 2);
      cvindx[0] = variable_db->registerVariableAndContext(
	 cvar[0], dummy, no_ghosts);
      cvar[1] = new pdat::CellVariable<3,double>("cvar1", 2);
      cvindx[1] = variable_db->registerVariableAndContext(
	 cvar[1], dummy, no_ghosts);

      tbox::Pointer< pdat::CellVariable<3,double> >
	 cwgt = new pdat::CellVariable<3,double>("cwgt", 1);
      int cwgt_id = variable_db->registerVariableAndContext(
	 cwgt, dummy, no_ghosts);

      tbox::Pointer< pdat::FaceVariable<3,double> > fvar[NFACE_VARS];
      int fvindx[NFACE_VARS];
      fvar[0] = new pdat::FaceVariable<3,double>("fvar0", 1);
      fvindx[0] = variable_db->registerVariableAndContext(
	 fvar[0], dummy, no_ghosts);
      fvar[1] = new pdat::FaceVariable<3,double>("fvar1", 1);
      fvindx[1] = variable_db->registerVariableAndContext(
	 fvar[1], dummy, no_ghosts);

      tbox::Pointer< pdat::FaceVariable<3,double> >
	 fwgt = new pdat::FaceVariable<3,double>("fwgt", 1);
      int fwgt_id = variable_db->registerVariableAndContext(
	 fwgt, dummy, no_ghosts);

      tbox::Pointer< pdat::NodeVariable<3,double> > nvar[NNODE_VARS];
      int nvindx[NNODE_VARS];
      nvar[0] = new pdat::NodeVariable<3,double>("nvar0", 1);
      nvindx[0] = variable_db->registerVariableAndContext(
	 nvar[0], dummy, no_ghosts);
      nvar[1] = new pdat::NodeVariable<3,double>("nvar1", 1);
      nvindx[1] = variable_db->registerVariableAndContext(
	 nvar[1], dummy, no_ghosts);
      nvar[2] = new pdat::NodeVariable<3,double>("nvar2", 1);
      nvindx[2] = variable_db->registerVariableAndContext(
	 nvar[2], dummy, no_ghosts);
      nvar[3] = new pdat::NodeVariable<3,double>("nvar3", 1);
      nvindx[3] = variable_db->registerVariableAndContext(
	 nvar[3], dummy, no_ghosts);

      tbox::Pointer< pdat::NodeVariable<3,double> > 
	 nwgt = new pdat::NodeVariable<3,double>("nwgt", 1);
      int nwgt_id = variable_db->registerVariableAndContext(
	 nwgt, dummy, no_ghosts);

      for (ln = 0; ln < 2; ln++) {
	 hierarchy->getPatchLevel(ln)->allocatePatchData(cwgt_id);
	 hierarchy->getPatchLevel(ln)->allocatePatchData(fwgt_id);
	 hierarchy->getPatchLevel(ln)->allocatePatchData(nwgt_id);
      }

      tbox::Pointer< math::HierarchyCellDataOpsReal<3,double> > cell_ops = 
	 math::HierarchyDataOpsManager<3>::getManager()->getOperationsDouble(cwgt, 
									     hierarchy);
      tbox::Pointer< math::HierarchyFaceDataOpsReal<3,double> > face_ops = 
	 math::HierarchyDataOpsManager<3>::getManager()->getOperationsDouble(fwgt,
									     hierarchy);
      tbox::Pointer< math::HierarchyNodeDataOpsReal<3,double> > node_ops = 
	 math::HierarchyDataOpsManager<3>::getManager()->getOperationsDouble(nwgt,
									     hierarchy);

      cell_ops->resetLevels(0,1); 
      face_ops->resetLevels(0,1); 
      node_ops->resetLevels(0,1); 

      tbox::Pointer<hier::Patch<3> > patch;
      tbox::Pointer<geom::CartesianPatchGeometry<3> > pgeom;
      tbox::Pointer< pdat::CellData<3,double> > cdata; 
      tbox::Pointer< pdat::FaceData<3,double> > fdata; 
      tbox::Pointer< pdat::NodeData<3,double> > ndata; 

      // Set control volume data for vector components
      hier::Box<3> coarse_fine = fine_level_box;
      coarse_fine.coarsen(ratio);

      // Initialize control volume data for cell-centered components

      for (ln = 0; ln < 2; ln++) {
	 tbox::Pointer<hier::PatchLevel<3> > level = hierarchy->getPatchLevel(ln);
	 for (hier::PatchLevel<3>::Iterator ip(level); ip; ip++) {
	    patch = level->getPatch(ip());
	    pgeom = patch->getPatchGeometry();
	    const double* dx = pgeom->getDx();
	    const double cell_vol = dx[0]*dx[1]*dx[2];
	    tbox::Pointer< pdat::CellData<3,double> > cvdata = patch->getPatchData(cwgt_id);
	    cvdata->fillAll(cell_vol);
	    if (ln == 0) cvdata->fillAll(0.0, (coarse_fine * patch->getBox()) );
	 }
      }

      // Initialize control volume data for face-centered components
      for (ln = 0; ln < 2; ln++) {

	 tbox::Pointer<hier::PatchLevel<3> > level = hierarchy->getPatchLevel(ln);
	 for (hier::PatchLevel<3>::Iterator ip(level); ip; ip++) {
	    tbox::Pointer< pdat::FaceData<3,double> > data;
	    patch = level->getPatch(ip());
	    pgeom = patch->getPatchGeometry();
	    const double* dx = pgeom->getDx();
	    const double face_vol = dx[0]*dx[1]*dx[2];
	    data = patch->getPatchData(fwgt_id);
	    data->fillAll(face_vol);
	    pdat::FaceIndex<3> fi;
	    int plo0 = patch->getBox().lower(0);
	    int phi0 = patch->getBox().upper(0);
	    int plo1 = patch->getBox().lower(1);
	    int phi1 = patch->getBox().upper(1);
	    int plo2 = patch->getBox().lower(2);
	    int phi2 = patch->getBox().upper(2);
	    int ic, jc, kc;
	    double bdry_face_factor;
	    hier::Box<3> level_box;
 
	    if (ln == 0) {
	       data->fillAll(0.0, (coarse_fine * patch->getBox()) );
	       bdry_face_factor = 0.5;
	       level_box = coarse_domain_box;
	    } else {
	       bdry_face_factor = 1.5;
	       level_box = fine_level_box; 
	    }
	    //X face boundaries
	    if (plo0 == level_box.lower(0)) {
	       for (kc = plo2; kc <= phi2; kc++) {
		  for (jc = plo1; jc <= phi1; jc++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(plo0,jc,kc),
					     pdat::FaceIndex<3>::X,
					     pdat::FaceIndex<3>::Lower);
		     (*data)(fi) *= bdry_face_factor;
		  }
	       }
	    } else {
	       for (kc = plo2; kc <= phi2; kc++) {
		  for (jc = plo1; jc <= phi1; jc++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(plo0,jc,kc),
					     pdat::FaceIndex<3>::X,
					     pdat::FaceIndex<3>::Lower);
		     (*data)(fi) = 0.0;
		  }
	       }
	    }
	    if (phi0 == level_box.upper(0)) {
	       for (kc = plo2; kc <= phi2; kc++) {
		  for (jc = plo1; jc <= phi1; jc++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(phi0,jc,kc),
					     pdat::FaceIndex<3>::X,
					     pdat::FaceIndex<3>::Upper);
		     (*data)(fi) *= bdry_face_factor;
		  }
	       }
	    }

	    //Y face boundaries
	    if (plo1 == level_box.lower(1)) {
	       for (kc = plo2; kc <= phi2; kc++) {
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(ic,plo1,kc),
					     pdat::FaceIndex<3>::Y,
					     pdat::FaceIndex<3>::Lower);
		     (*data)(fi) *= bdry_face_factor;
		  }
	       }
	    } else {
	       for (kc = plo2; kc <= phi2; kc++) {
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(ic,plo1,kc),
					     pdat::FaceIndex<3>::Y,
					     pdat::FaceIndex<3>::Lower);
		     (*data)(fi) = 0.0;
		  }
	       }
	    }
	    if (phi1 == level_box.upper(1)) {
	       for (kc = plo2; kc <= phi2; kc++) {
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(ic,phi1,kc),
					     pdat::FaceIndex<3>::Y,
					     pdat::FaceIndex<3>::Upper);
		     (*data)(fi) *= bdry_face_factor;
		  }
	       }
	    }

	    //Z face boundaries
	    if (plo2 == level_box.lower(2)) {
	       for (jc = plo1; jc <= phi1; jc++) {
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(ic,jc,plo2),
					     pdat::FaceIndex<3>::Z,
					     pdat::FaceIndex<3>::Lower);
		     (*data)(fi) *= bdry_face_factor;
		  }
	       }
	    } else {
	       for (jc = plo1; jc <= phi1; jc++) {
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(ic,jc,plo2),
					     pdat::FaceIndex<3>::Z,
					     pdat::FaceIndex<3>::Lower);
		     (*data)(fi) = 0.0;
		  }
	       }
	    }
	    if (phi2 == level_box.upper(2)) {
	       for (jc = plo1; jc <= phi1; jc++) {
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::FaceIndex<3>(hier::Index<3>(ic,jc,phi2),
					     pdat::FaceIndex<3>::Z,
					     pdat::FaceIndex<3>::Upper);
		     (*data)(fi) *= bdry_face_factor;
		  }
	       }
	    }
	 }
      }

      for (ln = 0; ln < 2; ln++) {

	 tbox::Pointer<hier::PatchLevel<3> > level = hierarchy->getPatchLevel(ln);
	 for (hier::PatchLevel<3>::Iterator ip(level); ip; ip++) {
	    tbox::Pointer< pdat::NodeData<3,double> > data;
	    patch = level->getPatch(ip());
	    pgeom = patch->getPatchGeometry();
	    const double* dx = pgeom->getDx();
	    const double node_vol = dx[0]*dx[1]*dx[2];
	    data = patch->getPatchData(nwgt_id);
	    data->fillAll(node_vol);
	    pdat::NodeIndex<3> ni;
	    hier::Index<3> plo = patch->getBox().lower();
	    hier::Index<3> phi = patch->getBox().upper();
	    int ic,jc,kc;
	    double bdry_face_factor;
	    double bdry_edge_factor;
	    double bdry_node_factor;
	    hier::Box<3> level_box;
  
	    if (ln == 0) {
	       data->fillAll(0.0, (coarse_fine * patch->getBox()) );
	       bdry_face_factor = 0.5;
	       bdry_edge_factor = 0.25;
	       bdry_node_factor = 0.125;
	       level_box = coarse_domain_box;
	    } else {
	       bdry_face_factor = 1.5;
	       bdry_edge_factor = 2.25;
	       bdry_node_factor = 3.375;
	       level_box = fine_level_box;
	    }

	    //X faces
	    if (plo(0) == level_box.lower(0)) {
	       for (kc = plo(2); kc < phi(2); kc++) {
		  for (jc = plo(1); jc < phi(1); jc++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),jc,kc), pdat::NodeIndex<3>::LUU);
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    } else {
	       for (kc = plo(2); kc < phi(2); kc++) {
		  for (jc = plo(1); jc < phi(1); jc++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),jc,kc), pdat::NodeIndex<3>::LUU);
		     (*data)(ni) = 0.0;
		  }
	       }
	    }
            if (phi(0) == level_box.upper(0)) {
	       for (kc = plo(2); kc < phi(2); kc++) {
		  for (jc = plo(1); jc < phi(1); jc++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),jc,kc), pdat::NodeIndex<3>::UUU);
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }

	    //Y faces
	    if (plo(1) == level_box.lower(1)) {
	       for (kc = plo(2); kc < phi(2); kc++) {
		  for (ic = plo(0); ic < phi(0); ic++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(ic,plo(1),kc), pdat::NodeIndex<3>::ULU);
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    } else {
	       for (kc = plo(2); kc < phi(2); kc++) {
		  for (ic = plo(0); ic < phi(0); ic++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(ic,plo(1),kc), pdat::NodeIndex<3>::ULU);
		     (*data)(ni) = 0.0;
		  }
	       }
	    }
	    if (phi(1) == level_box.upper(1)) {
	       for (kc = plo(2); kc < phi(2); kc++) {
		  for (ic = plo(0); ic < phi(0); ic++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(ic,phi(1),kc), pdat::NodeIndex<3>::UUU);
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }

	    //Z faces
	    if (plo(2) == level_box.lower(2)) {
	       for (jc = plo(1); jc < phi(1); jc++) {
		  for (ic = plo(0); ic < phi(0); ic++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(ic,jc,plo(2)), pdat::NodeIndex<3>::UUL);
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    } else {
	       for (jc = plo(1); jc < phi(1); jc++) {
		  for (ic = plo(0); ic < phi(0); ic++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(ic,jc,plo(2)), pdat::NodeIndex<3>::UUL);
		     (*data)(ni) = 0.0;
		  }
	       }
	    }
	    if (phi(2) == level_box.upper(2)) {
	       for (jc = plo(1); jc < phi(1); jc++) {
		  for (ic = plo(0); ic < phi(0); ic++) {
		     ni = pdat::NodeIndex<3>(hier::Index<3>(ic,jc,phi(2)), pdat::NodeIndex<3>::UUU);
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }

	    // edge boundaries
	    for (ic = plo(0); ic < phi(0); ic++) {
	       ni = pdat::NodeIndex<3>(hier::Index<3>(ic,plo(1),plo(2)), pdat::NodeIndex<3>::ULL);
	       if (plo(1) == level_box.lower(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(ic,phi(1),plo(2)), pdat::NodeIndex<3>::UUL);
	       if (phi(1) == level_box.upper(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_face_factor;
		  } else { 
		     (*data)(ni) *= 0.0;
		  }
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(ic,plo(1),phi(2)), pdat::NodeIndex<3>::ULU);
	       if (plo(1) == level_box.lower(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(ic,phi(1),phi(2)), pdat::NodeIndex<3>::UUU);
	       if (phi(1) == level_box.upper(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }

	    for (jc = plo(1); jc < phi(1); jc++) {
	       ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),jc,plo(2)), pdat::NodeIndex<3>::LUL);
	       if (plo(0) == level_box.lower(0)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),jc,plo(2)), pdat::NodeIndex<3>::UUL);
	       if (phi(0) == level_box.upper(0)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_face_factor;
		  } else { 
		     (*data)(ni) *= 0.0;
		  }
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),jc,phi(2)), pdat::NodeIndex<3>::LUU);
	       if (plo(0) == level_box.lower(0)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),jc,phi(2)), pdat::NodeIndex<3>::UUU);
	       if (phi(0) == level_box.upper(0)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }

	    for (kc = plo(2); kc < phi(2); kc++) {
	       ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),plo(1),kc), pdat::NodeIndex<3>::LLU);
	       if (plo(0) == level_box.lower(0)) {
		  if (plo(1) == level_box.lower(1)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),plo(1),kc), pdat::NodeIndex<3>::ULU);
	       if (phi(0) == level_box.upper(0)) {
		  if (plo(1) == level_box.lower(1)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  if (plo(1) == level_box.lower(1)) {
		     (*data)(ni) *= bdry_face_factor;
		  } else { 
		     (*data)(ni) *= 0.0;
		  }
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),phi(1),kc), pdat::NodeIndex<3>::LUU);
	       if (plo(0) == level_box.lower(0)) {
		  if (phi(1) == level_box.upper(1)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	       ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),phi(1),kc), pdat::NodeIndex<3>::UUU);
	       if (phi(0) == level_box.upper(0)) {
		  if (phi(1) == level_box.upper(1)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  if (phi(1) == level_box.upper(1)) {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }

	    //corner boundaries
	    ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),plo(1),plo(2)), pdat::NodeIndex<3>::LLL);
	    if (plo(0) == level_box.lower(0)) {
	       if (plo(1) == level_box.lower(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	    } else {
	       (*data)(ni) = 0.0;
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),plo(1),phi(2)), pdat::NodeIndex<3>::LLU);
	    if (plo(0) == level_box.lower(0)) {
	       if (plo(1) == level_box.lower(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) *= bdry_edge_factor;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	    } else {
	       (*data)(ni) = 0.0;
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),phi(1),plo(2)), pdat::NodeIndex<3>::LUL);
	    if (plo(0) == level_box.lower(0)) {
	       if (phi(1) == level_box.upper(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       }
	    } else {
	       (*data)(ni) = 0.0;
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(plo(0),phi(1),phi(2)), pdat::NodeIndex<3>::LUU);
	    if (plo(0) == level_box.lower(0)) {
	       if (phi(1) == level_box.upper(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) *= bdry_edge_factor;
		  }
	       } else {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    } else {
	       (*data)(ni) = 0.0;
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),plo(1),plo(2)), pdat::NodeIndex<3>::ULL);
	    if (phi(0) == level_box.upper(0)) {
	       if (plo(1) == level_box.lower(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) *= 0.0;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	    } else {
	       if (plo(1) == level_box.lower(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),plo(1),phi(2)), pdat::NodeIndex<3>::ULU);
	    if (phi(0) == level_box.upper(0)) {
	       if (plo(1) == level_box.lower(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) *= bdry_edge_factor;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	    } else {
	       if (plo(1) == level_box.lower(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  (*data)(ni) = 0.0;
	       }
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),phi(1),plo(2)), pdat::NodeIndex<3>::UUL);
	    if (phi(0) == level_box.upper(0)) {
	       if (phi(1) == level_box.upper(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       }
	    } else {
	       if (phi(1) == level_box.upper(1)) {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       } else {
		  if (plo(2) == level_box.lower(2)) {
		     (*data)(ni) *= bdry_face_factor;
		  } else {
		     (*data)(ni) = 0.0;
		  }
	       }
	    }

	    ni = pdat::NodeIndex<3>(hier::Index<3>(phi(0),phi(1),phi(2)), pdat::NodeIndex<3>::UUU);
	    if (phi(0) == level_box.upper(0)) {
	       if (phi(1) == level_box.upper(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_node_factor;
		  } else {
		     (*data)(ni) *= bdry_edge_factor;
		  }
	       } else {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    } else {
	       if (phi(1) == level_box.upper(1)) {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_edge_factor;
		  } else {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       } else {
		  if (phi(2) == level_box.upper(2)) {
		     (*data)(ni) *= bdry_face_factor;
		  }
	       }
	    }
	 }
      }

      // Create SAMRAI vectors:
      // Each vector has four components (1 cell component with depth = 2,
      // 1 face component with depth = 1, and 2 node components with depth = 1). 
      tbox::Pointer< solv::SAMRAIVectorReal<3,double> > my_vec0 = 
	 new solv::SAMRAIVectorReal<3,double>("my_vec0", hierarchy, 0, 1);
      my_vec0->addComponent(cvar[0], cvindx[0], cwgt_id); 
      my_vec0->addComponent(fvar[0], fvindx[0], fwgt_id); 
      my_vec0->addComponent(nvar[0], nvindx[0], nwgt_id); 
      my_vec0->addComponent(nvar[1], nvindx[1], nwgt_id); 

      tbox::Pointer< solv::SAMRAIVectorReal<3,double> > my_vec1 = 
	 new solv::SAMRAIVectorReal<3,double>("my_vec1", hierarchy, 0, 1);
      my_vec1->addComponent(cvar[1], cvindx[1], cwgt_id);
      my_vec1->addComponent(fvar[1], fvindx[1], fwgt_id);
      my_vec1->addComponent(nvar[2], nvindx[2], nwgt_id);
      my_vec1->addComponent(nvar[3], nvindx[3], nwgt_id);

      my_vec0->allocateVectorData();
      my_vec1->allocateVectorData();


      // Print out control volume data and compute integrals...
 
      tbox::plog << "cell control volume data" << endl;
      cell_ops->printData(cwgt_id, tbox::plog);
      tbox::plog << "face control volume data" << endl;
      face_ops->printData(fwgt_id, tbox::plog);
      tbox::plog << "node control volume data" << endl;
      node_ops->printData(nwgt_id, tbox::plog);

      double norm;
      //pout << "sum of each control volume is " << endl;

      norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm,(double)0.5)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #1, norm != 0.5\n";
      }
      //pout << "Component 0 : " << norm << " = 0.5?" << endl;

      norm = face_ops->sumControlVolumes(fvindx[0], fwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm,(double)0.75)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #2, norm != 0.75\n";
      }
      //pout << "Component 1 : " << norm << " = 0.75?" << endl;

      norm = node_ops->sumControlVolumes(nvindx[0], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm,(double)0.25)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #3, norm != 0.25\n";
      }
      //pout << "Component 2 : " << norm << " = 0.25?" << endl;

      norm = node_ops->sumControlVolumes(nvindx[1], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm,(double)0.25)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #4, norm != 0.25\n";
      }
      //pout << "Component 3 : " << norm << " = 0.25?\n" << endl;

      // Simple tests of SAMRAI vector operations

      // Construct SAMRAI vector wrappers to test operations via Sundials calls
      N_Vector kvec0 = solv::Sundials_SAMRAIVector<3>::createSundialsVector(my_vec0) -> getNVector();
      N_Vector kvec1 = solv::Sundials_SAMRAIVector<3>::createSundialsVector(my_vec1) -> getNVector();

      double zero = 0.0; 
      double half = 0.5; 
      double one = 1.0; 
      double two = 2.0; 
      double three = 3.0; 
      double four = 4.0; 
      double twelve = 12.0; 

      // my_vec0 = 2.0
      my_vec0->setToScalar(2.0);
      N_VPrint_SAMRAI(kvec0);

      double my_norm;
      double p_norm;
      my_norm = my_vec0->L1Norm();
      //pout << "L1-norm of my_vec0 is " << norm << " = 6.0?\n" << endl;
      p_norm = N_VL1Norm(kvec0);
      //pout << "L1-norm of kvec0 is " << norm << " = 6.0?\n" << endl;
      if (!tbox::MathUtilities<double>::equalEps(my_norm, p_norm)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #5, L1-norm calculation\n";
      }

      // Set fine data in my_vec0 = 3.0
      cell_ops->resetLevels(1,1);
      face_ops->resetLevels(1,1);
      node_ops->resetLevels(1,1);
      cell_ops->setToScalar(cvindx[0], 3.0);
      face_ops->setToScalar(fvindx[0], 3.0);
      node_ops->setToScalar(nvindx[0], 3.0);
      node_ops->setToScalar(nvindx[1], 3.0);

      tbox::plog << "CHECK my_vec0" << endl;
      my_vec0->print(tbox::plog); 

      double my_min_val = my_vec0->min();
 
      double p_min_val;
      p_min_val = N_VMin(kvec0);
      tbox::plog << "min of kvec0 is " << p_min_val << " = 2.0?\n" << endl; 
      if (!tbox::MathUtilities<double>::equalEps(my_min_val, p_min_val)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #9, min val calculation\n";
      }

      // my_vec1 = 1/my_vec0
      my_vec1->reciprocal(my_vec0);

      double my_max_val = my_vec1->max();

      double p_max_val;
      p_max_val = N_VMaxNorm(kvec1);
      if (!tbox::MathUtilities<double>::equalEps(my_max_val, p_max_val)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #10, reciprocal max val calculation\n";
      }

      // Manipulate patch data on vector and test norms

      // Set some bogus values on Level in my_vec1 that should be masked out
      // in ensuing vector norm calculations

      tbox::Pointer<hier::PatchLevel<3> > level_zero 
	 = hierarchy->getPatchLevel(0);
      for (hier::PatchLevel<3>::Iterator ip(level_zero); ip; ip++) {
	 patch = level_zero->getPatch(ip());

	 cdata = patch->getPatchData(cvindx[1]);
	 hier::Index<3> cindex0(2,2,2);
	 hier::Index<3> cindex1(5,3,2);
	 hier::Index<3> cindex2(4,2,2);
	 hier::Index<3> cindex3(6,3,2);
	 if (patch->getBox().contains(cindex0)) {
	    (*cdata)(pdat::CellIndex<3>(cindex0), 0) = 100.0;
	 }
	 if (patch->getBox().contains(cindex1)) {
	    (*cdata)(pdat::CellIndex<3>(cindex1), 0) = -1000.0;
	 }
	 if (patch->getBox().contains(cindex2)) {
	    (*cdata)(pdat::CellIndex<3>(cindex2), 1) = 1100.0;
	 }
	 if (patch->getBox().contains(cindex3)) {
	    (*cdata)(pdat::CellIndex<3>(cindex3), 1) = -10.0;
	 }

	 fdata = patch->getPatchData(fvindx[1]);
	 hier::Index<3> findex0(2,2,2);
	 hier::Index<3> findex1(5,3,2);
	 if (patch->getBox().contains(findex0)) {
	    (*fdata)
	       (pdat::FaceIndex<3>(findex0, pdat::FaceIndex<3>::X, pdat::FaceIndex<3>::Lower)) = 200.0;
	 }
	 if (patch->getBox().contains(findex1)) {
	    (*fdata)
	       (pdat::FaceIndex<3>(findex1,pdat::FaceIndex<3>::Y,pdat::FaceIndex<3>::Upper)) = -2000.0;
	 }

	 hier::Index<3> nindex0(2,2,2);
	 hier::Index<3> nindex1(5,3,2);
	 if (patch->getBox().contains(nindex0)) {
	    ndata = patch->getPatchData(nvindx[2]);
	    (*ndata)(pdat::NodeIndex<3>(nindex0, pdat::NodeIndex<3>::LLL)) = 300.0;
	    ndata = patch->getPatchData(nvindx[3]);
	    (*ndata)(pdat::NodeIndex<3>(nindex0, pdat::NodeIndex<3>::LUL)) = 30.0;
	 }
	 if (patch->getBox().contains(nindex1)) {
	    ndata = patch->getPatchData(nvindx[2]);
	    (*ndata)(pdat::NodeIndex<3>(nindex1, pdat::NodeIndex<3>::UUL)) = -300.0;
	    ndata = patch->getPatchData(nvindx[3]);
	    (*ndata)(pdat::NodeIndex<3>(nindex1, pdat::NodeIndex<3>::ULL)) = -3300.0;
	 }
      }

      tbox::plog << "my_vec1 = 0.5 (& bogus values) on L0, = 0.3333 on L1?" << endl; 
      my_vec1->print(tbox::plog);

      double max_val = my_vec1->max();
      if (!tbox::MathUtilities<double>::equalEps(max_val, (double)1100.0)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #11, max bogus value\n";
      }

      double min_val = my_vec1->min();
      if (!tbox::MathUtilities<double>::equalEps(min_val, (double)-3300.0)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #12, min bogus value\n";
      }

      double my_dot = my_vec1->dot(my_vec1);

      double p_dot;
      p_dot = N_VDotProd(kvec1, kvec1);
      if (!tbox::MathUtilities<double>::equalEps(my_dot, p_dot)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #14, dot product calculation\n";
	 std::cout << "SGS " << my_dot << "," << p_dot <<std::endl;
      }

      my_norm = my_vec1->maxNorm();
      p_norm = N_VMaxNorm(kvec1);
      if (!tbox::MathUtilities<double>::equalEps(my_norm, p_norm)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #16, max norm calculation\n";
      }

      N_VConst(twelve, kvec0);
      norm = my_vec1->weightedL2Norm(my_vec0);
      if (!tbox::MathUtilities<double>::equalEps(norm, (double)7.6393717)) {
	 fail_count++;
      }
      norm = N_VWL2Norm(kvec1, kvec0);
      if (!tbox::MathUtilities<double>::equalEps(norm, (double)7.6393717)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #17, weighted L2 norm calculation\n";
      }

      norm = my_vec0->RMSNorm();
      if (!tbox::MathUtilities<double>::equalEps(norm, (double)12.0)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #18, RMS norm calculation\n";
      }

      norm = my_vec0->weightedRMSNorm(my_vec1);
      if (!tbox::MathUtilities<double>::equalEps(norm, (double)5.77482219887084)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #19, weighted RMS norm calculation\n";
      }
      norm = N_VWrmsNorm(kvec0, kvec1);
      if (!tbox::MathUtilities<double>::equalEps(norm, (double)5.77482219887084)) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #19, weighted RMS norm calculation\n";
      }

      // Vector test routines
      int test = my_vec0->testReciprocal(my_vec1);
      if (test != 1) {
	 fail_count++;
	 tbox::perr << "FAILED: - Test #20, reciprocal\n";
      }

      // NOTE: -1's on face and node components on level 1 are due to fact
      //       that these entries are redundant and are ignored in the
      //       test calculation.  Thus, they should be ignored in the
      //       operations that follow.

      // Duplicate vectors
      tbox::Pointer< solv::SAMRAIVectorReal<3,double> > my_vec2 = 
	 my_vec0->cloneVector("my_vec2");
      my_vec2->allocateVectorData();

      my_vec2->setRandomValues(1.0, 0.0); 

      N_Vector kvec2 = solv::Sundials_SAMRAIVector<3>::createSundialsVector(my_vec2) -> getNVector();

      tbox::plog << "\n\nPRINTING VARIABLE DATABASE before adding new vector" << endl;
      variable_db->printClassData(tbox::plog);

      N_Vector kvec3 = N_VClone(kvec2);

      tbox::Pointer< solv::SAMRAIVectorReal<3,double> > sam_vec3 =
	 solv::Sundials_SAMRAIVector<3>::getSAMRAIVector(kvec3);

      tbox::plog << "\nVariables and data components in new vector...";
      int ncomp = sam_vec3->getNumberOfComponents();
      for (int ic = 0; ic < ncomp; ic++) {
	 tbox::plog << "\n   Comp id, variable, data id = "
		    << ic << ", "
		    << sam_vec3->getComponentVariable(ic)->getName() << ", "
		    << sam_vec3->getComponentDescriptorIndex(ic);
      }

      tbox::plog << "\n\nPRINTING VARIABLE DATABASE after adding new vector" << endl;
      variable_db->printClassData(tbox::plog);

      N_VScale(one, kvec2, kvec3);

      tbox::plog << "kvec3 = Random....?" << endl;
//      N_VPrint(kvec3);

      N_VScale(one, kvec0, kvec3);
      N_VScale(one, kvec2, kvec0);

      N_VLinearSum(one, kvec0, four, kvec3, kvec0);

      N_VScale(half, kvec3, kvec3);

      N_VConst(one, kvec0);
      N_VConst(two, kvec1);
      N_VConst(three, kvec2);

      N_VLinearSum(three, kvec1, half, kvec0, kvec0);
      tbox::plog << "kvec0 = 3 * 2 + 0.5 * 1 = 6.5?" << endl;
//      N_VPrint(kvec0);

      N_VLinearSum(one, kvec0, twelve, kvec1, kvec1);
      tbox::plog << "kvec1 = 6.5 + 12 * 2 = 30.5?" << endl;
//      N_VPrint(kvec1);

      N_VLinearSum(zero, kvec0, one, kvec0, kvec0);
      tbox::plog << "kvec0 = 0 * 6.5 + 6.5 = 6.5?" << endl;
//      N_VPrint(kvec0);

      N_VProd(kvec2, kvec1, kvec0);
      tbox::plog << "kvec0 = 3.0 * 30.5 = 91.5?" << endl;
//      N_VPrint(kvec0);

      N_VDiv(kvec2, kvec1, kvec0);
      tbox::plog << "kvec0 = 3.0 / 30.5 = 0.098360656?" << endl;
//      N_VPrint(kvec0);

      my_vec0->setRandomValues(7.0, -5.0);
      tbox::plog << "kvec0 has random values between -5.0 and 2.0" << endl;
//      N_VPrint(kvec0);

      N_VAbs(kvec0, kvec0);
      tbox::plog << "all negative values from kvec0 have been made positive" << endl;
//      N_VPrint(kvec0);

      N_VInv(kvec1, kvec0);
      tbox::plog << "kvec0 = 1.0 / 30.5 = 0.032786885?" << endl;
//      N_VPrint(kvec0);

      N_VAddConst(kvec0, twelve, kvec0);
      tbox::plog << "kvec0 = .032786885 + 12.0 = 12.032786885?" << endl;
      N_VPrint_SAMRAI(kvec0);

      my_vec0->setRandomValues(5.0, 0.0);
      N_VCompare(2.5, kvec0, kvec1);
      tbox::plog << "entries of kvec1 are all either 1.0 or 0.0 " << endl <<
	 "at points with non-zero control volume" << endl;
      N_VPrint_SAMRAI(kvec1);

      test = N_VInvTest(kvec2, kvec0);
      tbox::plog << "is kvec0 pointwise inverted kvec2 (1.0/3.0)" << endl <<
	 "at points with non-zero control volume?" <<  endl;
      if (test == 0) tbox::plog << "kvec1 has at least one zero element" << endl;
      N_VPrint_SAMRAI(kvec0);

      my_vec0->setRandomValues(2.0, -1.0);
      N_VConst(two, kvec1);

      // No more tests....Destroy vectors and data...

      N_VDestroy(kvec3);

      tbox::plog << "\n\nPRINTING VARIABLE DATABASE after freeing new vector" << endl;
      variable_db->printClassData(tbox::plog);
   
      N_VDestroy(kvec0);
      N_VDestroy(kvec1);
      N_VDestroy(kvec2);

      // Deallocate vector data and control volumes
      my_vec0->freeVectorComponents();
      my_vec1->freeVectorComponents();
      my_vec2->freeVectorComponents();

      for (ln = 0; ln < 2; ln++) {
	 hierarchy->getPatchLevel(ln)->deallocatePatchData(cwgt_id);
	 hierarchy->getPatchLevel(ln)->deallocatePatchData(fwgt_id);
	 hierarchy->getPatchLevel(ln)->deallocatePatchData(nwgt_id);
      }

      for (iv = 0; iv < NCELL_VARS; iv++) {
	 cvar[iv].setNull();
      }
      for (iv = 0; iv < NFACE_VARS; iv++) {
	 fvar[iv].setNull();
      }
      for (iv = 0; iv < NNODE_VARS; iv++) {
	 nvar[iv].setNull();
      }
      cwgt.setNull();
      fwgt.setNull();
      nwgt.setNull();
      cell_ops.setNull();
      face_ops.setNull();
      node_ops.setNull();

      geometry.setNull();
      hierarchy.setNull();

#endif

      if ( fail_count == 0 ) {
	 tbox::pout << "\nPASSED:  kvtest" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(fail_count); 
}
