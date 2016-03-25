//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/trunk/source/geometry/skeleton/grid_geom/BlockGridGeometry.C $
// Package:	SAMRAI geometry package
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 878 $
// Modified:	$LastChangedDate: 2006-01-09 16:55:30 -0800 (Mon, 09 Jan 2006) $
// Description: Simple Block grid geometry for an AMR hierarchy.
//

#ifndef included_geom_BlockGridGeometry_C
#define included_geom_BlockGridGeometry_C

#include "BlockGridGeometry.h"
#include <stdlib.h>

#include <fstream>

#include "BlockPatchGeometry.h"
#include "BoundaryLookupTable.h"
#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

// Time interpolation operators
#ifdef HAVE_DCOMPLEX
#include "CellComplexLinearTimeInterpolateOp.h"
#include "SideComplexLinearTimeInterpolateOp.h"
#include "FaceComplexLinearTimeInterpolateOp.h"
#include "NodeComplexLinearTimeInterpolateOp.h"
#include "OuterfaceComplexLinearTimeInterpolateOp.h"
#include "OutersideComplexLinearTimeInterpolateOp.h"
#include "SideComplexLinearTimeInterpolateOp.h"
#endif

#ifdef HAVE_FLOAT
#include "CellFloatLinearTimeInterpolateOp.h"
#include "SideFloatLinearTimeInterpolateOp.h"
#include "FaceFloatLinearTimeInterpolateOp.h"
#include "NodeFloatLinearTimeInterpolateOp.h"
#include "OuterfaceFloatLinearTimeInterpolateOp.h"
#include "OutersideFloatLinearTimeInterpolateOp.h"
#include "SideFloatLinearTimeInterpolateOp.h"
#endif

#include "CellDoubleLinearTimeInterpolateOp.h"
#include "SideDoubleLinearTimeInterpolateOp.h"
#include "FaceDoubleLinearTimeInterpolateOp.h"
#include "NodeDoubleLinearTimeInterpolateOp.h"
#include "OuterfaceDoubleLinearTimeInterpolateOp.h"
#include "OutersideDoubleLinearTimeInterpolateOp.h"
#include "SideDoubleLinearTimeInterpolateOp.h"


#include "PatchLevel.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"


#define GEOM_SKELETON_GRID_GEOMETRY_VERSION (2)

namespace SAMRAI {
    namespace geom {

/*
*************************************************************************
*                                                                       *
* Constructors for BlockGridGeometry.  Both set up operator     *
* handlers.  However, one initializes data members based on arguments.  *
* The other initializes the object based on input file information.     *
*                                                                       *
*************************************************************************
*/
template<int DIM>  BlockGridGeometry<DIM>::BlockGridGeometry(
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   const int block_number,
   bool register_for_restart)
:  xfer::Geometry<DIM>(object_name)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( !object_name.empty() );
   TBOX_ASSERT( !input_db.isNull() );
#endif

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   d_block_number = block_number;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }

   makeStandardOperators();

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if ( is_from_restart ) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

   hier::BoundaryLookupTable<DIM>
      ::setUsingOriginalLocations(d_using_original_locations);

}

template<int DIM>  BlockGridGeometry<DIM>::BlockGridGeometry(
   const std::string& object_name, 
   const hier::BoxArray<DIM>& domain,
   bool register_for_restart) 
:  xfer::Geometry<DIM>(object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
#endif

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
         registerRestartItem(d_object_name, this);
   }

   setPhysicalDomain(domain);

   makeStandardOperators();

   d_using_original_locations = true;
   hier::BoundaryLookupTable<DIM>::
      setUsingOriginalLocations(d_using_original_locations);

}

/*
*************************************************************************
*                                                                       *
* Destructor for BlockGridGeometry deallocates grid storage.    *
* Note that operator handlers that are created in constructor are       *
* deallocated in xfer::Geometry<DIM> destructor.                             *
*                                                                       *
*************************************************************************
*/

template<int DIM>  BlockGridGeometry<DIM>::~BlockGridGeometry()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}

/*
*************************************************************************
*                                                                       *
* Create and return pointer to refined version of this Cartesian        *
* grid geometry object refined by the given ratio.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::GridGeometry<DIM> > 
BlockGridGeometry<DIM>::makeRefinedGridGeometry(
   const std::string& fine_geom_name,
   const hier::IntVector<DIM>& refine_ratio,
   bool register_for_restart) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!fine_geom_name.empty());
   TBOX_ASSERT(fine_geom_name != d_object_name);
   TBOX_ASSERT(refine_ratio > hier::IntVector<DIM>(0));
#endif

   hier::BoxArray<DIM> fine_domain(this -> getPhysicalDomain());
   fine_domain.refine(refine_ratio);

   geom::BlockGridGeometry<DIM>* fine_geometry =
      new geom::BlockGridGeometry<DIM>(fine_geom_name,
                                       fine_domain,
                                       register_for_restart);

   fine_geometry->initializePeriodicShift(this->getPeriodicShift());

   return(tbox::Pointer<hier::GridGeometry<DIM> >(fine_geometry));
}

/*
*************************************************************************
*                                                                       *
* Create and return pointer to coarsened version of this Cartesian      *
* grid geometry object coarsened by the given ratio.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> tbox::Pointer<hier::GridGeometry<DIM> > 
BlockGridGeometry<DIM>::makeCoarsenedGridGeometry(
   const std::string& coarse_geom_name,
   const hier::IntVector<DIM>& coarsen_ratio,
   bool register_for_restart) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!coarse_geom_name.empty());
   TBOX_ASSERT(coarse_geom_name != d_object_name);
   TBOX_ASSERT(coarsen_ratio > hier::IntVector<DIM>(0));
#endif

   hier::BoxArray<DIM> coarse_domain(this -> getPhysicalDomain());
   coarse_domain.coarsen(coarsen_ratio);

   /*
    * Need to check that domain can be coarsened by given ratio.
    */
   const hier::BoxArray<DIM>& fine_domain = this -> getPhysicalDomain();
   const int nboxes = fine_domain.getNumberOfBoxes();
   for (int ib = 0; ib < nboxes; ib++) {
      hier::Box<DIM> testbox = hier::Box<DIM>::refine(coarse_domain[ib], coarsen_ratio);
      if (testbox != fine_domain[ib]) {
#ifdef DEBUG_CHECK_ASSERTIONS
         tbox::plog << "BlockGridGeometry::makeCoarsenedGridGeometry : Box # " << ib << std::endl;
         tbox::plog << "      fine box = " << fine_domain[ib] << std::endl;
         tbox::plog << "f      coarse box = " << coarse_domain[ib] << std::endl;
         tbox::plog << "      refined coarse box = " << testbox << std::endl;
#endif
         TBOX_ERROR("geom::BlockGridGeometry::makeCoarsenedGridGeometry() error...\n"
                    << "    geometry object with name = " << d_object_name
                    << "\n    Cannot be coarsened by ratio " << coarsen_ratio << std::endl);
      }
   }

   geom::BlockGridGeometry<DIM>* coarse_geometry =
      new geom::BlockGridGeometry<DIM>(coarse_geom_name,
                                     coarse_domain,
                                     register_for_restart);

   coarse_geometry->initializePeriodicShift(this->getPeriodicShift());

   return(tbox::Pointer<hier::GridGeometry<DIM> >(coarse_geometry));
}

/*
*************************************************************************
*                                                                       *
* Create default interlevel transfer operator handlers and time         *
* interpolation operator handlers.  Add them to appropriate chains.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BlockGridGeometry<DIM>::makeStandardOperators()
{
   /*
    * Standard linear time interpolation operators.
    */
#ifdef HAVE_DCOMPLEX
   addTimeInterpolateOperator(new pdat::CellComplexLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::SideComplexLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::FaceComplexLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::NodeComplexLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::OuterfaceComplexLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::OutersideComplexLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::SideComplexLinearTimeInterpolateOp<DIM>());
#endif

#ifdef HAVE_FLOAT
   addTimeInterpolateOperator(new pdat::CellFloatLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::SideFloatLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::FaceFloatLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::NodeFloatLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::OuterfaceFloatLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::OutersideFloatLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::SideFloatLinearTimeInterpolateOp<DIM>());
#endif

   addTimeInterpolateOperator(new pdat::CellDoubleLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::SideDoubleLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::FaceDoubleLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::NodeDoubleLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::OuterfaceDoubleLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::OutersideDoubleLinearTimeInterpolateOp<DIM>());
   addTimeInterpolateOperator(new pdat::SideDoubleLinearTimeInterpolateOp<DIM>());


}

/*
*************************************************************************
*                                                                       *
* Create BlockPatchGeometry geometry object, initializing its   *
* boundary and assigning it to the given patch.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BlockGridGeometry<DIM>::setGeometryDataOnPatch(
   hier::Patch<DIM>& patch, 
   const hier::IntVector<DIM>& ratio_to_level_zero, 
   const tbox::Array< tbox::Array<bool> >& touches_regular_bdry,
   const tbox::Array< tbox::Array<bool> >& touches_periodic_bdry) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < DIM; i++) {
      TBOX_ASSERT( ratio_to_level_zero(i) != 0 );
   }
   if (DIM > 1) {
      for (i = 0; i < DIM; i++) {
	 TBOX_ASSERT( (ratio_to_level_zero(i)*ratio_to_level_zero((i+1)%DIM) > 0)
		 || (ratio_to_level_zero(i) == 1)
		 || (ratio_to_level_zero((i+1)%DIM) == 1) );
      }
   }
#endif

   tbox::Pointer<BlockPatchGeometry<DIM> >
      geometry = new BlockPatchGeometry<DIM>(ratio_to_level_zero,
                                             d_block_number,
                                             touches_regular_bdry,
                                             touches_periodic_bdry);

   patch.setPatchGeometry(geometry); 

}


/*
*************************************************************************
*                                                                       *
* Writes out version number and data members for the class.		*
*                                                                       *
*************************************************************************
*/

template<int DIM> void BlockGridGeometry<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif
   db->putInteger("GEOM_SKELETON_GRID_GEOMETRY_VERSION",
      GEOM_SKELETON_GRID_GEOMETRY_VERSION);
   tbox::Array<tbox::DatabaseBox> temp_box_array = this -> getPhysicalDomain();
   db->putDatabaseBoxArray("d_physical_domain", temp_box_array);

   hier::IntVector<DIM> level0_shift = getPeriodicShift(hier::IntVector<DIM>(1));
   int* temp_shift = level0_shift;
   db->putIntegerArray("d_periodic_shift", temp_shift, DIM);

   db->putBool("d_using_original_locations", d_using_original_locations);

}


/*
*************************************************************************
*                                                                       *
* Data is read from input only if the simulation is not from restart.   *
* Otherwise, all values specifed in the input database are ignored.	*
* In this method data from the database are read to local 		*
* variables and the setPhysicalDomain() method is called. 		*
*                                                                       *
*************************************************************************
*/

template<int DIM> void BlockGridGeometry<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db, 
   bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif
   
   if (!is_from_restart) {

      hier::BoxArray<DIM> domain;
      if (db->keyExists("domain_boxes")) {
         domain = db->getDatabaseBoxArray("domain_boxes");
         if (domain.getNumberOfBoxes() == 0) {
            TBOX_ERROR(d_object_name << ":  "
                       << "Block `domain_boxes' array found in input.");
         }
      } else {
         TBOX_ERROR(d_object_name << ":  "
                    << "Key data `domain_boxes' not found in input.");
      }

      int pbc[DIM];
      hier::IntVector<DIM> per_bc(0);
      if (db->keyExists("periodic_dimension")) {
         db->getIntegerArray("periodic_dimension", pbc, DIM);
         for (int i = 0; i < DIM; i++) {
            per_bc(i) = ((pbc[i] == 0) ? 0 : 1);
         }
      }

      if (DIM > 3) {
         d_using_original_locations = false;
      } else {
         d_using_original_locations = true;
      }
      if (db->keyExists("use_original_location_indices")) {
         d_using_original_locations =
            db->getBool("use_original_location_indices");
      }

      setPhysicalDomain(domain);

      initializePeriodicShift(per_bc);

     
   }
}

/*
*************************************************************************
*                                                                       *
* Checks to see if the version number for the class is the same as	*
* as the version number of the restart file.				*
* If they are equal, then the data from the database are read to local	*
* variables and the setPhysicalDomain() method is called. 		*
*                                                                       *
*************************************************************************
*/
template<int DIM> void BlockGridGeometry<DIM>::getFromRestart()
{
   tbox::Pointer<tbox::Database> restart_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;

   if ( restart_db ->isDatabase(d_object_name) ) {
      db = restart_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in the restart file.");
   }

   int ver = db->getInteger("GEOM_SKELETON_GRID_GEOMETRY_VERSION");
   if (ver != GEOM_SKELETON_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
              << "Restart file version is different than class version.");
   }
   hier::BoxArray<DIM> domain = db->getDatabaseBoxArray("d_physical_domain");

   setPhysicalDomain(domain);

   hier::IntVector<DIM> periodic_shift;
   int* temp_shift = periodic_shift;
   db->getIntegerArray("d_periodic_shift", temp_shift, DIM);
   initializePeriodicShift(periodic_shift);
 
   d_using_original_locations = db->getBool("d_using_original_locations");
}

/*
*************************************************************************
*                                                                       *
* Print BlockGridGeometry class data.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BlockGridGeometry<DIM>::printClassData(std::ostream& os) const
{
   os << "Printing BlockGridGeometry data: this = "
      << (BlockGridGeometry<DIM>*)this << std::endl;

   xfer::Geometry<DIM>::printClassData(os);
}

}
}
#endif
