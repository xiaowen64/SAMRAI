/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple Skeleton grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_SkeletonGridGeometry_C
#define included_geom_SkeletonGridGeometry_C

#include "SAMRAI/geom/SkeletonGridGeometry.h"
#include "SAMRAI/geom/SkeletonPatchGeometry.h"
#include "SAMRAI/geom/SAMRAITransferOperatorRegistry.h"

// Skeleton no-operation refine and coarsen operators
#include "SAMRAI/geom/SkeletonCoarsen.h"
#include "SAMRAI/geom/SkeletonRefine.h"

// Time interpolation operators
#include "SAMRAI/pdat/CellComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideComplexLinearTimeInterpolateOp.h"

#include "SAMRAI/pdat/CellFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideFloatLinearTimeInterpolateOp.h"

#include "SAMRAI/pdat/CellDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/hier/BoxContainerConstIterator.h"
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cstdlib>
#include <fstream>

namespace SAMRAI {
namespace geom {

const int SkeletonGridGeometry::GEOM_SKELETON_GRID_GEOMETRY_VERSION = 2;

/*
 *************************************************************************
 *                                                                       *
 * Constructors for SkeletonGridGeometry.  Both set up operator     *
 * handlers.  However, one initializes data members based on arguments.  *
 * The other initializes the object based on input file information.     *
 *                                                                       *
 *************************************************************************
 */
SkeletonGridGeometry::SkeletonGridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   bool register_for_restart):
   hier::GridGeometry(dim, object_name,
                      tbox::Pointer<hier::TransferOperatorRegistry>(
                         new SAMRAITransferOperatorRegistry(dim)))
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!input_db.isNull());

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart && d_registered_for_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

}

SkeletonGridGeometry::SkeletonGridGeometry(
   const std::string& object_name,
   const hier::BoxList& domain,
   bool register_for_restart):
   hier::GridGeometry(domain.getDim(), object_name,
                      tbox::Pointer<hier::TransferOperatorRegistry>(
                         new SAMRAITransferOperatorRegistry(domain.getDim())))
{
   TBOX_ASSERT(!object_name.empty());

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   tbox::Array<hier::BoxList> domain_array(1, domain);
   this->setPhysicalDomain(domain_array);

}

/*
 *************************************************************************
 *                                                                       *
 * Destructor for SkeletonGridGeometry deallocates grid storage.         *
 *                                                                       *
 *************************************************************************
 */

SkeletonGridGeometry::~SkeletonGridGeometry()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(getObjectName());
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

tbox::Pointer<hier::GridGeometry>
SkeletonGridGeometry::makeRefinedGridGeometry(
   const std::string& fine_geom_name,
   const hier::IntVector& refine_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!fine_geom_name.empty());
   TBOX_ASSERT(fine_geom_name != getObjectName());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, refine_ratio);
   TBOX_ASSERT(refine_ratio > hier::IntVector::getZero(dim));

   hier::BoxList fine_domain(this->getPhysicalDomain(hier::BlockId(0)));
   fine_domain.refine(refine_ratio);

   geom::SkeletonGridGeometry* fine_geometry =
      new geom::SkeletonGridGeometry(fine_geom_name,
         fine_domain,
         register_for_restart);

   fine_geometry->initializePeriodicShift(this->getPeriodicShift(hier::
         IntVector::getOne(dim)));

   return tbox::Pointer<hier::GridGeometry>(fine_geometry);
}

/*
 *************************************************************************
 *                                                                       *
 * Create and return pointer to coarsened version of this Cartesian      *
 * grid geometry object coarsened by the given ratio.                    *
 *                                                                       *
 *************************************************************************
 */

tbox::Pointer<hier::GridGeometry>
SkeletonGridGeometry::makeCoarsenedGridGeometry(
   const std::string& coarse_geom_name,
   const hier::IntVector& coarsen_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!coarse_geom_name.empty());
   TBOX_ASSERT(coarse_geom_name != getObjectName());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, coarsen_ratio);
   TBOX_ASSERT(coarsen_ratio > hier::IntVector::getZero(dim));

   hier::BoxList coarse_domain(this->getPhysicalDomain(hier::BlockId(0)));
   coarse_domain.coarsen(coarsen_ratio);

   /*
    * Need to check that domain can be coarsened by given ratio.
    */
   const hier::BoxList& fine_domain = this->getPhysicalDomain(hier::BlockId(0));
   const int nboxes = fine_domain.size();
   hier::BoxList::ConstIterator fine_domain_itr(fine_domain);
   hier::BoxList::Iterator coarse_domain_itr(coarse_domain);
   for (int ib = 0; ib < nboxes; ib++, fine_domain_itr++, coarse_domain_itr++) {
      hier::Box testbox = hier::Box::refine(*coarse_domain_itr, coarsen_ratio);
      if (!testbox.isSpatiallyEqual(*fine_domain_itr)) {
         tbox::plog
         << "SkeletonGridGeometry::makeCoarsenedGridGeometry : Box # "
         << ib << std::endl;
         tbox::plog << "      fine box = " << *fine_domain_itr << std::endl;
         tbox::plog << "      coarse box = " << *coarse_domain_itr << std::endl;
         tbox::plog << "      refined coarse box = " << testbox << std::endl;

         TBOX_ERROR(
            "geom::SkeletonGridGeometry::makeCoarsenedGridGeometry() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Cannot be coarsened by ratio " << coarsen_ratio
            << std::endl);
      }
   }

   geom::SkeletonGridGeometry* coarse_geometry =
      new geom::SkeletonGridGeometry(coarse_geom_name,
         coarse_domain,
         register_for_restart);

   coarse_geometry->initializePeriodicShift(this->getPeriodicShift(hier::
         IntVector::getOne(dim)));

   return tbox::Pointer<hier::GridGeometry>(coarse_geometry);
}

/*
 *************************************************************************
 *                                                                       *
 * Create SkeletonPatchGeometry geometry object, initializing its   *
 * boundary and assigning it to the given patch.                         *
 *                                                                       *
 *************************************************************************
 */

void SkeletonGridGeometry::setGeometryDataOnPatch(
   hier::Patch& patch,
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const TwoDimBool& touches_periodic_bdry) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim, patch, ratio_to_level_zero);

   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */

   int i;
   for (i = 0; i < dim.getValue(); i++) {
      TBOX_ASSERT(ratio_to_level_zero(i) != 0);
   }
   if (dim > tbox::Dimension(1)) {
      for (i = 0; i < dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % dim.getValue()) == 1));
      }
   }
#endif

   tbox::Pointer<SkeletonPatchGeometry>
   geometry(new SkeletonPatchGeometry(ratio_to_level_zero,
               touches_regular_bdry,
               touches_periodic_bdry));

   patch.setPatchGeometry(geometry);
}

/*
 *************************************************************************
 *                                                                       *
 * Writes out version number and data members for the class.             *
 *                                                                       *
 *************************************************************************
 */

void SkeletonGridGeometry::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   const tbox::Dimension& dim(getDim());

   db->putInteger("GEOM_SKELETON_GRID_GEOMETRY_VERSION",
      GEOM_SKELETON_GRID_GEOMETRY_VERSION);
   tbox::Array<tbox::DatabaseBox> temp_box_array = this->getPhysicalDomain(hier::BlockId(0));
   db->putDatabaseBoxArray("d_physical_domain", temp_box_array);

   hier::IntVector level0_shift = this->getPeriodicShift(
         hier::IntVector::getOne(dim));
   int* temp_shift = &level0_shift[0];
   db->putIntegerArray("d_periodic_shift", temp_shift, dim.getValue());

}

/*
 *************************************************************************
 *                                                                       *
 * Data is read from input only if the simulation is not from restart.   *
 * Otherwise, all values specifed in the input database are ignored.     *
 * In this method data from the database are read to local               *
 * variables and the setPhysicalDomain() method is called.               *
 *                                                                       *
 *************************************************************************
 */

void SkeletonGridGeometry::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
   TBOX_ASSERT(!db.isNull());

   tbox::Dimension dim(getDim());

   if (!is_from_restart) {

      hier::BoxList domain(dim);
      if (db->keyExists("domain_boxes")) {
         domain = db->getDatabaseBoxArray("domain_boxes");
         if (domain.size() == 0) {
            TBOX_ERROR(
               getObjectName() << ":  "
                               << "Skeleton `domain_boxes' array found in input.");
         }
      } else {
         TBOX_ERROR(
            getObjectName() << ":  "
                            << "Key data `domain_boxes' not found in input.");
      }

      int pbc[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      hier::IntVector per_bc(dim, 0);
      if (db->keyExists("periodic_dimension")) {
         db->getIntegerArray("periodic_dimension", pbc, dim.getValue());
         for (int i = 0; i < dim.getValue(); i++) {
            per_bc(i) = ((pbc[i] == 0) ? 0 : 1);
         }
      }

      tbox::Array<hier::BoxList> domain_array(1, domain);

      this->setPhysicalDomain(domain_array);

      this->initializePeriodicShift(per_bc);
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Checks to see if the version number for the class is the same as      *
 * as the version number of the restart file.                            *
 * If they are equal, then the data from the database are read to local  *
 * variables and the setPhysicalDomain() method is called.               *
 *                                                                       *
 *************************************************************************
 */
void SkeletonGridGeometry::getFromRestart()
{
   const tbox::Dimension dim(getDim());

   tbox::Pointer<tbox::Database> restart_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;

   if (restart_db->isDatabase(getObjectName())) {
      db = restart_db->getDatabase(getObjectName());
   } else {
      TBOX_ERROR("Restart database corresponding to "
         << getObjectName() << " not found in the restart file.");
   }

   int ver = db->getInteger("GEOM_SKELETON_GRID_GEOMETRY_VERSION");
   if (ver != GEOM_SKELETON_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR(
         getObjectName() << ":  "
                         << "Restart file version is different than class version.");
   }
   hier::BoxList domain(db->getDatabaseBoxArray("d_physical_domain"));

   tbox::Array<hier::BoxList> domain_array(1, domain);

   this->setPhysicalDomain(domain_array);

   hier::IntVector periodic_shift(dim);
   int* temp_shift = &periodic_shift[0];
   db->getIntegerArray("d_periodic_shift", temp_shift, dim.getValue());
   this->initializePeriodicShift(periodic_shift);

}

/*
 *************************************************************************
 *                                                                       *
 * Print SkeletonGridGeometry class data.                           *
 *                                                                       *
 *************************************************************************
 */

void SkeletonGridGeometry::printClassData(
   std::ostream& os) const
{
   os << "Printing SkeletonGridGeometry data: this = "
      << (SkeletonGridGeometry *)this << std::endl;

   hier::GridGeometry::printClassData(os);
}

}
}
#endif
