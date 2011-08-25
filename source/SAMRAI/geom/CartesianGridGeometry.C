/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_CartesianGridGeometry_C
#define included_geom_CartesianGridGeometry_C

#include "SAMRAI/geom/CartesianGridGeometry.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/SAMRAITransferOperatorRegistry.h"

// Cell data coarsen operators
#include "SAMRAI/geom/CartesianCellComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellFloatWeightedAverage.h"

// Cell data refine operators
#include "SAMRAI/geom/CartesianCellComplexConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatLinearRefine.h"
#include "SAMRAI/pdat/CellComplexConstantRefine.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellFloatConstantRefine.h"
#include "SAMRAI/pdat/CellIntegerConstantRefine.h"

// Edge data coarsen operators
#include "SAMRAI/geom/CartesianEdgeComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeFloatWeightedAverage.h"

// Edge data refine operators
#include "SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianEdgeFloatConservativeLinearRefine.h"
#include "SAMRAI/pdat/EdgeComplexConstantRefine.h"
#include "SAMRAI/pdat/EdgeDoubleConstantRefine.h"
#include "SAMRAI/pdat/EdgeFloatConstantRefine.h"
#include "SAMRAI/pdat/EdgeIntegerConstantRefine.h"

// Face data coarsen operators
#include "SAMRAI/geom/CartesianFaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceFloatWeightedAverage.h"

// Face data refine operators
#include "SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianFaceFloatConservativeLinearRefine.h"
#include "SAMRAI/pdat/FaceComplexConstantRefine.h"
#include "SAMRAI/pdat/FaceDoubleConstantRefine.h"
#include "SAMRAI/pdat/FaceFloatConstantRefine.h"
#include "SAMRAI/pdat/FaceIntegerConstantRefine.h"

// Node data coarsen operators
#include "SAMRAI/pdat/NodeComplexInjection.h"
#include "SAMRAI/pdat/NodeDoubleInjection.h"
#include "SAMRAI/pdat/NodeFloatInjection.h"
#include "SAMRAI/pdat/NodeIntegerInjection.h"

// Node data refine operators
#include "SAMRAI/geom/CartesianNodeComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeFloatLinearRefine.h"

// Outerface data coarsen operators
#include "SAMRAI/geom/CartesianOuterfaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceFloatWeightedAverage.h"

// Outerface data refine operators
#include "SAMRAI/pdat/OuterfaceComplexConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceDoubleConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceFloatConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceIntegerConstantRefine.h"

// Outernode data coarsen operators
#include "SAMRAI/pdat/OuternodeDoubleConstantCoarsen.h"

// Outerside data coarsen operators
#include "SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h"

// Side data coarsen operators
#include "SAMRAI/geom/CartesianSideComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideFloatWeightedAverage.h"

// Side data refine operators
#include "SAMRAI/geom/CartesianSideDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianSideFloatConservativeLinearRefine.h"
#include "SAMRAI/pdat/SideComplexConstantRefine.h"
#include "SAMRAI/pdat/SideDoubleConstantRefine.h"
#include "SAMRAI/pdat/SideFloatConstantRefine.h"
#include "SAMRAI/pdat/SideIntegerConstantRefine.h"

// Time interpolation operators
#include "SAMRAI/pdat/CellComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/CellDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/CellFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideFloatLinearTimeInterpolateOp.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cstdlib>
#include <fstream>

#ifndef SAMRAI_INLINE
#include "SAMRAI/geom/CartesianGridGeometry.I"
#endif

namespace SAMRAI {
namespace geom {

const int CartesianGridGeometry::GEOM_CARTESIAN_GRID_GEOMETRY_VERSION = 2;

// using namespace std;

/*
 *************************************************************************
 *                                                                       *
 * Constructors for CartesianGridGeometry.  Both set up operator    *
 * handlers and register the geometry object with the RestartManager.    *
 * However, one initializes data members based on arguments.             *
 * The other initializes the object based on input file information.     *
 *                                                                       *
 *************************************************************************
 */
CartesianGridGeometry::CartesianGridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   tbox::Pointer<tbox::Database> input_db,
   bool register_for_restart):
   hier::GridGeometry(dim, object_name,
                      tbox::Pointer<hier::TransferOperatorRegistry>(
                         new SAMRAITransferOperatorRegistry(dim))),
   d_domain_box(dim)
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

CartesianGridGeometry::CartesianGridGeometry(
   const std::string& object_name,
   const double* x_lo,
   const double* x_up,
   const hier::BoxList& domain,
   bool register_for_restart):
   hier::GridGeometry(domain.getDim(), object_name,
                      tbox::Pointer<hier::TransferOperatorRegistry>(
                         new SAMRAITransferOperatorRegistry(domain.getDim()))),
   d_domain_box(domain.getDim())
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!(x_lo == (double *)NULL));
   TBOX_ASSERT(!(x_up == (double *)NULL));

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   setGeometryData(x_lo, x_up, domain);

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
CartesianGridGeometry::makeRefinedGridGeometry(
   const std::string& fine_geom_name,
   const hier::IntVector& refine_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension dim(getDim());

   TBOX_ASSERT(!fine_geom_name.empty());
   TBOX_ASSERT(fine_geom_name != getObjectName());
   TBOX_ASSERT(refine_ratio > hier::IntVector::getZero(dim));

   hier::BoxList fine_domain(this->getPhysicalDomain(hier::BlockId(0)));
   fine_domain.refine(refine_ratio);

   CartesianGridGeometry* fine_geometry =
      new CartesianGridGeometry(fine_geom_name,
         d_x_lo,
         d_x_up,
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

tbox::Pointer<hier::GridGeometry> CartesianGridGeometry::
makeCoarsenedGridGeometry(
   const std::string& coarse_geom_name,
   const hier::IntVector& coarsen_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!coarse_geom_name.empty());
   TBOX_ASSERT(coarse_geom_name != getObjectName());
   TBOX_ASSERT(coarsen_ratio > hier::IntVector::getZero(dim));

   hier::BoxList coarse_domain(this->getPhysicalDomain(hier::BlockId(0)));
   coarse_domain.coarsen(coarsen_ratio);

   /*
    * Need to check that domain can be coarsened by given ratio.
    */
   const hier::BoxList& fine_domain = this->getPhysicalDomain(hier::BlockId(0));
   const int nboxes = fine_domain.getNumberOfBoxes();
   hier::BoxList::Iterator fine_domain_itr(fine_domain);
   hier::BoxList::Iterator coarse_domain_itr(coarse_domain);
   for (int ib = 0; ib < nboxes; ib++, fine_domain_itr++, coarse_domain_itr++) {
      hier::Box testbox = hier::Box::refine(*coarse_domain_itr, coarsen_ratio);
      if (!testbox.isSpatiallyEqual(*fine_domain_itr)) {
#ifdef DEBUG_CHECK_ASSERTIONS
         tbox::plog
         << "CartesianGridGeometry::makeCoarsenedGridGeometry : Box # "
         << ib << std::endl;
         tbox::plog << "      fine box = " << *fine_domain_itr << std::endl;
         tbox::plog << "      coarse box = " << *coarse_domain_itr << std::endl;
         tbox::plog << "      refined coarse box = " << testbox << std::endl;
#endif
         TBOX_ERROR(
            "geom::CartesianGridGeometry::makeCoarsenedGridGeometry() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Cannot be coarsened by ratio " << coarsen_ratio
            << std::endl);
      }
   }

   hier::GridGeometry* coarse_geometry =
      new geom::CartesianGridGeometry(coarse_geom_name,
         d_x_lo,
         d_x_up,
         coarse_domain,
         register_for_restart);

   coarse_geometry->initializePeriodicShift(this->getPeriodicShift(hier::
         IntVector::getOne(dim)));

   return tbox::Pointer<hier::GridGeometry>(coarse_geometry);
}

/*
 *************************************************************************
 *                                                                       *
 * Destructor for CartesianGridGeometry deallocates grid storage.        *
 *                                                                       *
 *************************************************************************
 */

CartesianGridGeometry::~CartesianGridGeometry()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(getObjectName());
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Set data members for this geometry object based on arguments.         *
 *                                                                       *
 *************************************************************************
 */

void CartesianGridGeometry::setGeometryData(
   const double* x_lo,
   const double* x_up,
   const hier::BoxList& domain)
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!(x_lo == (double *)NULL));
   TBOX_ASSERT(!(x_up == (double *)NULL));
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, domain);

   for (int id = 0; id < dim.getValue(); id++) {
      d_x_lo[id] = x_lo[id];
      d_x_up[id] = x_up[id];
   }

   tbox::Array<hier::BoxList> domain_array(1, domain);
   this->setPhysicalDomain(domain_array);

   hier::Box bigbox(dim);
   for (hier::BoxList::Iterator k(getPhysicalDomain(hier::BlockId(0))); k; k++)
      bigbox += *k;

   d_domain_box = bigbox;

   hier::IntVector ncells = d_domain_box.numberCells();
   for (int id2 = 0; id2 < dim.getValue(); id2++) {
      double length = d_x_up[id2] - d_x_lo[id2];
      d_dx[id2] = length / ((double)ncells(id2));
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Create CartesianPatchGeometry geometry object, initializing its  *
 * boundary and grid information and assign it to the given patch.       *
 *                                                                       *
 *************************************************************************
 */

void CartesianGridGeometry::setGeometryDataOnPatch(
   hier::Patch& patch,
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const TwoDimBool& touches_periodic_bdry) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim, patch, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != hier::IntVector::getZero(dim));

   if (dim > tbox::Dimension(1)) {
      for (int i = 0; i < dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % dim.getValue()) == 1));
      }
   }
#endif

   double dx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   bool coarsen = false;
   if (ratio_to_level_zero(0) < 0) coarsen = true;
   hier::IntVector tmp_rat = (ratio_to_level_zero);
   for (int id2 = 0; id2 < dim.getValue(); id2++) {
      tmp_rat(id2) = abs(ratio_to_level_zero(id2));
   }

   hier::Box index_box = d_domain_box;
   hier::Box box = patch.getBox();

   if (coarsen) {
      index_box.coarsen(tmp_rat);
      for (int id3 = 0; id3 < dim.getValue(); id3++) {
         dx[id3] = d_dx[id3] * ((double)tmp_rat(id3));
      }
   } else {
      index_box.refine(tmp_rat);
      for (int id4 = 0; id4 < dim.getValue(); id4++) {
         dx[id4] = d_dx[id4] / ((double)tmp_rat(id4));
      }
   }

   for (int id5 = 0; id5 < dim.getValue(); id5++) {
      x_lo[id5] = d_x_lo[id5]
         + ((double)(box.lower(id5) - index_box.lower(id5))) * dx[id5];
      x_up[id5] = x_lo[id5] + ((double)box.numberCells(id5)) * dx[id5];
   }

   tbox::Pointer<CartesianPatchGeometry> geom(
      new CartesianPatchGeometry(ratio_to_level_zero,
         touches_regular_bdry,
         touches_periodic_bdry,
         dx, x_lo, x_up));

   patch.setPatchGeometry(geom);

}

/*
 *************************************************************************
 *                                                                       *
 * Print CartesianGridGeometry class data.                          *
 *                                                                       *
 *************************************************************************
 */

void CartesianGridGeometry::printClassData(
   std::ostream& os) const
{
   const tbox::Dimension& dim(getDim());

   os << "Printing CartesianGridGeometry data: this = "
      << (CartesianGridGeometry *)this << std::endl;
   os << "d_object_name = " << getObjectName() << std::endl;

   int id;
   os << "d_x_lo = ";
   for (id = 0; id < dim.getValue(); id++) {
      os << d_x_lo[id] << "   ";
   }
   os << std::endl;
   os << "d_x_up = ";
   for (id = 0; id < dim.getValue(); id++) {
      os << d_x_up[id] << "   ";
   }
   os << std::endl;
   os << "d_dx = ";
   for (id = 0; id < dim.getValue(); id++) {
      os << d_dx[id] << "   ";
   }
   os << std::endl;

   os << "d_domain_box = " << d_domain_box << std::endl;

   hier::GridGeometry::printClassData(os);
}

/*
 *************************************************************************
 *                                                                       *
 * Write class version number and object state to database.              *
 *                                                                       *
 *************************************************************************
 */

void CartesianGridGeometry::putToDatabase(
   tbox::Pointer<tbox::Database> db)
{
   TBOX_ASSERT(!db.isNull());

   const tbox::Dimension& dim(getDim());

   db->putInteger("GEOM_CARTESIAN_GRID_GEOMETRY_VERSION",
      GEOM_CARTESIAN_GRID_GEOMETRY_VERSION);
   tbox::Array<tbox::DatabaseBox> temp_box_array = this->getPhysicalDomain(hier::BlockId(0));
   db->putDatabaseBoxArray("d_physical_domain", temp_box_array);

   db->putDoubleArray("d_dx", d_dx, dim.getValue());
   db->putDoubleArray("d_x_lo", d_x_lo, dim.getValue());
   db->putDoubleArray("d_x_up", d_x_up, dim.getValue());

   hier::IntVector level0_shift(this->getPeriodicShift(
                                   hier::IntVector::getOne(dim)));
   int* temp_shift = &level0_shift[0];
   db->putIntegerArray("d_periodic_shift", temp_shift, dim.getValue());

}

/*
 *************************************************************************
 *                                                                       *
 * Data is read from input only if the simulation is not from restart.   *
 * Otherwise, all values specifed in the input database are ignored.     *
 * In this method data from the database are read to local               *
 * variables and the setGeometryData() method is called to               *
 * initialize the data members.                                          *
 *                                                                       *
 *************************************************************************
 */

void CartesianGridGeometry::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart)
{
   TBOX_ASSERT(!db.isNull());

   const tbox::Dimension& dim(getDim());

   if (!is_from_restart) {

      hier::BoxList domain(dim);
      if (db->keyExists("domain_boxes")) {
         domain = db->getDatabaseBoxArray("domain_boxes");
         if (domain.getNumberOfBoxes() == 0) {
            TBOX_ERROR(
               "CartesianGridGeometry::getFromInput() error...\n"
               << "    geometry object with name = " << getObjectName()
               << "\n    Empty `domain_boxes' array found in input."
               << std::endl);
         }
      } else {
         TBOX_ERROR("CartesianGridGeometry::getFromInput() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Key data `domain_boxes' not found in input." << std::endl);
      }

      double x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE],
             x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

      if (db->keyExists("x_lo")) {
         db->getDoubleArray("x_lo", x_lo, dim.getValue());
      } else {
         TBOX_ERROR("CartesianGridGeometry::getFromInput() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Key data `x_lo' not found in input." << std::endl);
      }
      if (db->keyExists("x_up")) {
         db->getDoubleArray("x_up", x_up, dim.getValue());
      } else {
         TBOX_ERROR("CartesianGridGeometry::getFromInput() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n   Key data `x_up' not found in input." << std::endl);
      }

      int pbc[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      hier::IntVector per_bc(dim, 0);
      if (db->keyExists("periodic_dimension")) {
         db->getIntegerArray("periodic_dimension", pbc, dim.getValue());
         for (int i = 0; i < dim.getValue(); i++) {
            per_bc(i) = ((pbc[i] == 0) ? 0 : 1);
         }
      }

      setGeometryData(x_lo, x_up, domain);

      this->initializePeriodicShift(per_bc);

   }
}

/*
 *************************************************************************
 *                                                                       *
 * Checks to see if the version number for the class is the same as      *
 * as the version number of the restart file.                            *
 * If they are equal, then the data from the database are read to local  *
 * variables and the setGeometryData() method is called to               *
 * initialize the data members.                                          *
 *                                                                       *
 *************************************************************************
 */
void CartesianGridGeometry::getFromRestart()
{
   const tbox::Dimension& dim(getDim());

   tbox::Pointer<tbox::Database> restart_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;

   if (restart_db->isDatabase(getObjectName())) {
      db = restart_db->getDatabase(getObjectName());
   } else {
      TBOX_ERROR("CartesianGridGeometry::getFromRestart() error...\n"
         << "    database with name " << getObjectName()
         << " not found in the restart file" << std::endl);
   }

   int ver = db->getInteger("GEOM_CARTESIAN_GRID_GEOMETRY_VERSION");
   if (ver != GEOM_CARTESIAN_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR("CartesianGridGeometry::getFromRestart() error...\n"
         << "    geometry object with name = " << getObjectName()
         << "Restart file version is different than class version" << std::endl);
   }
   hier::BoxList domain(db->getDatabaseBoxArray("d_physical_domain"));
   double x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE],
          x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   db->getDoubleArray("d_x_lo", x_lo, dim.getValue());
   db->getDoubleArray("d_x_up", x_up, dim.getValue());

   setGeometryData(x_lo, x_up, domain);

   hier::IntVector periodic_shift(dim);
   int* temp_shift = &periodic_shift[0];
   db->getIntegerArray("d_periodic_shift", temp_shift, dim.getValue());
   this->initializePeriodicShift(periodic_shift);

}

}
}
#endif
