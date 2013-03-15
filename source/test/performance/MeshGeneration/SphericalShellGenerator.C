/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   SphericalShellGenerator class implementation
 *
 ************************************************************************/
#include "SphericalShellGenerator.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iomanip>

#include <cmath>

using namespace SAMRAI;



SphericalShellGenerator::SphericalShellGenerator(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database> &database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(),
   d_radii(0)
{
   tbox::Array<double> init_disp;
   tbox::Array<double> velocity;
   tbox::Array<double> period;


   if (database) {

      if (database->isDouble("radii")) {

         d_radii = database->getDoubleVector("radii");

         if ( static_cast<int>(d_radii.size())%2 != 0 ) {
            d_radii.push_back(tbox::MathUtilities<double>::getMax());
         }

         tbox::plog << "SphericalShellGenerator radii:\n";
         for ( int i=0; i<static_cast<int>(d_radii.size()); ++i ) {
            tbox::plog << "\tradii[" << i << "] = " << d_radii[i] << '\n';
         }
      }

      if (database->isDouble("center")) {
         std::vector<double> tmpa = database->getDoubleVector("center");
         for ( int d=0; d<d_dim.getValue(); ++d ) {
            d_center[d] = tmpa[d];
         }
      }
      else {
         for (int d = 0; d < d_dim.getValue(); ++d) {
            d_center[d] = 0.0;
         }
      }

      /*
       * Input parameters to determine whether to tag by buffering
       * fronts or shrinking level, and by how much.
       */
      const std::string sname( "shrink_distance_" );
      const std::string bname( "buffer_distance_" );
      std::string aname;
      for ( int ln=0; ; ++ln ) {
         const std::string lnstr( tbox::Utilities::intToString(ln) );

         // Look for buffer input first, then shrink input.
         const std::string bnameln = bname + lnstr;
         const std::string snameln = sname + lnstr;

         std::vector<double> tmpa;

         if ( database->isDouble(bnameln) ) {
            tmpa = database->getDoubleVector(bnameln);
            if ( static_cast<int>(tmpa.size()) != dim.getValue() ) {
               TBOX_ERROR(bnameln << " input parameter must have " << dim << " values");
            }
            d_buffer_shrink.push_back('b');
         }

         if ( database->isDouble(snameln) ) {
            if ( !tmpa.empty() ) {
               TBOX_ERROR("Cannot specify both " << bnameln << " and " << snameln);
            }
            tmpa = database->getDoubleVector(snameln);
            if ( static_cast<int>(tmpa.size()) != dim.getValue() ) {
               TBOX_ERROR(snameln << " input parameter must have " << dim << " values");
            }
            d_buffer_shrink.push_back('s');
         }

         if ( !tmpa.empty() ) {
            d_buffer_shrink_distance.resize(d_buffer_shrink_distance.size() + 1);
            d_buffer_shrink_distance.back().insert( d_buffer_shrink_distance.back().end(),
                                                    &tmpa[0],
                                                    &tmpa[0]+static_cast<int>(tmpa.size()) );
         }
         else {
            break;
         }

      }

   }

   if ( d_buffer_shrink.empty() ) {
      TBOX_ERROR("SphericalShellGenerator: You must specify either\n"
                 << "buffer_distance_# or shrink_distance_# for each level\n"
                 << "you plan to use, except the finest level.");
   }

}



SphericalShellGenerator::~SphericalShellGenerator()
{
}




/*
 ***********************************************************************
 ***********************************************************************
 */
void SphericalShellGenerator::setTags(
   bool &exact_tagging,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int tag_ln,
   int tag_data_id )
{
   if ( d_buffer_shrink[tag_ln] == 's' ) {
      setTagsByShrinkingLevel(
         hierarchy,
         tag_ln,
         tag_data_id,
         hier::IntVector::getZero(d_dim),
         &d_buffer_shrink_distance[1][0]);
      exact_tagging = true;
      return;
   }

   const boost::shared_ptr<hier::PatchLevel> &tag_level(
      hierarchy->getPatchLevel(tag_ln));

   resetHierarchyConfiguration(hierarchy, 0, 1);

   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      boost::shared_ptr<hier::Patch> patch = *pi;

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         patch->getPatchGeometry(),
         BOOST_CAST_TAG);

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_data_id),
         BOOST_CAST_TAG);

      TBOX_ASSERT(patch_geom);
      TBOX_ASSERT(tag_data);

      tagShells( *tag_data, *patch_geom, d_buffer_shrink_distance[tag_ln] );

   }

   exact_tagging = false;

   return;
}



void SphericalShellGenerator::setDomain(
   hier::BoxContainer &domain,
   double xlo[],
   double xhi[],
   int autoscale_base_nprocs,
   const tbox::SAMRAI_MPI &mpi)
{
   TBOX_ASSERT( !domain.isEmpty() );
   NULL_USE(xlo);
   NULL_USE(xhi);

   if ( domain.size() != 1 ) {
      TBOX_ERROR("SphericalShellGenerator only supports single-box domains.");
   }

   hier::BoxContainer::const_iterator ii = domain.begin();
   hier::Box domain_box = *ii;
   hier::IntVector tmp_intvec = ii->numberCells();
   const tbox::Dimension &dim = domain.begin()->getDim();

   double scale_factor = static_cast<double>(mpi.getSize()) / autoscale_base_nprocs;
   double linear_scale_factor = pow( scale_factor, 1.0/dim.getValue() );

   for ( int d=0; d<dim.getValue(); ++d ) {
      // xhi[d] = xlo[d] + linear_scale_factor*(xhi[d]-xlo[d]);
      tmp_intvec(d) = static_cast<int>(0.5 + tmp_intvec(d)*linear_scale_factor);
   }
   tmp_intvec -= hier::IntVector::getOne(domain_box.getDim());
   tbox::plog << "SphericalShellGenerator::setDomain changing domain from "
              << domain_box << " to ";
   domain_box.upper() = domain_box.lower() + tmp_intvec;
   tbox::plog << domain_box << '\n';
   domain.pushBack(domain_box);

}



void SphericalShellGenerator::resetHierarchyConfiguration(
   /*! New hierarchy */ const boost::shared_ptr<hier::PatchHierarchy> &new_hierarchy,
   /*! Coarsest level */ const int coarsest_level,
   /*! Finest level */ const int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
   TBOX_ASSERT(new_hierarchy->getDim() == d_dim);
   d_hierarchy = new_hierarchy;
   TBOX_ASSERT(d_hierarchy);
}



/*
 * Compute the various data due to the fronts.
 */
void SphericalShellGenerator::tagShells(
   pdat::CellData<int> &tag_data,
   const geom::CartesianPatchGeometry &patch_geom,
   const std::vector<double> &buffer_distance ) const
{
   const tbox::Dimension &dim(tag_data.getDim());

   const hier::Box& pbox = tag_data.getBox();
   const hier::BlockId& block_id = pbox.getBlockId();
   const int tag_val = 1;

   const double *dx = patch_geom.getDx();
   const double *xlo = patch_geom.getXLower();

   // Compute the buffer in terms of cells.
   hier::IntVector buffer_cells(dim);
   for ( int i=0; i<dim.getValue(); ++i ) {
      buffer_cells(i) = static_cast<int>(0.5 + buffer_distance[i]/dx[i]);
   }

   /*
    * Compute tags on the nodes.  Tag node if it falls within one of
    * the shells.
    */
   pdat::NodeData<double> node_tag_data(pbox, 1, tag_data.getGhostCellWidth()+buffer_cells);
   node_tag_data.getArrayData().fillAll(0);

   pdat::NodeData<int>::iterator niend(pdat::NodeGeometry::end(node_tag_data.getGhostBox()));
   for (pdat::NodeData<int>::iterator ni(pdat::NodeGeometry::begin(node_tag_data.getGhostBox()));
        ni != niend; ++ni) {
      const pdat::NodeIndex& idx = *ni;
      double r[SAMRAI::MAX_DIM_VAL];
      double rr = 0;
      for (int d = 0; d < dim.getValue(); ++d) {
         r[d] = xlo[d] + dx[d] * ( idx(d) - pbox.lower()(d) ) - d_center[d];
         rr += r[d] * r[d];
      }
      rr = sqrt(rr);
      for (int i = 0; i < static_cast<int>(d_radii.size()); i += 2) {
         if (d_radii[i] <= rr && rr < d_radii[i + 1]) {
            node_tag_data(idx) = tag_val;
            break;
         }
      }
   }

   /*
    * Initialize tag_data to zero then tag specific cells.
    * Tag cell if any node in its buffered neighborhood is tagged.
    */
   tag_data.getArrayData().fillAll(0);

   pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(tag_data.getGhostBox()));
   for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(tag_data.getGhostBox()));
        ci != ciend; ++ci) {
      const pdat::CellIndex& cid = *ci;

      hier::Box check_box(cid,cid, block_id);
      check_box.grow(buffer_cells);
      pdat::NodeIterator node_itr_end(pdat::NodeGeometry::end(check_box));
      for ( pdat::NodeIterator node_itr(pdat::NodeGeometry::begin(check_box));
            node_itr != node_itr_end; ++node_itr ) {
         if ( node_tag_data(*node_itr) == tag_val ) {
            tag_data(cid) = tag_val;
            break;
         }
      }
   }

}



bool SphericalShellGenerator::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index) const
{
   (void)region;
   (void)depth_index;

   if (variable_name == "Tag value") {

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         patch.getPatchGeometry(),
         BOOST_CAST_TAG);
      TBOX_ASSERT(patch_geom);

      pdat::CellData<int> tag_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      tagShells(tag_data, *patch_geom, d_buffer_shrink_distance[patch.getPatchLevelNumber()]);
      pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(patch.getBox()));
      for (pdat::CellData<double>::iterator ci(pdat::CellGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = tag_data(*ci);
      }

   } else {
      TBOX_ERROR("Unrecognized name " << variable_name);
   }
   return true;
}
