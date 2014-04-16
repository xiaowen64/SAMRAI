/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   SphericalShellGenerator class declaration
 *
 ************************************************************************/
#ifndef included_SphericalShellGenerator
#define included_SphericalShellGenerator

#include "MeshGenerationStrategy.h"

#include <string>

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"

#include "boost/shared_ptr.hpp"

using namespace SAMRAI;

/*!
 * @brief Class to tag spherical shell patternse in given domain.
 *
 * Inputs:
 *
 * radii:
 * Starting at shell origin, tag cells that intersect regions defined by
 * radii[0]<r<radii[1], radii[2]<r<radii[3], radii[2i]<r<radii[2i+1] and so on.
 *
 * buffer_distance_0, buffer_distance_1, ...:
 * buffer_distance[ln] is the buffer distance when tagging ON
 * level ln.  We tag the shells and buffer the tags by this amount.
 */
class SphericalShellGenerator:
   public MeshGenerationStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   SphericalShellGenerator(
      /*! Ojbect name */
      const std::string& object_name,
      const tbox::Dimension& dim,
      /*! Input database */
      const boost::shared_ptr<tbox::Database> &database = boost::shared_ptr<tbox::Database>() );

   ~SphericalShellGenerator();

   /*!
    * @brief Set tas on the tag level.
    */
   virtual void setTags(
      bool &exact_tagging,
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int tag_ln,
      int tag_data_id);

   //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy virtuals

public:

   /*!
    * @brief Set the domain, possibly scaling up the specifications.
    *
    * Take the domain_boxes, xlo and xhi to be the size for the
    * (integer) value of autoscale_base_nprocs.  Scale the problem
    * from there to the number of process running by doubling the
    * size starting with the j direction.
    *
    * The number of processes must be a power of 2 times the value
    * of autoscale_base_nprocs.
    */
   void setDomain(
      hier::BoxContainer &domain,
      double xlo[],
      double xhi[],
      int autoscale_base_nprocs,
      const tbox::SAMRAI_MPI &mpi);

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const boost::shared_ptr<hier::PatchHierarchy> &new_hierarchy,
      /*! Coarsest level */ const int coarsest_level,
      /*! Finest level */ const int finest_level);

   //@}

   void
   initializePatchData(
      hier::Patch& patch,
      const double init_data_time,
      const bool initial_time,
      const bool allocate_data)
      {
         NULL_USE(patch);
         NULL_USE(init_data_time);
         NULL_USE(initial_time);
         NULL_USE(allocate_data);
         TBOX_ERROR("Should not be here.");
      }

   bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_index) const;

public:

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a VisIt plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& writer);
#endif

private:

   void tagShells(
      pdat::CellData<int> &tag_data,
      const geom::CartesianPatchGeometry &patch_geom,
      const std::vector<double> &buffer_distance ) const;

   std::string d_name;

   const tbox::Dimension d_dim;

   /*!
    * @brief PatchHierarchy for use in implementations of some
    * abstract interfaces that do not specify a hierarch.
    */
   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Radii of shells.
    */
   double d_center[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Radii of shells.
    */
   std::vector<double> d_radii;

   /*!
    * @brief Buffer distances for generating tags.
    */
   std::vector<std::vector<double> >  d_buffer_distance;

   /*!
    * @brief Whether to allocate data on the mesh.
    */
   bool d_allocate_data;

};

#endif  // included_SphericalShellGenerator
