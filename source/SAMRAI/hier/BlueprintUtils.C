/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2020 Lawrence Livermore National Security, LLC
 * Description:   Blueprint utilities
 *
 ************************************************************************/
#include "SAMRAI/hier/BlueprintUtils.h"

#ifdef SAMRAI_HAVE_CONDUIT
#include "SAMRAI/hier/PatchHierarchy.h"

#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io_hdf5.hpp"

namespace SAMRAI {
namespace hier {

/*
 * Constructor does nothing because the objects are stateless.
 */

BlueprintUtils::BlueprintUtils(BlueprintUtilsStrategy* strategy)
: d_strategy(strategy)
{
}

BlueprintUtils::~BlueprintUtils()
{
}

/*
 ***************************************************************************
 *
 * Put topology/coordinate data into a database
 *
 ***************************************************************************
 */
void BlueprintUtils::putTopologyAndCoordinatesToDatabase(
   const std::shared_ptr<tbox::Database>& blueprint_db,
   const PatchHierarchy& hierarchy,
   const std::string& topology_name) const
{
   TBOX_ASSERT(blueprint_db);

   std::vector<int> first_patch_id;
   first_patch_id.push_back(0);

   int patch_count = 0;
   for (int i = 1; i < hierarchy.getNumberOfLevels(); ++i) {
      patch_count += hierarchy.getPatchLevel(i-1)->getNumberOfPatches();
      first_patch_id.push_back(patch_count);
   }

   for (int i = 0; i < hierarchy.getNumberOfLevels(); ++i) {
      const std::shared_ptr<PatchLevel>& level(
         hierarchy.getPatchLevel(i));

      for (PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {

         const std::shared_ptr<Patch>& patch = *p;
         const Box& patch_box = patch->getBox();
         const BoxId& box_id = patch_box.getBoxId();
         const LocalId& local_id = box_id.getLocalId();

         int domain_id = first_patch_id[i] + local_id.getValue();
         std::string domain_name =
            "domain_" + tbox::Utilities::intToString(domain_id, 6);

         std::shared_ptr<tbox::Database> domain_db;
         if (blueprint_db->isDatabase(domain_name)) {
            domain_db = blueprint_db->getDatabase(domain_name);
         } else {
            domain_db = blueprint_db->putDatabase(domain_name);
         }

         std::shared_ptr<tbox::Database> coordsets_db;
         if (domain_db->isDatabase("coordsets")) {
            coordsets_db = domain_db->getDatabase("coordsets");
         } else {
            coordsets_db = domain_db->putDatabase("coordsets");
         }

         std::shared_ptr<tbox::Database> coords_db;
         if (coordsets_db->isDatabase("coords")) {
            coords_db = coordsets_db->getDatabase("coords");
         } else {
            coords_db = coordsets_db->putDatabase("coords");
         }

         std::shared_ptr<tbox::Database> topologies_db(
            domain_db->putDatabase("topologies"));

         std::shared_ptr<tbox::Database> topo_db(
            topologies_db->putDatabase(topology_name));

         std::shared_ptr<tbox::Database> elem_db(
            topo_db->putDatabase("elements"));
         std::shared_ptr<tbox::Database> origin_db(
            elem_db->putDatabase("origin"));
         origin_db->putInteger("i0", patch_box.lower(0));
         if (patch_box.getDim().getValue() > 1) {
            origin_db->putInteger("j0", patch_box.lower(1));
         }
         if (patch_box.getDim().getValue() > 2) {
            origin_db->putInteger("k0", patch_box.lower(2));
         }

         if (d_strategy) {
            d_strategy->putCoordinatesToDatabase(
               coords_db, *patch);
         }

         topo_db->putString("coordset", "coords");

         std::string coords_type = coords_db->getString("type");
         if (coords_type == "explicit") {
            std::shared_ptr<tbox::Database> elements_db;
            if (topo_db->isDatabase("elements")) {
               elements_db = topo_db->getDatabase("elements");
            } else {
               elements_db = topo_db->putDatabase("elements");
            }

            std::shared_ptr<tbox::Database> dims_db;
            if (elements_db->isDatabase("dims")) {
               dims_db = elements_db->getDatabase("dims");
            } else {
               dims_db = elements_db->putDatabase("dims");
            }

            const tbox::Dimension& dim = patch->getDim();
            dims_db->putInteger("i", patch->getBox().numberCells(0));
            if (dim.getValue() > 1) { 
               dims_db->putInteger("j", patch->getBox().numberCells(1));
            }
            if (dim.getValue() > 2) { 
               dims_db->putInteger("k", patch->getBox().numberCells(2));
            }

            topo_db->putString("type", "structured");
         } else {
            topo_db->putString("type", coords_type);
         }
      }
   }
}

/*
 ***************************************************************************
 *
 * Write a blueprint to files
 *
 ***************************************************************************
 */
void BlueprintUtils::writeBlueprintMesh(
   const conduit::Node& blueprint,
   const tbox::SAMRAI_MPI& samrai_mpi,
   const int num_global_domains,
   const std::string& mesh_name,
   const std::string& data_name,
   const std::string& rootfile_name,
   const std::string& io_protocol) const
{
   conduit::Node index;
   int my_rank = samrai_mpi.getRank();


   std::vector<int> partition_map(num_global_domains, 0);
   conduit::NodeConstIterator itr = blueprint.children();
   while (itr.has_next()) {
      const conduit::Node &child = itr.next();
      if (child.has_child("state")) {
         const conduit::Node &state = child["state"];
         if (state.has_child("domain_id")) {
            conduit::int64 dom_id = state["domain_id"].as_int64();
            partition_map[dom_id] = my_rank;
         }
      }
   }

   samrai_mpi.AllReduce(&partition_map[0], num_global_domains, MPI_MAX);

   if (my_rank == 0) {
      conduit::Node &bpindex = index["blueprint_index"];
      conduit::blueprint::mesh::generate_index(
         blueprint["domain_000000"], "",
         num_global_domains, bpindex[mesh_name]);

      std::string file_pattern = data_name + "%06d." + io_protocol;
      index["protocol/name"].set(io_protocol);
      index["protocol/version"].set(CONDUIT_VERSION);

#if 1
      index["number_of_files"].set(samrai_mpi.getSize());
      index["number_of_trees"].set(num_global_domains);
      index["file_pattern"].set(file_pattern);
      index["tree_pattern"].set("/domain_%06d");
#else
      //This is the new Blueprint index format that isn't official yet.

      conduit::Node &mesh = bpindex[mesh_name]; 
      mesh["state/partition_size"].set(samrai_mpi.getSize());
      mesh["state/number_of_domains"].set(num_global_domains);
      mesh["state/partition_pattern"].set(file_pattern);
      mesh["state/domain_pattern"].set("/domain_%06d");
      mesh["state/domain_to_partition_map"].set(partition_map);
#endif

      if (io_protocol != "hdf5") {
         index.save(rootfile_name, io_protocol);
      } else {
         conduit::relay::io::hdf5_save(index, rootfile_name);
      }
   }

   std::string file_name = data_name + tbox::Utilities::intToString(my_rank, 6) + "." + io_protocol;

   if (io_protocol != "hdf5") {
      blueprint.save(file_name, io_protocol);
   } else {
      conduit::relay::io::hdf5_save(blueprint, file_name);
   }
}

}
}

#endif // SAMRAI_HAVE_CONDUIT
