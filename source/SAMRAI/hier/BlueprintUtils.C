/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2017 Lawrence Livermore National Security, LLC
 * Description:   BlueprintUtilities
 *
 ************************************************************************/
#include "SAMRAI/hier/BlueprintUtils.h"

#ifdef HAVE_CONDUIT
#include "SAMRAI/hier/PatchHierarchy.h"

#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

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
 ***************************************************************************
 */

void BlueprintUtils::putTopologyAndCoordinatesToDatabase(
   const std::shared_ptr<tbox::Database>& database,
   const PatchHierarchy& hierarchy) const
{
   TBOX_ASSERT(database);

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
         const BoxId& box_id = patch->getBox().getBoxId();
         const LocalId& local_id = box_id.getLocalId();

         int domain_id = first_patch_id[i] + local_id.getValue();
         std::string domain_name =
            "domain_" + tbox::Utilities::intToString(domain_id, 6);

         std::shared_ptr<tbox::Database> domain_db;
         if (database->isDatabase(domain_name)) {
            domain_db = database->getDatabase(domain_name);
         } else {
            domain_db = database->putDatabase(domain_name);
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

         std::shared_ptr<tbox::Database> topologies_db;
         if (domain_db->isDatabase("topologies")) {
            topologies_db = domain_db->getDatabase("topologies");
         } else {
            topologies_db = domain_db->putDatabase("topologies");
         }

         std::shared_ptr<tbox::Database> mesh_db;
         if (topologies_db->isDatabase("mesh")) {
            mesh_db = topologies_db->getDatabase("mesh");
         } else {
            mesh_db = topologies_db->putDatabase("mesh");
         }

         d_strategy->putCoordinatesToDatabase(
            coords_db, *patch);

         mesh_db->putString("coordset", "coords");

         std::string coords_type = coords_db->getString("type");
         if (coords_type == "explicit") {
            std::shared_ptr<tbox::Database> elements_db;
            if (mesh_db->isDatabase("elements")) {
               elements_db = mesh_db->getDatabase("elements");
            } else {
               elements_db = mesh_db->putDatabase("elements");
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

            mesh_db->putString("type", "structured");
         } else {
            mesh_db->putString("type", coords_db->getString("type"));
         }
      }
   }
}

void BlueprintUtils::writeBlueprintMesh(
   const conduit::Node& reference,
   const tbox::SAMRAI_MPI& samrai_mpi,
   const int num_global_domains,
   const std::string& mesh_name,
   const std::string& data_name,
   const std::string& rootfile_name,
   const std::string& io_protocol) const
{
   conduit::Node index;
   int my_rank = samrai_mpi.getRank();

   if (my_rank == 0) {
      conduit::Node &bpindex = index["blueprint_index"];
         conduit::blueprint::mesh::generate_index(
            reference["domain_000000"], "",
            num_global_domains, bpindex[mesh_name]);

      std::string file_pattern = data_name + "%06d." + io_protocol;
      index["protocol/name"].set(io_protocol);
      index["protocol/version"].set(CONDUIT_VERSION);

      index["number_of_files"].set(samrai_mpi.getSize());
      index["number_of_trees"].set(num_global_domains);
      index["file_pattern"].set(file_pattern);
      index["tree_pattern"].set("/domain_%06d");

      index.save(rootfile_name, io_protocol);
   }

   std::string file_name = data_name + tbox::Utilities::intToString(my_rank, 6) + "." + io_protocol;

   reference.save(file_name, io_protocol);
}

}
}

#endif // HAVE_CONDUIT
