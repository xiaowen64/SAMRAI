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
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay.hpp"

namespace SAMRAI {
namespace hier {

/*
 * Constructor does nothing because the objects are stateless.
 */

BlueprintUtils::BlueprintUtils(BlueprintUtilsStrategy* strategy, bool do_uniform)
: d_strategy(strategy),
  d_do_uniform(do_uniform)
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

         std::shared_ptr<tbox::Database> dims_db(
            elem_db->putDatabase("dims"));

         const tbox::Dimension& dim = patch->getDim();
         dims_db->putInteger("i", patch->getBox().numberCells(0));
         if (dim.getValue() > 1) {
            dims_db->putInteger("j", patch->getBox().numberCells(1));
         }
         if (dim.getValue() > 2) {
            dims_db->putInteger("k", patch->getBox().numberCells(2));
         }




         if (d_strategy) {
            d_strategy->putCoordinatesToDatabase(
               coords_db, *patch, patch->getBox(), d_do_uniform);
         }

         topo_db->putString("coordset", "coords");

         std::string coords_type = coords_db->getString("type");
         if (coords_type == "explicit") {
            topo_db->putString("type", "structured");
         } else {
            topo_db->putString("type", coords_type);
         }
      }
   }
}

void BlueprintUtils::setFlattenedCoordset(
   std::shared_ptr<tbox::Database>& coords_db,
   const std::shared_ptr<tbox::Database>& topo_db,
   const Patch& patch,
   const Box& box) const
{
         if (d_strategy) {
            //Flattened coordest has to be explicit because of coarse-fine shared nodes
            d_strategy->putCoordinatesToDatabase(
               coords_db, patch, box, false);
         }

         topo_db->putString("coordset", "coords");

         std::string coords_type = coords_db->getString("type");
         if (coords_type == "explicit") {
            topo_db->putString("type", "structured");
         } else {
            topo_db->putString("type", coords_type);
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

   conduit::Node &bpindex = index["blueprint_index"];
   conduit::blueprint::mpi::mesh::generate_index(
            blueprint, "",
            bpindex[mesh_name],
            samrai_mpi.getCommunicator());

   if (my_rank == 0) {
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

   blueprint.save(file_name, io_protocol);
}

void BlueprintUtils::flattenFields(conduit::Node& flat_node,
                   const conduit::Node& amr_node)
{

   conduit::NodeIterator itr = flat_node.children();

   while(itr.has_next())
   {
      conduit::Node& domain = itr.next();
//      std::string domain_name = itr.name();

      int64_t amr_domain_id = domain["state/amr_domain_id"].as_int64();

      std::ostringstream dom_oss;
      dom_oss << "domain_" << std::setw(6) << std::setfill('0') << amr_domain_id;
      std::string amr_domain_name = dom_oss.str();

      const conduit::Node& amr_domain = amr_node[amr_domain_name];

      conduit::Node& fields = domain["fields"];
      const conduit::Node& amr_fields = amr_domain["fields"];

      conduit::NodeConstIterator flds_itr = amr_fields.children();

      while(flds_itr.has_next())
      {
         const conduit::Node& amr_field = flds_itr.next();
         std::string fld_name = flds_itr.name();

         conduit::NodeConstIterator fld_itr = amr_field.children();
         while(fld_itr.has_next())
         {
            const conduit::Node& subfld = fld_itr.next();
            std::string sub_name = fld_itr.name();

            if (sub_name != "values") {
               fields[fld_name][sub_name] = subfld;
            }
         }

         std::string topo_name = fields[fld_name]["topology"].as_string();
         const conduit::Node& topo = domain["topologies"][topo_name];
         const conduit::Node& amr_topo = amr_domain["topologies"][topo_name];

         int64_t amr_isize = amr_topo["elements/dims/i"].as_int64();
         int64_t amr_jsize = amr_topo["elements/dims/j"].as_int64();
         int64_t amr_ksize = 1;
         int64_t dom_isize = topo["elements/dims/i"].as_int64();
         int64_t dom_jsize = topo["elements/dims/j"].as_int64();
         int64_t dom_ksize = 1;
         if (amr_topo.has_child("elements/dims/k")) {
            amr_ksize = amr_topo["elements/dims/k"].as_int64();
            dom_ksize = topo["elements/dims/k"].as_int64(); 
         }

         int64_t array_size = dom_isize*dom_jsize*dom_ksize;
         std::vector<double> fld_vec(array_size);

         int64_t amr_i0 = amr_topo["elements/origin/i0"].as_int64();
         int64_t amr_j0 = amr_topo["elements/origin/j0"].as_int64();
         int64_t amr_k0 = 1;
         int64_t dom_i0 = topo["elements/origin/i0"].as_int64();
         int64_t dom_j0 = topo["elements/origin/j0"].as_int64();
         int64_t dom_k0 = 1;
         if (amr_topo.has_child("elements/origin/k0")) {
            amr_k0 = amr_topo["elements/origin/k0"].as_int64();
            dom_k0 = topo["elements/origin/k0"].as_int64();
         }

         conduit::double_array amr_vals = amr_field["values"].as_double_array();

         int64_t ioffset = dom_i0 - amr_i0;
         int64_t joffset = dom_j0 - amr_j0;
         int64_t koffset = dom_k0 - amr_k0;

         for (int64_t k = 0; k < dom_ksize; ++k) {
            for (int64_t j = 0; j < dom_jsize; ++j) {
               for (int64_t i = 0; i < dom_isize; ++i) {

                  int64_t dom_offset = i + j*dom_isize + k*dom_isize*dom_ksize; 
                  int64_t amr_offset = (i+ioffset) +
                                       (j+joffset)*amr_isize +
                                       (k+koffset)*amr_isize*amr_jsize;


                  fld_vec[dom_offset] = amr_vals[amr_offset];
               }
            }
         }

         fields[fld_name]["values"].set(fld_vec); 

      }
   }
}


}
}

#endif // SAMRAI_HAVE_CONDUIT
