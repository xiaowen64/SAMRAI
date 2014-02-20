/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Class for finding multiblockSingularities
 *
 ************************************************************************/
#include "SAMRAI/hier/SingularityFinder.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace hier {

std::vector< std::vector<int> > SingularityFinder::s_face_edges;
std::vector< std::vector<int> > SingularityFinder::s_face_nodes;

/*
 * ************************************************************************
 *
 * Constructors
 *
 * ************************************************************************
 */

SingularityFinder::SingularityFinder(const tbox::Dimension& dim)
: d_dim(dim)
{

   if (d_dim.getValue() == 3 && s_face_edges.empty()) {
      std::vector<int> edges;
      edges.resize(4);
      edges[0] = 0;
      edges[1] = 2;
      edges[2] = 4;
      edges[3] = 6;
      s_face_edges.push_back(edges);
      edges[0] = 1;
      edges[1] = 3;
      edges[2] = 5;
      edges[3] = 7;
      s_face_edges.push_back(edges);
      edges[0] = 0;
      edges[1] = 1;
      edges[2] = 8;
      edges[3] = 10;
      s_face_edges.push_back(edges);
      edges[0] = 2;
      edges[1] = 3;
      edges[2] = 9;
      edges[3] = 11;
      s_face_edges.push_back(edges);
      edges[0] = 4;
      edges[1] = 5;
      edges[2] = 8;
      edges[3] = 9;
      s_face_edges.push_back(edges);
      edges[0] = 6;
      edges[1] = 7;
      edges[2] = 10;
      edges[3] = 11;
      s_face_edges.push_back(edges);
   }

   if (s_face_nodes.empty()) {
      std::vector<int> nodes;
      if (d_dim.getValue() == 3) {
         nodes.resize(4);
         nodes[0] = 0;
         nodes[1] = 2;
         nodes[2] = 4;
         nodes[3] = 6;
         s_face_nodes.push_back(nodes);
         nodes[0] = 1;
         nodes[1] = 3;
         nodes[2] = 5;
         nodes[3] = 7;
         s_face_nodes.push_back(nodes);
         nodes[0] = 0;
         nodes[1] = 1;
         nodes[2] = 4;
         nodes[3] = 5;
         s_face_nodes.push_back(nodes);
         nodes[0] = 2;
         nodes[1] = 3;
         nodes[2] = 6;
         nodes[3] = 7;
         s_face_nodes.push_back(nodes);
         nodes[0] = 0;
         nodes[1] = 1;
         nodes[2] = 2;
         nodes[3] = 3;
         s_face_nodes.push_back(nodes);
         nodes[0] = 4;
         nodes[1] = 5;
         nodes[2] = 6;
         nodes[3] = 7;
         s_face_nodes.push_back(nodes);
      } else if (d_dim.getValue() == 2) {
         nodes.resize(2);
         nodes[0] = 0;
         nodes[1] = 2;
         s_face_nodes.push_back(nodes);
         nodes[0] = 1;
         nodes[1] = 3;
         s_face_nodes.push_back(nodes);
         nodes[0] = 0;
         nodes[1] = 1;
         s_face_nodes.push_back(nodes);
         nodes[0] = 2;
         nodes[1] = 3;
         s_face_nodes.push_back(nodes);
      }
   }
}

/*
 * ************************************************************************
 *
 * Destructor
 *
 * ************************************************************************
 */

SingularityFinder::~SingularityFinder()
{
}

void
SingularityFinder::findSingularities(
   std::vector<std::set<int> >& singularity_blocks,
   BaseGridGeometry* grid_geometry,
   const std::map< BlockId, std::set<BlockId> > face_neighbors)
{
   if (face_neighbors.empty()) return;

   std::map< BlockId, std::set<BlockId> > unprocessed = face_neighbors;
   std::set<BlockId> can_be_processed;
   can_be_processed.insert(face_neighbors.begin()->first);

   do {
      std::map<BlockId, std::set<BlockId> >::iterator un_itr = unprocessed.begin();
      BlockId base_block = un_itr->first;

      while (can_be_processed.find(base_block) == can_be_processed.end()) {
         ++un_itr;
         base_block = un_itr->first;
      }

      const std::set<BlockId>& nbr_blocks = face_neighbors.find(base_block)->second;
      for (std::set<BlockId>::const_iterator nbr_itr = nbr_blocks.begin();
           nbr_itr != nbr_blocks.end(); ++nbr_itr) {

         const BlockId& nbr_block = *nbr_itr;

         if (unprocessed[base_block].find(nbr_block) != unprocessed[base_block].end()) {
            connect(base_block, nbr_block, grid_geometry);

            can_be_processed.insert(nbr_block); 

            unprocessed[base_block].erase(nbr_block);
            if (unprocessed[base_block].empty()) {
               unprocessed.erase(base_block);
            }
            unprocessed[nbr_block].erase(base_block);
            if (unprocessed[nbr_block].empty()) {
               unprocessed.erase(nbr_block);
            }
         }
      }
   } while (!unprocessed.empty());

   analyzeConnections();

   for (int i = 0; i < d_edges.size(); i++) {
      if (d_edges[i]->d_bdry) continue;
      int nblocks = d_edges[i]->d_blocks.size();
      bool enhanced = false;
      bool reduced = false;
      if (nblocks < 4) {
         reduced = true;
      } else if (nblocks > 4) {
         enhanced = true;
      }

      if (reduced || enhanced) {
         singularity_blocks.push_back(d_edges[i]->d_blocks); 
      }
   }

   for (int i = 0; i < d_points.size(); i++) {
      if (d_points[i]->d_bdry) continue;
      int nblocks = d_points[i]->d_blocks.size();
      bool enhanced = false;
      bool reduced = false;
      if (nblocks < (1<<d_dim.getValue())) {
         reduced = true;
      } else if (nblocks > (1<<d_dim.getValue())) {
         enhanced = true;
      }

      if (reduced || enhanced) {
         singularity_blocks.push_back(d_points[i]->d_blocks);
      }
   }

}

void
SingularityFinder::connect(const BlockId& block_a,
                           const BlockId& block_b,
                           const BaseGridGeometry* grid_geometry)
{

   int facea = -1;
   int faceb = -1; 

   const std::map<BlockId,BaseGridGeometry::Neighbor>& nbrs_of_a =
      grid_geometry->getNeighbors(block_a);

   std::map<BlockId,BaseGridGeometry::Neighbor>::const_iterator itr = nbrs_of_a.find(block_b);
   if (itr != nbrs_of_a.end()) {
      const BaseGridGeometry::Neighbor& neighbor = itr->second;
      TBOX_ASSERT(neighbor.getBlockId() == block_b); 

      BoxContainer a_domain(grid_geometry->getPhysicalDomain(), block_a);
      Box node_box_a(a_domain.getBoundingBox());
      node_box_a.upper() += IntVector::getOne(d_dim);

      Box node_box_b(neighbor.getTransformedDomain().getBoundingBox());
      node_box_b.upper() += IntVector::getOne(d_dim);

      Box shared_face(node_box_a * node_box_b);

      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (shared_face.lower()(d) == shared_face.upper()(d)) {
            if (shared_face.lower(d) == node_box_a.lower(d)) {
               facea = 2*d;
            } else {
               facea = 2*d + 1;
            }
         }
      } 
   }

   const std::map<BlockId,BaseGridGeometry::Neighbor>& nbrs_of_b =
      grid_geometry->getNeighbors(block_b);

   itr = nbrs_of_b.find(block_a);
   if (itr != nbrs_of_b.end()) {
      const BaseGridGeometry::Neighbor& neighbor = itr->second;
      TBOX_ASSERT(neighbor.getBlockId() == block_a);

      BoxContainer b_domain(grid_geometry->getPhysicalDomain(), block_b);
      Box node_box_b(b_domain.getBoundingBox());
      node_box_b.upper() += IntVector::getOne(d_dim);

      Box node_box_a(neighbor.getTransformedDomain().getBoundingBox());
      node_box_a.upper() += IntVector::getOne(d_dim);

      Box shared_face(node_box_a * node_box_b);

      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (shared_face.lower()(d) == shared_face.upper()(d)) {
            if (shared_face.lower(d) == node_box_b.lower(d)) {
               faceb = 2*d;
            } else {
               faceb = 2*d + 1;
            }
         }
      }
   }

   boost::shared_ptr<Face> face = boost::make_shared<Face>();
   d_faces.push_back(face);

   if (d_blocks.empty()) {
      d_blocks.resize(grid_geometry->getNumberBlocks());
   }

   int a = block_a.getBlockValue();
   int b = block_b.getBlockValue();
   if (d_blocks[a].get() == 0) {
      d_blocks[a] = boost::make_shared<Block>(d_dim);
   }
   if (d_blocks[b].get() == 0) {
      d_blocks[b] = boost::make_shared<Block>(d_dim);
   }

   d_blocks[a]->d_face[facea] = face;
   d_blocks[b]->d_face[faceb] = face;
   face->d_block_to_face[a] = facea; 
   face->d_block_to_face[b] = faceb;

   int nedges_per_face = d_dim.getValue() == 3 ? 4 : 0;

   /*
    * Map from edge on 'a' to edge on 'b'
    */ 
   std::map<int,int> map_of_edges;
   findCoincidentEdges(map_of_edges, block_a, block_b, facea, grid_geometry);

   for (std::map<int,int>::const_iterator e_itr = map_of_edges.begin();
        e_itr != map_of_edges.end(); ++e_itr) {

      boost::shared_ptr<Edge>& edgea = d_blocks[a]->d_edge[ e_itr->first ];
      boost::shared_ptr<Edge>& edgeb = d_blocks[b]->d_edge[ e_itr->second ];

      if (edgea.get() == 0 && edgeb.get() == 0) {
         edgea.reset(new Edge());
         edgeb = edgea;
         d_edges.push_back(edgea);
      } else if (edgea.get() != 0) {
         edgeb = edgea;
      }
      else {
         edgea = edgeb;
      }

      edgea->d_blocks.insert(a);
      edgea->d_blocks.insert(b);
      edgea->d_block_to_edge[a] = e_itr->first;
      edgea->d_block_to_edge[b] = e_itr->second;
   }

   std::map<int,int> map_of_points;
   findCoincidentPoints(map_of_points, block_a, block_b, facea, grid_geometry);

   for (std::map<int,int>::const_iterator e_itr = map_of_points.begin();
        e_itr != map_of_points.end(); ++e_itr) {

      boost::shared_ptr<Point>& pointa = d_blocks[a]->d_point[ e_itr->first ];
      boost::shared_ptr<Point>& pointb = d_blocks[b]->d_point[ e_itr->second ];

      if (pointa.get() == 0 && pointb.get() == 0) {
         pointa.reset(new Point());
         pointb = pointa;
         d_points.push_back(pointa);
      } else if (pointa.get() != 0) {
         pointb = pointa;
      }
      else {
         pointa = pointb;
      }

      pointa->d_blocks.insert(a);
      pointa->d_blocks.insert(b);
      pointa->d_block_to_point[a] = e_itr->first;
      pointa->d_block_to_point[b] = e_itr->second;
   }

}

void
SingularityFinder::analyzeConnections()
{
   for (int iblock = 0; iblock < d_blocks.size(); iblock++) {

      boost::shared_ptr<Block>& block = d_blocks[iblock];

      for (int iface = 0; iface < 2*d_dim.getValue(); iface++) {

         boost::shared_ptr<Face>& face = block->d_face[iface];

         if (face.get() == 0) {

            face.reset(new Face());
            d_faces.push_back(face);
            face->d_bdry = true;
            face->d_blocks.insert(iblock);

            if (d_dim.getValue() > 2) {
               for (int iedge = 0; iedge < 4; iedge++) {
                  int edge_idx = s_face_edges[iface][iedge];
                  if (!block->d_edge[edge_idx]) {
                     boost::shared_ptr<Edge> edge = boost::make_shared<Edge>();
                     d_edges.push_back(edge);
                     edge->d_blocks.insert(iblock);
                     block->d_edge[edge_idx] = edge;
                  }
                  block->d_edge[edge_idx]->d_bdry = true;
               }
            }

            int num_pts = 1 << (d_dim.getValue()-1);
            for (int ipoint = 0; ipoint < num_pts; ++ipoint) {
               int point_idx = s_face_nodes[iface][ipoint];
               if (!block->d_point[point_idx]) {
                  boost::shared_ptr<Point> point = boost::make_shared<Point>();
                  point->d_blocks.insert(iblock);
                  block->d_point[point_idx] = point;
               }
               block->d_point[point_idx]->d_bdry = true;
            }
         }

      }
   }

}

void
SingularityFinder::findCoincidentEdges(
   std::map<int,int>& map_of_edges,
   const BlockId& block_a,
   const BlockId& block_b,
   int facea,
   const BaseGridGeometry* grid_geometry)
{
   if (d_dim.getValue() != 3) return;

   BoxContainer a_domain(grid_geometry->getPhysicalDomain(), block_a);
   Box a_box(a_domain.getBoundingBox());

   BoxContainer b_domain(grid_geometry->getPhysicalDomain(), block_b);
   Box b_node_box(b_domain.getBoundingBox());
   b_node_box.upper() += IntVector::getOne(d_dim);
   IntVector b_box_size(b_node_box.numberCells());

   int nedges_per_face = 4;

   for (int i = 0; i < nedges_per_face; i++) {

      int edgea_idx = s_face_edges[facea][i];
      int edgeb_idx = -1;

      Box edge_box(a_box);

      switch (edgea_idx) {

         case 0:
            edge_box.upper()(0) = edge_box.lower(0);
            edge_box.upper()(1) = edge_box.lower(1);
            break;
         case 1:
            edge_box.lower()(0) = edge_box.upper(0);
            edge_box.upper()(1) = edge_box.lower(1);
            break;
         case 2:
            edge_box.upper()(0) = edge_box.lower(0);
            edge_box.lower()(1) = edge_box.upper(1);
            break;
         case 3:
            edge_box.lower()(0) = edge_box.upper(0);
            edge_box.lower()(1) = edge_box.upper(1);
            break;
         case 4:
            edge_box.upper()(0) = edge_box.lower(0);
            edge_box.upper()(2) = edge_box.lower(2);
            break;
         case 5:
            edge_box.lower()(0) = edge_box.upper(0);
            edge_box.upper()(2) = edge_box.lower(2);
            break;
         case 6:
            edge_box.upper()(0) = edge_box.lower(0);
            edge_box.lower()(2) = edge_box.upper(2);
            break;
         case 7:
            edge_box.lower()(0) = edge_box.upper(0);
            edge_box.lower()(2) = edge_box.upper(2);
            break;
         case 8:
            edge_box.upper()(1) = edge_box.lower(1);
            edge_box.upper()(2) = edge_box.lower(2);
            break;
         case 9:
            edge_box.lower()(1) = edge_box.upper(1);
            edge_box.upper()(2) = edge_box.lower(2);
            break;
         case 10:
            edge_box.upper()(1) = edge_box.lower(1);
            edge_box.lower()(2) = edge_box.upper(2);
            break;
         case 11:
            edge_box.lower()(1) = edge_box.upper(1);
            edge_box.lower()(2) = edge_box.upper(2);
            break;
         default:
            break;
      }

      bool transformed = grid_geometry->transformBox(edge_box,
                                                     IntVector::getOne(d_dim),
                                                     block_b,
                                                     block_a);
      TBOX_ASSERT(transformed);
      edge_box.upper() += IntVector::getOne(d_dim);
      Box b_edge(edge_box*b_node_box);
      
      IntVector b_edge_dirs(b_edge.numberCells());
      for (int d = 0; d < d_dim.getValue(); ++d) {
         TBOX_ASSERT(b_edge_dirs[d] >= 1);
         if (b_edge_dirs[d] == b_box_size[d]) {
            b_edge_dirs[d] = 0;
         } else if (b_edge.lower()(d) == b_node_box.lower(d)) {
            b_edge_dirs[d] = -1;
         } else if (b_edge.upper()(d) == b_node_box.upper(d)) {
            b_edge_dirs[d] = 1;
         } else {
            TBOX_ERROR("   ");
         }
      }

      if (b_edge_dirs[0] == 0) {
         TBOX_ASSERT(b_edge_dirs[1] != 0 && b_edge_dirs[2] != 0);
         if (b_edge_dirs[1] == -1) {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 8;       
            } else {
               edgeb_idx = 10;
            }
         } else {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 9;
            } else {
               edgeb_idx = 11;
            }
         }
      } else if (b_edge_dirs[1] == 0) {
         TBOX_ASSERT(b_edge_dirs[0] != 0 && b_edge_dirs[2] != 0);
         if (b_edge_dirs[0] == -1) {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 4;
            } else {
               edgeb_idx = 6;
            }
         } else {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 5;
            } else {
               edgeb_idx = 7;
            }
         }
      } else if (b_edge_dirs[2] == 0) {
         TBOX_ASSERT(b_edge_dirs[0] != 0 && b_edge_dirs[1] != 0);
         if (b_edge_dirs[0] == -1) {
            if (b_edge_dirs[1] == -1) {
               edgeb_idx = 0;
            } else {
               edgeb_idx = 2;
            }
         } else {
            if (b_edge_dirs[1] == -1) {
               edgeb_idx = 1;
            } else {
               edgeb_idx = 3;
            }
         }
      } else {
         TBOX_ERROR("  "); 
      }

      TBOX_ASSERT(edgeb_idx >= 0);

      map_of_edges[edgea_idx] = edgeb_idx; 

   }
   

}


void
SingularityFinder::findCoincidentPoints(
   std::map<int,int>& map_of_points,
   const BlockId& block_a,
   const BlockId& block_b,
   int facea,
   const BaseGridGeometry* grid_geometry)
{
   BoxContainer a_domain(grid_geometry->getPhysicalDomain(), block_a);
   Box a_box(a_domain.getBoundingBox());

   BoxContainer b_domain(grid_geometry->getPhysicalDomain(), block_b);
   Box b_node_box(b_domain.getBoundingBox());
   b_node_box.upper() += IntVector::getOne(d_dim);
   IntVector b_box_size(b_node_box.numberCells());

   int npoints_per_face = 1 << (d_dim.getValue()-1); 

   for (int i = 0; i < npoints_per_face; i++) {

      int pointa_idx = s_face_nodes[facea][i];

      Box point_box(a_box);

      IntVector corner_dirs(d_dim, 0);
      corner_dirs[0] = (pointa_idx % 2 == 0) ? -1 : 1;
      if (d_dim.getValue() > 1) {
         corner_dirs[1] = ((pointa_idx/2) % 2 == 0) ? -1 : 1;
      }
      if (d_dim.getValue() > 2) {
         corner_dirs[2] = ((pointa_idx/4) % 2 == 0) ? -1 : 1;
      }

      for (int d = 0; d < d_dim.getValue(); ++d) {
         TBOX_ASSERT(corner_dirs[d] == -1 || corner_dirs[d] == 1);
         if (corner_dirs[d] == -1) {
            point_box.upper()(d) = point_box.lower(d);
         } else {
            point_box.lower()(d) = point_box.upper(d);
         }
      }

      bool transformed = grid_geometry->transformBox(point_box,
                                                     IntVector::getOne(d_dim),
                                                     block_b,
                                                     block_a);
      TBOX_ASSERT(transformed);
      point_box.upper() += IntVector::getOne(d_dim);
      Box b_point(point_box*b_node_box);

      IntVector b_point_dirs(d_dim, 0);
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (b_point.lower()(d) == b_node_box.lower()(d)) {
            b_point_dirs[d] = -1;
         } else if (b_point.upper()(d) == b_node_box.upper()(d)) {
            b_point_dirs[d] = 1;
         } else {
            TBOX_ERROR("   ");
         }
      }

      int pointb_idx = 0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (b_point_dirs[d] == 1) {
            pointb_idx += (1 << d);
         }
      }

      TBOX_ASSERT(pointb_idx >= 0);

      map_of_points[pointa_idx] = pointb_idx; 

   }
   

}




}
}
