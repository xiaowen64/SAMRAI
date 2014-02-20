/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Class for finding multiblock singularities
 *
 ************************************************************************/

#ifndef included_hier_SingularityFinder
#define included_hier_SingularityFinder

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"

#include "boost/make_shared.hpp"


namespace SAMRAI {
namespace hier {


/*!
 */

class SingularityFinder
{
public:

   SingularityFinder(const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   ~SingularityFinder();

   void findSingularities(
      std::vector<std::set<int> >& singularity_blocks,
      BaseGridGeometry* grid_geometry,
      const std::map< BlockId, std::set<BlockId> > face_neighbors);
 

private:

   struct Block;
   struct Face;
   struct Edge;
   struct Point;

   struct Block {

      Block(const tbox::Dimension& dim) {

         if (dim.getValue() == 1) {
            d_nfaces = 2;
            d_nedges = 0;
            d_npoints = 0;
         } else if (dim.getValue() == 2) {
            d_nfaces = 4;
            d_nedges = 0;
            d_npoints = 4;
         } else if (dim.getValue() == 3) { 
            d_nfaces = 6;
            d_nedges = 12;
            d_npoints = 8;
         } else {
            d_nfaces = 0;
            d_nedges = 0;
            d_npoints = 0;
         }

         if (d_nfaces) {
            d_face.resize(d_nfaces);
         }
         if (d_nedges) {
            d_edge.resize(d_nedges);
         }
         if (d_npoints) {
            d_point.resize(d_npoints);
         }

      }

      int d_nfaces;
      int d_nedges;
      int d_npoints;
      std::vector<boost::shared_ptr<Face> > d_face;
      std::vector<boost::shared_ptr<Edge> > d_edge;
      std::vector<boost::shared_ptr<Point> > d_point;
   };

   struct Face {

      Face() {
         d_bdry = false;
      }

      bool d_bdry;
      std::set<int> d_blocks;
      std::map<int,int> d_block_to_face;
   };

   struct Edge {

      Edge() {
         d_bdry = false;
      }

      bool d_bdry;
      std::set<int> d_blocks;
      std::map<int,int> d_block_to_edge;
   };

   struct Point {

      Point() {
         d_bdry = false;
      }

      bool d_bdry;
      std::set<int> d_blocks;
      std::map<int,int> d_block_to_point;
   };

   SingularityFinder();

   void connect(const BlockId& block_a,
                const BlockId& block_b,
                const BaseGridGeometry* grid_geometry);

   void analyzeConnections();

   void findCoincidentEdges(std::map<int,int>& map_of_edges,
                            const BlockId& block_a,
                            const BlockId& block_b,
                            int facea,
                            const BaseGridGeometry* grid_geometry);

   void findCoincidentPoints(std::map<int,int>& map_of_points,
                             const BlockId& block_a,
                             const BlockId& block_b,
                             int facea,
                             const BaseGridGeometry* grid_geometry);

   tbox::Dimension d_dim;

   std::vector<boost::shared_ptr<Block> > d_blocks;
   std::vector<boost::shared_ptr<Face> > d_faces;
   std::vector<boost::shared_ptr<Edge> > d_edges;
   std::vector<boost::shared_ptr<Point> > d_points;

   static std::vector< std::vector<int> > s_face_edges;
   static std::vector< std::vector<int> > s_face_nodes;



/*
int face_edges[6][4] = {{ 0, 2, 4, 6}, // face 0
                        { 1, 3, 5, 7}, // face 1
                        { 0, 1, 8,10}, // face 2
                        { 2, 3, 9,11}, // face 3
                        { 4, 5, 8, 9}, // face 4
                        { 6, 7,10,11}  // face 5
};

int face_points[6][4] = {{ 0, 2, 4, 6}, // face 0
                         { 1, 3, 5, 7}, // face 1
                         { 0, 1, 4, 5}, // face 2
                         { 2, 3, 6, 7}, // face 3
                         { 0, 1, 2, 3}, // face 4
                         { 4, 5, 6, 7}  // face 5
*/
};

}
}

#endif
