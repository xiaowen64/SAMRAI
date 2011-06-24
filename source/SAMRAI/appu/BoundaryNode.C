/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Boundary node struct for embedded boundary implementations. 
 *
 ************************************************************************/

#ifndef included_appu_BoundaryNode_C
#define included_appu_BoundaryNode_C

#include "SAMRAI/appu/BoundaryNode.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/appu/EmbeddedBoundaryDefines.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace appu {

const int BoundaryNode::BOUNDARYNODE_VERSION = 1;
const int BoundaryNode::BOUNDARYNODE_LOC_UNDEFINED = -99999;

/*
 *************************************************************************
 *                                                                       *
 * Initialization for static data members.                               *
 *                                                                       *
 *************************************************************************
 */

double BoundaryNode::s_on_boundary_threshold = 1.e-6;

/*
 *************************************************************************
 *                                                                       *
 * Static function to set static data members.                           *
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setOnBoundaryThreshold(
   const double th)
{
   s_on_boundary_threshold = th;
}

/*
 *************************************************************************
 *                                                                       *
 * Default Constructor                                                   *
 *                                                                       *
 *************************************************************************
 */

BoundaryNode::BoundaryNode(
   const tbox::Dimension& dim):
   d_index(dim)
{
   hier::Index dummy(dim, BOUNDARYNODE_LOC_UNDEFINED);
   hier::IntVector loc(dim, 0);
   d_index = pdat::NodeIndex(dummy, loc);
   initializeBoundaryNodeData();
}

/*
 *************************************************************************
 *                                                                       *
 * Construct a boundary cell, given the cell index.
 *                                                                       *
 *************************************************************************
 */

BoundaryNode::BoundaryNode(
   const pdat::NodeIndex& in):
   d_index(in)
{
   initializeBoundaryNodeData();
}

/*
 *************************************************************************
 *                                                                       *
 * Copy Constructor                                                      *
 *                                                                       *
 *************************************************************************
 */

BoundaryNode::BoundaryNode(
   const appu::BoundaryNode& bdry_node):
   d_index(bdry_node.d_index)
{
   copyBoundaryNodeData(bdry_node);
}

/*
 *************************************************************************
 *                                                                       *
 * Assignment operator                                                   *
 *                                                                       *
 *************************************************************************
 */

BoundaryNode& BoundaryNode::operator = (
   const appu::BoundaryNode& bdry_node)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_index, bdry_node.d_index);

   d_index = bdry_node.d_index;
   copyBoundaryNodeData(bdry_node);
   return *this;
}

/*
 *************************************************************************
 *                                                                       *
 * Destructor
 *                                                                       *
 *************************************************************************
 */

BoundaryNode::~BoundaryNode()
{
}

/*
 *************************************************************************
 *
 *  Return cell index
 *                                                                       *
 *************************************************************************
 */
pdat::NodeIndex
BoundaryNode::getIndex() const
{
   return d_index;
}

/*
 *************************************************************************
 *
 *  Return whether the node is on the physical boundary of the shape.
 *                                                                       *
 *************************************************************************
 */
bool
BoundaryNode::getNodeOnBoundary() const
{
   return d_on_boundary;
}

/*
 *************************************************************************
 *
 *  Return number of nearest neighbors
 *                                                                       *
 *************************************************************************
 */
int
BoundaryNode::getNumberOfNearestNeighborNodes() const
{
   return d_num_nearest_neighbors;
}

/*
 *************************************************************************
 *
 *  Return number of nearest neighbors
 *                                                                       *
 *************************************************************************
 */
int
BoundaryNode::getNumberOfOutsideNeighborNodes() const
{
   return d_num_outside_neighbors;
}

/*
 *************************************************************************
 *
 *  Return array of nearest neighbors
 *                                                                       *
 *************************************************************************
 */
tbox::Array<pdat::NodeIndex>
BoundaryNode::getNearestNeighborNodes() const
{
   return d_nearest_neighbors;
}

/*
 *************************************************************************
 *
 *  Return particular nearest neighbor
 *                                                                       *
 *************************************************************************
 */
pdat::NodeIndex
BoundaryNode::getNearestNeighborNode(
   const int i) const
{
   TBOX_ASSERT(i < d_num_nearest_neighbors);

   return d_nearest_neighbors[i];
}

/*
 *************************************************************************
 *
 *  Return location of closest boundary point
 *                                                                       *
 *************************************************************************
 */
const double *
BoundaryNode::getClosestBoundaryPoint() const
{
   return d_closest_boundary_point;
}

double
BoundaryNode::getClosestBoundaryPoint(
   const int i) const
{
   TBOX_ASSERT(i < d_index.getDim().getValue());

   return d_closest_boundary_point[i];
}

/*
 *************************************************************************
 *
 *  Return the distance & normal vector to the embedded boundary
 *                                                                       *
 *************************************************************************
 */

double
BoundaryNode::getDistanceToBoundary() const
{
   if (d_distance_to_boundary < 0.) {
      TBOX_ERROR("BoundaryNode::getDistanceToBoundary()"
         << "\nYou must first set the distance before accessing it."
         << "\nCall 'setDistanceToBoundary(patch)'." << std::endl);
   }
   return d_distance_to_boundary;
}

const double *
BoundaryNode::getNormalToBoundary() const
{
   if (d_distance_to_boundary < 0.) {
      TBOX_ERROR("BoundaryNode::getNormalToBoundary()"
         << "\nYou must first set the normal before accessing it."
         << "\nCall 'setNormalToBoundary(patch)'." << std::endl);
   }
   return d_normal_to_boundary;
}

double
BoundaryNode::getNormalToBoundary(
   const int i) const
{
   if (d_distance_to_boundary < 0.) {
      TBOX_ERROR("BoundaryNode::getNormalToBoundary()"
         << "\nYou must first set the normal before accessing it."
         << "\nEither call 'setNormalToBoundary(patch)' or"
         << "\nuse the method 'getNormalToBoundary(patch)'."
         << std::endl);
   }
   return d_normal_to_boundary[i];
}

/*
 *************************************************************************
 *
 * Set whether the node is on the embedded boundary.
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setNodeOnBoundary()
{
   d_on_boundary = true;
}

/*
 *************************************************************************
 *
 *  Compute the number of outside neighbors to the boundary node.
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setNumOutsideNeighborNodes(
   tbox::Pointer<pdat::NodeData<int> >& node_flag,
   hier::Index& cut_cell_index)
{
   NULL_USE(cut_cell_index);

   TBOX_DIM_ASSERT_CHECK_ARGS3(d_index, *node_flag, cut_cell_index);

   const tbox::Dimension dim(d_index.getDim());

   /*
    * Form a 2-cell box with the boundary node at the center.  Look for
    * neighboring nodes that are OUTSIDE.
    */
   hier::Index lo(dim);
   hier::Index hi(dim);
   for (int i = 0; i < dim.getValue(); i++) {
      lo(i) = d_index(i) - 1;
      hi(i) = d_index(i);
   }
   hier::Box two_cell_bn_box(lo, hi);

   d_num_outside_neighbors = 0;
   for (pdat::NodeIterator on(two_cell_bn_box);
        on; on++) {
      pdat::NodeIndex outside_node = on();
      if ((*node_flag)(outside_node) == EmbeddedBoundaryDefines::OUTSIDE) {
         d_num_outside_neighbors++;
      }
      if (d_num_outside_neighbors >= dim.getValue() * 4) {
         TBOX_ERROR("BoundaryNode::calculateBoundaryNodeInfo()"
            << "\nMore than 2*DIM outside neighbors were found!"
            << std::endl);
      }
   }
}

/*
 *************************************************************************
 *
 *  Set nearest neighbor node
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setNearestNeighborNode(
   pdat::NodeIndex& index)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_index, index);

   /*
    * We assume the number of nearest neighbor nodes is DIM, but allow
    * for there to be more.  Reset the size of the array accordingly.
    */
   if (d_num_nearest_neighbors < index.getDim().getValue()) {
      d_nearest_neighbors[d_num_nearest_neighbors] = index;
      d_num_nearest_neighbors++;
   } else {
      TBOX_ERROR("BoundaryNode: There have already been "
         << d_num_nearest_neighbors
         << "\nregistered with boundary node: "
         << d_index << std::endl);
   }
}

/*
 *************************************************************************
 *
 *  Compute the set of nearest neighbor nodes, given the patch and
 *  flag array index.
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setNearestNeighborNodes(
   tbox::Pointer<pdat::NodeData<int> >& node_flag,
   hier::Index& cut_cell_index)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(d_index, *node_flag, cut_cell_index);

   const tbox::Dimension dim(d_index.getDim());

   int n;

   if (d_num_outside_neighbors < 0) {
      setNumOutsideNeighborNodes(node_flag,
         cut_cell_index);
   }

   /*
    * There should be at least DIM outside neighbors.  If not, we have
    * some sort of convex shaped boundary.
    */
   if (d_num_outside_neighbors < dim.getValue()) {
      TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
         << "\nLess than DIM outside neighbors were found."
         << "\nCannot compute nearest neighbor nodes with fewer"
         << "\nthan DIM outside neighbors." << std::endl);
   }

   /*
    * Form a 2-cell box with the boundary node at the center.  Look for
    * neighboring nodes that are OUTSIDE.
    */
   hier::Index lo(dim);
   hier::Index hi(dim);
   for (int i = 0; i < dim.getValue(); i++) {
      lo(i) = d_index(i) - 1;
      hi(i) = d_index(i);
   }
   hier::Box two_cell_bn_box(lo, hi);

   /*
    * Compute outside neighbors array.
    */
   tbox::Array<pdat::NodeIndex> outside_neighbors(dim.getValue() * 4, pdat::NodeIndex(dim));
   int num_outside_neighbors = 0;
   for (pdat::NodeIterator on(two_cell_bn_box);
        on; on++) {
      pdat::NodeIndex outside_node = on();
      if ((*node_flag)(outside_node) == EmbeddedBoundaryDefines::OUTSIDE) {
         outside_neighbors[num_outside_neighbors] = outside_node;
         num_outside_neighbors++;
      }
      if (num_outside_neighbors >= dim.getValue() * 4) {
         TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
            << "\nMore than 2*DIM outside neighbors were found!"
            << std::endl);
      }
   }

   if (num_outside_neighbors != d_num_outside_neighbors) {
      TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
         << "\nThe number of outside neighbors computed does not"
         << "\ncorrespond with d_num_outside_neighbors.  There"
         << "\nis a bug." << std::endl);
   }

   /*
    * From the array of outside neighbors, determine the DIM
    * outside neighbors that are nearest to bn.  If there are
    * only DIM outside neighbors found, just use those.  Otherwise,
    * find those closest to bn.
    */
   d_num_nearest_neighbors = 0;
   if (d_num_outside_neighbors == dim.getValue()) {
      for (int i = 0; i < dim.getValue(); i++) {
         pdat::NodeIndex neighbor = outside_neighbors[i];
         setNearestNeighborNode(neighbor);
      }
   } else {
      double dist = 0.;
      // find distance == 1 cases
      for (n = 0; n < num_outside_neighbors; n++) {
         dist = 0.;
         hier::Index diff_index(dim, 0);
         for (int i = 0; i < dim.getValue(); i++) {
            diff_index(i) = d_index(i) - outside_neighbors[n](i);
            dist += (double)diff_index(i) * (double)diff_index(i);
         }
         if (tbox::MathUtilities<double>::equalEps(dist, 1.0)) {
            setNearestNeighborNode(outside_neighbors[n]);
         }
         if (d_num_nearest_neighbors == dim.getValue())
            break;
      }

      if (d_num_nearest_neighbors < dim.getValue()) {
         // if they are greater than distance == 1, it doesn't
         // matter which one we pick.
         for (n = 0; n < num_outside_neighbors; n++) {
            dist = 0.;
            hier::Index diff_index(dim, 0);
            for (int i = 0; i < dim.getValue(); i++) {
               diff_index(i) = d_index(i)
                  - outside_neighbors[n](i);
               dist += (double)diff_index(i) * (double)diff_index(i);
            }
            if (dist > 1.0) {
               setNearestNeighborNode(outside_neighbors[n]);
            }
            if (d_num_nearest_neighbors == dim.getValue()) break;
         }
      }

      if (d_num_nearest_neighbors < dim.getValue()) {
         TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
            << "\nDid not find DIM nearest neighbors!"
            << std::endl);
      }
   }
}

/*
 *************************************************************************
 *
 *  Set location of closest boundary point
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setClosestBoundaryPoint(
   const double* location)
{
   for (int i = 0; i < d_index.getDim().getValue(); i++) {
      d_closest_boundary_point[i] = location[i];
   }
}

void
BoundaryNode::setClosestBoundaryPoint(
   const double location,
   const int i)
{
   TBOX_ASSERT(i < d_index.getDim().getValue());

   d_closest_boundary_point[i] = location;
}

/*
 *************************************************************************
 *                                                                       *
 * Set distance and normal to boundary.
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::setDistanceToBoundary(
   tbox::Pointer<hier::Patch>& patch)
{
   TBOX_ASSERT(!patch.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_index, *patch);

   if (d_distance_to_boundary < 0.) {
      const tbox::Pointer<geom::CartesianPatchGeometry> pgeom =
         patch->getPatchGeometry();
      const double* dx = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const hier::Index ifirst = patch->getBox().lower();

      double node_loc[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      double offset;
      double dist = 0.;
      double distsq = 0.;

      for (int i = 0; i < d_index.getDim().getValue(); i++) {
         offset = (double)(d_index(i) - ifirst(i));
         node_loc[i] = xlo[i] + offset * dx[i];
         dist = node_loc[i] - d_closest_boundary_point[i];
         distsq += dist * dist;
      }

      d_distance_to_boundary = sqrt(distsq);

      if (d_distance_to_boundary < s_on_boundary_threshold) {
         d_on_boundary = true;
      }
   }
}

void
BoundaryNode::setDistanceToBoundary(
   const double dist)
{
   d_distance_to_boundary = dist;
}

void
BoundaryNode::setNormalToBoundary(
   tbox::Pointer<hier::Patch>& patch)
{
   TBOX_ASSERT(!patch.isNull());
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_index, *patch);

   if (d_distance_to_boundary < 0.) {
      setDistanceToBoundary(patch);
   }

   const tbox::Pointer<geom::CartesianPatchGeometry> pgeom =
      patch->getPatchGeometry();
   const double* dx = pgeom->getDx();
   const double* xlo = pgeom->getXLower();
   const hier::Index ifirst = patch->getBox().lower();

   double node_loc[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double offset;
   double xdist;

   for (int i = 0; i < d_index.getDim().getValue(); i++) {
      offset = (double)(d_index(i) - ifirst(i));
      node_loc[i] = xlo[i] + offset * dx[i];
      xdist = d_closest_boundary_point[i] - node_loc[i];
      d_normal_to_boundary[i] = xdist / d_distance_to_boundary;
   }
}

void
BoundaryNode::setNormalToBoundary(
   const double* normal)
{
   for (int i = 0; i < d_index.getDim().getValue(); i++) {
      d_normal_to_boundary[i] = normal[i];
   }
}

void
BoundaryNode::setNormalToBoundary(
   const double normal,
   const int i)
{
   d_normal_to_boundary[i] = normal;
}

/*
 *************************************************************************
 *                                                                       *
 * Initialize data in a new boundary node
 *                                                                       *
 *************************************************************************
 */
void
BoundaryNode::initializeBoundaryNodeData()
{
   const tbox::Dimension& dim(d_index.getDim());

   d_num_nearest_neighbors = -1;
   d_num_outside_neighbors = -1;
   d_nearest_neighbors.resizeArray(dim.getValue(), pdat::NodeIndex(dim));
   for (int i = 0; i < dim.getValue(); i++) {
      d_closest_boundary_point[i] = BOUNDARYNODE_LOC_UNDEFINED;
      d_normal_to_boundary[i] = BOUNDARYNODE_LOC_UNDEFINED;
   }
   d_distance_to_boundary = BOUNDARYNODE_LOC_UNDEFINED;
   d_on_boundary = false;
}

/*
 *************************************************************************
 *                                                                       *
 * Copy data from supplied boundary node                                 *
 *                                                                       *
 *************************************************************************
 */

void
BoundaryNode::copyBoundaryNodeData(
   const appu::BoundaryNode& bdry_node)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_index, bdry_node.d_index);

   const tbox::Dimension dim(d_index.getDim());

   d_on_boundary = bdry_node.d_on_boundary;
   d_num_nearest_neighbors = bdry_node.d_num_nearest_neighbors;
   d_num_outside_neighbors = bdry_node.d_num_outside_neighbors;
   d_nearest_neighbors.resizeArray(d_num_nearest_neighbors, pdat::NodeIndex(dim));

   int i;
   for (i = 0; i < d_num_nearest_neighbors; i++) {
      d_nearest_neighbors[i] = bdry_node.d_nearest_neighbors[i];
   }
   for (i = 0; i < dim.getValue(); i++) {
      d_closest_boundary_point[i] = bdry_node.d_closest_boundary_point[i];
      d_normal_to_boundary[i] = bdry_node.d_normal_to_boundary[i];
   }
   d_distance_to_boundary = bdry_node.d_distance_to_boundary;

}

}
}
#endif
