/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Boundary node class for embedded boundary implementations 
 *
 ************************************************************************/

#ifndef included_appu_BoundaryNode
#define included_appu_BoundaryNode

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/IOStream.h"

namespace SAMRAI {
namespace appu {

/*!
 * @brief The BoundaryNode struct holds data and methods to define a boundary
 * node (i.e. the first node inside the boundary) on an irregular boundary.
 * An array of boundary nodes is maintained by each "CutCell" object,
 * if the appropriate functions are called to enable boundary node storage.
 * For more information, see the CutCell class documentation.
 *
 * Information maintained by the struct includes the following:
 *
 *     - INDEX                 - node index (i,j,k) of the boundary node
 *     - NEAREST_NBR[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]      - indices (i,j,k) of the nearest neighbor
 *                             nodes that are OUTSIDE the boundary.
 *
 * @see appu::CutCell
 */

class BoundaryNode
{
public:
   /*!
    * Set threshold for determining whether a node is on the boundary
    * (if not set, the default is 1.e-6).
    */
   static void
   setOnBoundaryThreshold(
      const double th);

   /*!
    * Create a new ``empty'' BoundaryNode.
    */
   BoundaryNode(
      const tbox::Dimension& dim);

   /*!
    * Create a new cut cell with specified node index.
    */
   BoundaryNode(
      const pdat::NodeIndex& in);

   /*!
    * The copy constructor copies the data of the argument cell.
    */
   BoundaryNode(
      const appu::BoundaryNode& bdry_node);

   /*!
    * The assignment operator copies the data of the argument cell.
    */
   BoundaryNode&
   operator = (
      const appu::BoundaryNode& bdry_node);

   /*!
    * The destructor for BoundaryNode.
    */
   ~BoundaryNode();

   /*!
    * Returns the index (i,j,k) of the node.
    */
   pdat::NodeIndex
   getIndex() const;

   /*!
    * Returns whether the boundary node is on the embedded boundary.
    */
   bool
   getNodeOnBoundary() const;

   /*!
    * Return the number of nearest neighbor nodes.
    */
   int
   getNumberOfNearestNeighborNodes() const;

   /*!
    * Return the number of outside neighbor nodes.
    */
   int
   getNumberOfOutsideNeighborNodes() const;

   /*!
    * Returns the array of nearest neighbor nodes.
    */
   tbox::Array<pdat::NodeIndex>
   getNearestNeighborNodes() const;

   /*!
    * Returns the designated neighbor node.
    */
   pdat::NodeIndex
   getNearestNeighborNode(
      const int i) const;

   /*!
    * Returns the location of the closest point on the boundary to
    * the node.
    */
   const double *
   getClosestBoundaryPoint() const;

   /*!
    * Returns the ith element of the location of the closest point on
    * the boundary to the node.
    */
   double
   getClosestBoundaryPoint(
      const int i) const;

   /*!
    * Returns the distance to the embedded boundary.
    */
   double
   getDistanceToBoundary() const;

   /*!
    * Returns the normal vector to the boundary.
    */
   const double *
   getNormalToBoundary() const;

   /*!
    * Returns the ith component of the normal vector to the boundary.
    */
   double
   getNormalToBoundary(
      const int i) const;

   /*!
    * Returns whether the boundary node is on the embedded boundary.
    */
   void
   setNodeOnBoundary();

   /*!
    * Set the number of outside neighbor nodes for the boundary node.
    */
   void
   setNumOutsideNeighborNodes(
      tbox::Pointer<pdat::NodeData<int> >& node_flag,
      hier::Index& cut_cell_index);

   /*!
    * Sets the nearest neighbor node.
    */
   void
   setNearestNeighborNode(
      pdat::NodeIndex& index);

   /*!
    * Set the nearest neighbor nodes for the boundary node.
    */
   void
   setNearestNeighborNodes(
      tbox::Pointer<pdat::NodeData<int> >& node_flag,
      hier::Index& cut_cell_index);

   /*!
    * Sets the location of the closest point on the b oundary to the node.
    */
   void
   setClosestBoundaryPoint(
      const double* location);

   /*!
    * Sets the ith element of the location of the closest point on the
    * boundary to the node.
    */
   void
   setClosestBoundaryPoint(
      const double location,
      const int i);

   /*!
    * Sets the distance to the embedded boundary. If the patch
    * is provided as an argument, and the distance will be computed.
    * Otherwise, it will be set to the supplied value.
    */
   void
   setDistanceToBoundary(
      tbox::Pointer<hier::Patch>& patch);

   void
   setDistanceToBoundary(
      const double dist);

   /*!
    * Sets the normal vector to the embedded boundary. If the patch
    * is provided as an argument, and the normal will be computed.
    * Otherwise, it will be set to the supplied value.
    */
   void
   setNormalToBoundary(
      tbox::Pointer<hier::Patch>& patch);

   void
   setNormalToBoundary(
      const double* normal);

   void
   setNormalToBoundary(
      const double normal,
      const int i);

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int BOUNDARYNODE_VERSION;

   /*
    * Static integer constant describing undefined boundary node location.
    */
   static const int BOUNDARYNODE_LOC_UNDEFINED;

   /*
    * Initialize data in a new boundary node.
    */
   void
   initializeBoundaryNodeData();

   /*
    * Copy data from supplied boundary node.
    */
   void
   copyBoundaryNodeData(
      const appu::BoundaryNode& bdry_node);

   /*
    * Threshold used to determine whether a node is on the boundary
    */
   static double s_on_boundary_threshold;

   /*
    * Index of BoundaryNode
    */
   pdat::NodeIndex d_index;

   /*
    * Number of nearest and outside neighbor nodes
    */
   int d_num_nearest_neighbors;
   int d_num_outside_neighbors;

   /*
    * Array of nearest neighbor nodes.
    */
   tbox::Array<pdat::NodeIndex> d_nearest_neighbors;

   /*
    * Location of closest boundary point, along with distance
    * and normal.
    */
   double d_closest_boundary_point[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double d_distance_to_boundary;
   double d_normal_to_boundary[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Boolean specifying whether the boundary node is actually ON the
    * embedded_boundary.
    */
   bool d_on_boundary;
};

}
}
#endif
