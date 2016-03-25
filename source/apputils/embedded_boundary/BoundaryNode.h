//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/embedded_boundary/BoundaryNode.h $
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Boundary node class for embedded boundary implementations
//

#ifndef included_appu_BoundaryNode
#define included_appu_BoundaryNode

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_pdat_NodeIndex
#include "NodeIndex.h"
#endif
#ifndef included_pdat_NodeData
#include "NodeData.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif

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
 *     - NEAREST_NBR[DIM]      - indices (i,j,k) of the nearest neighbor
 *                             nodes that are OUTSIDE the boundary.
 *
 * @see appu::CutCell  
 */


template<int DIM> class BoundaryNode
{
public:
   /*!
    * Set threshold for determining whether a node is on the boundary
    * (if not set, the default is 1.e-6).
    */
   static void setOnBoundaryThreshold(const double th);
   
   /*!
    * Create a new ``empty'' BoundaryNode.
    */
   BoundaryNode();

   /*!
    * Create a new cut cell with specified node index.
    */
   BoundaryNode(const pdat::NodeIndex<DIM>& in);

   /*!
    * The copy constructor copies the data of the argument cell.
    */
   BoundaryNode(const appu::BoundaryNode<DIM>& bdry_node);
 
   /*!
    * The assignment operator copies the data of the argument cell.
    */
   BoundaryNode& operator=(const appu::BoundaryNode<DIM>& bdry_node);

   /*!
    * The destructor for BoundaryNode.
    */
   ~BoundaryNode();

   /*!
    * Returns the index (i,j,k) of the node.
    */
   pdat::NodeIndex<DIM> getIndex() const;

   /*!
    * Returns whether the boundary node is on the embedded boundary.
    */
   bool getNodeOnBoundary() const;

   /*!
    * Return the number of nearest neighbor nodes.
    */
   int getNumberNearestNeighborNodes() const;

   /*!
    * Return the number of outside neighbor nodes.
    */
   int getNumberOutsideNeighborNodes() const;

   /*!
    * Returns the array of nearest neighbor nodes.
    */
   tbox::Array<pdat::NodeIndex<DIM> > getNearestNeighborNodes() const;

   /*!
    * Returns the designated neighbor node.
    */
   pdat::NodeIndex<DIM> getNearestNeighborNode(const int i) const;

   /*!
    * Returns the location of the closest point on the boundary to
    * the node.
    */
   const double* getClosestBoundaryPoint() const;

   /*!
    * Returns the ith element of the location of the closest point on 
    * the boundary to the node.
    */
   double getClosestBoundaryPoint(const int i) const;

   /*!
    * Returns the distance to the embedded boundary. 
    */
   double getDistanceToBoundary() const;

   /*!
    * Returns the normal vector to the boundary. 
    */
   const double* getNormalToBoundary() const;

   /*!
    * Returns the ith component of the normal vector to the boundary.
    */
   double getNormalToBoundary(
      const int i) const;

   /*!
    * Returns whether the boundary node is on the embedded boundary.
    */
   void setNodeOnBoundary();

   /*!
    * Set the number of outside neighbor nodes for the boundary node.
    */
   void setNumOutsideNeighborNodes(
      tbox::Pointer< pdat::NodeData<DIM,int> >& node_flag,
      hier::Index<DIM>& cut_cell_index);

   /*!
    * Sets the nearest neighbor node.
    */
   void setNearestNeighborNode(pdat::NodeIndex<DIM>& index);

   /*!
    * Set the nearest neighbor nodes for the boundary node.
    */
   void setNearestNeighborNodes(
      tbox::Pointer< pdat::NodeData<DIM,int> >& node_flag,
      hier::Index<DIM>& cut_cell_index);

   /*!
    * Sets the location of the closest point on the b oundary to the node.
    */
   void setClosestBoundaryPoint(const double* location);

   /*!
    * Sets the ith element of the location of the closest point on the 
    * boundary to the node.
    */
   void setClosestBoundaryPoint(const double location,
                                const int i);

   /*!
    * Sets the distance to the embedded boundary. If the patch
    * is provided as an argument, and the distance will be computed.
    * Otherwise, it will be set to the supplied value.
    */
   void setDistanceToBoundary(
      tbox::Pointer<hier::Patch<DIM> > &patch);

   void setDistanceToBoundary(const double dist);

   /*!
    * Sets the normal vector to the embedded boundary. If the patch
    * is provided as an argument, and the normal will be computed.
    * Otherwise, it will be set to the supplied value. 
    */
   void setNormalToBoundary(
      tbox::Pointer<hier::Patch<DIM> > &patch);

   void setNormalToBoundary(const double* normal);

   void setNormalToBoundary(const double normal,
                            const int i);

private:

   /*
    * Initialize data in a new boundary node.
    */
   void initializeBoundaryNodeData();

   /*
    * Copy data from supplied boundary node.
    */
   void copyBoundaryNodeData(const appu::BoundaryNode<DIM>& bdry_node);

   /*
    * Threshold used to determine whether a node is on the boundary
    */
   static double s_on_boundary_threshold;
   
   /*
    * Index of BoundaryNode
    */
   pdat::NodeIndex<DIM> d_index;

   /*
    * Number of nearest and outside neighbor nodes
    */
   int d_num_nearest_neighbors;
   int d_num_outside_neighbors;

   /*
    * Array of nearest neighbor nodes.
    */
   tbox::Array<pdat::NodeIndex<DIM> > d_nearest_neighbors;

   /*
    * Location of closest boundary point, along with distance
    * and normal.
    */
   double d_closest_boundary_point[DIM];
   double d_distance_to_boundary;
   double d_normal_to_boundary[DIM];
   
   /*
    * Boolean specifying whether the boundary node is actually ON the
    * embedded_boundary.
    */
   bool   d_on_boundary;
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoundaryNode.C"
#endif

