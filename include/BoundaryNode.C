//
// File:        BoundaryNode.C
// Package:     SAMRAI application
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 605 $
// Modified:    $Date: 2005-09-09 15:39:55 -0700 (Fri, 09 Sep 2005) $
// Description: Boundary node struct for embedded boundary implementations.
//

#ifndef included_appu_BoundaryNode_C
#define included_appu_BoundaryNode_C

#include "BoundaryNode.h"

using namespace std;
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#define included_assert
#endif

#include "CartesianPatchGeometry.h"
#include "EmbeddedBoundaryDefines.h"
#include "tbox/IEEE.h"
#include "tbox/IOStream.h"
#include "tbox/Utilities.h"


#define BOUNDARYNODE_VERSION 1
#define BOUNDARYNODE_LOC_UNDEFINED -99999


namespace SAMRAI {
   namespace appu {

/*
*************************************************************************
*                                                                       *
* Initialization for static data members.                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> double
BoundaryNode<DIM>::s_on_boundary_threshold = 1.e-6;


/*
*************************************************************************
*                                                                       *
* Static function to set static data members.                           *
*                                                                       *
*************************************************************************
*/
template<int DIM> void
BoundaryNode<DIM>::setOnBoundaryThreshold(const double th)
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

template<int DIM> BoundaryNode<DIM>::BoundaryNode()
{
   hier::Index<DIM> dummy(BOUNDARYNODE_LOC_UNDEFINED);
   hier::IntVector<DIM> loc(0);
   d_index = pdat::NodeIndex<DIM>(dummy,loc);
   initializeBoundaryNodeData();
}


/*
*************************************************************************
*                                                                       *
* Construct a boundary cell, given the cell index.
*                                                                       *
*************************************************************************
*/

template<int DIM> BoundaryNode<DIM>::BoundaryNode(
   const pdat::NodeIndex<DIM>& in) :
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

template<int DIM> BoundaryNode<DIM>::BoundaryNode(
   const appu::BoundaryNode<DIM>& bdry_node) :
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

template<int DIM> BoundaryNode<DIM>& BoundaryNode<DIM>::operator=(
   const appu::BoundaryNode<DIM>& bdry_node)
{
   d_index = bdry_node.d_index;
   copyBoundaryNodeData(bdry_node);
   return(*this);
}


/*
*************************************************************************
*                                                                       *
* Destructor
*                                                                       *
*************************************************************************
*/

template<int DIM> BoundaryNode<DIM>::~BoundaryNode()
{
}

/*
*************************************************************************
*  
*  Return cell index
*                                                                       *
*************************************************************************
*/
template<int DIM> pdat::NodeIndex<DIM> 
BoundaryNode<DIM>::getIndex() const 
{
   return(d_index);
}

/*
*************************************************************************
*  
*  Return whether the node is on the physical boundary of the shape.
*                                                                       *
*************************************************************************
*/
template<int DIM> bool
BoundaryNode<DIM>::getNodeOnBoundary() const 
{
   return(d_on_boundary);
}

/*
*************************************************************************
*  
*  Return number of nearest neighbors
*                                                                       *
*************************************************************************
*/
template<int DIM> int
BoundaryNode<DIM>::getNumberNearestNeighborNodes() const 
{
   return(d_num_nearest_neighbors);
}

/*
*************************************************************************
*  
*  Return number of nearest neighbors
*                                                                       *
*************************************************************************
*/
template<int DIM> int
BoundaryNode<DIM>::getNumberOutsideNeighborNodes() const 
{
   return(d_num_outside_neighbors);
}


/*
*************************************************************************
*  
*  Return array of nearest neighbors
*                                                                       *
*************************************************************************
*/
template<int DIM> tbox::Array<pdat::NodeIndex<DIM> >
BoundaryNode<DIM>::getNearestNeighborNodes() const 
{
   return(d_nearest_neighbors);
}

/*
*************************************************************************
*  
*  Return particular nearest neighbor
*                                                                       *
*************************************************************************
*/
template<int DIM> pdat::NodeIndex<DIM>
BoundaryNode<DIM>::getNearestNeighborNode(const int i) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(i < d_num_nearest_neighbors);
#endif
   return(d_nearest_neighbors[i]);
}

/*
*************************************************************************
*  
*  Return location of closest boundary point
*                                                                       *
*************************************************************************
*/
template<int DIM> const double*
BoundaryNode<DIM>::getClosestBoundaryPoint() const 
{
   return(d_closest_boundary_point);
}

template<int DIM> double
BoundaryNode<DIM>::getClosestBoundaryPoint(const int i) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(i < DIM);
#endif
   return(d_closest_boundary_point[i]);
}

/*
*************************************************************************
*  
*  Return the distance & normal vector to the embedded boundary
*                                                                       *
*************************************************************************
*/

template<int DIM> double
BoundaryNode<DIM>::getDistanceToBoundary() const
{
   if (d_distance_to_boundary < 0.) {
      TBOX_ERROR("BoundaryNode::getDistanceToBoundary()"
                 << "\nYou must first set the distance before accessing it."
                 << "\nCall 'setDistanceToBoundary(patch)'." << endl);      
   }      
   return(d_distance_to_boundary);
}

template<int DIM> const double*
BoundaryNode<DIM>::getNormalToBoundary() const
{
   if (d_distance_to_boundary < 0.) {
      TBOX_ERROR("BoundaryNode::getNormalToBoundary()"
                 << "\nYou must first set the normal before accessing it."
                 << "\nCall 'setNormalToBoundary(patch)'." << endl);
   }
   return(d_normal_to_boundary);
}


template<int DIM> double
BoundaryNode<DIM>::getNormalToBoundary(
   const int i) const
{
   if (d_distance_to_boundary < 0.) {
      TBOX_ERROR("BoundaryNode::getNormalToBoundary()"
                 << "\nYou must first set the normal before accessing it."
                 << "\nEither call 'setNormalToBoundary(patch)' or"
                 << "\nuse the method 'getNormalToBoundary(patch)'."
                 << endl);
   }
   return(d_normal_to_boundary[i]);
}

/*
*************************************************************************
*  
* Set whether the node is on the embedded boundary.
*                                                                       *
*************************************************************************
*/
template<int DIM> void
BoundaryNode<DIM>::setNodeOnBoundary()
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
template<int DIM> void 
BoundaryNode<DIM>::setNumOutsideNeighborNodes(
   tbox::Pointer< pdat::NodeData<DIM,int> >& node_flag,
   hier::Index<DIM>& cut_cell_index)
{
   
   /*
    * Form a 2-cell box with the boundary node at the center.  Look for
    * neighboring nodes that are OUTSIDE.
    */
   hier::Index<DIM> lo;
   hier::Index<DIM> hi;
   for (int i = 0; i < DIM; i++) {
      lo(i) = d_index(i) - 1;
      hi(i) = d_index(i);
   }
   hier::Box<DIM> two_cell_bn_box(lo,hi);

   d_num_outside_neighbors = 0;
   for (typename pdat::NodeIterator<DIM> on(two_cell_bn_box); 
        on; on++) {
      pdat::NodeIndex<DIM> outside_node = on();
      if ((*node_flag)(outside_node) == EmbeddedBoundaryDefines::OUTSIDE) {
         d_num_outside_neighbors++;
      }
      if (d_num_outside_neighbors >= DIM*4) {
         TBOX_ERROR("BoundaryNode::calculateBoundaryNodeInfo()"
                    << "\nMore than 2*DIM outside neighbors were found!"
                    << endl);
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
template<int DIM> void 
BoundaryNode<DIM>::setNearestNeighborNode(
   pdat::NodeIndex<DIM>& index)
{
   /*
    * We assume the number of nearest neighbor nodes is DIM, but allow
    * for there to be more.  Reset the size of the array accordingly.
    */
   if (d_num_nearest_neighbors < DIM) {
      d_nearest_neighbors[d_num_nearest_neighbors] = index;
      d_num_nearest_neighbors++;
   } else {
      TBOX_ERROR("BoundaryNode: There have already been " 
                 << d_num_nearest_neighbors
                 << "\nregistered with boundary node: "
                 << d_index << endl);
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
template<int DIM> void 
BoundaryNode<DIM>::setNearestNeighborNodes(
   tbox::Pointer< pdat::NodeData<DIM,int> >& node_flag,
   hier::Index<DIM>& cut_cell_index)
{

   int i,n;
   
   if (d_num_outside_neighbors < 0) {
      setNumOutsideNeighborNodes(node_flag,
                                 cut_cell_index);
   }
   
   
   /*
    * There should be at least DIM outside neighbors.  If not, we have
    * some sort of convex shaped boundary.
    */
   if (d_num_outside_neighbors < DIM) {
      TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
                 << "\nLess than DIM outside neighbors were found."
                 << "\nCannot compute nearest neighbor nodes with fewer"
                 << "\nthan DIM outside neighbors." << endl);
   }
   

   /*
    * Form a 2-cell box with the boundary node at the center.  Look for
    * neighboring nodes that are OUTSIDE.
    */
   hier::Index<DIM> lo;
   hier::Index<DIM> hi;
   for (int i = 0; i < DIM; i++) {
      lo(i) = d_index(i) - 1;
      hi(i) = d_index(i);
   }
   hier::Box<DIM> two_cell_bn_box(lo,hi);

   /*
    * Compute outside neighbors array.
    */
   tbox::Array<pdat::NodeIndex<DIM> > outside_neighbors(DIM*4);
   int num_outside_neighbors = 0;
   for (typename pdat::NodeIterator<DIM> on(two_cell_bn_box); 
        on; on++) {
      pdat::NodeIndex<DIM> outside_node = on();
      if ((*node_flag)(outside_node) == EmbeddedBoundaryDefines::OUTSIDE) {
         outside_neighbors[num_outside_neighbors] = outside_node;
         num_outside_neighbors++;
      }
      if (num_outside_neighbors >= DIM*4) {
         TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
                    << "\nMore than 2*DIM outside neighbors were found!"
                    << endl);
      }
   }

   if (num_outside_neighbors != d_num_outside_neighbors) {
      TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
                 << "\nThe number of outside neighbors computed does not"
                 << "\ncorrespond with d_num_outside_neighbors.  There"
                 << "\nis a bug." << endl);
   }

   /*
    * From the array of outside neighbors, determine the DIM 
    * outside neighbors that are nearest to bn.  If there are
    * only DIM outside neighbors found, just use those.  Otherwise,
    * find those closest to bn.
    */
   d_num_nearest_neighbors = 0;
   if (d_num_outside_neighbors == DIM) {
      for (i = 0; i < DIM; i++) {
         pdat::NodeIndex<DIM> neighbor = outside_neighbors[i];
         setNearestNeighborNode(neighbor);
      }
   } else {
      double dist = 0.;
      // find distance == 1 cases
      for (n = 0; n < num_outside_neighbors; n++) {
         dist = 0.;
         hier::Index<DIM> diff_index(0);
         for (i = 0; i < DIM; i++) {
            diff_index(i) = d_index(i) - outside_neighbors[n](i);
            dist += (double)diff_index(i) * (double)diff_index(i);
         }
         if (tbox::Utilities::deq(dist,1.0)) {
            setNearestNeighborNode(outside_neighbors[n]);
         }
         if (d_num_nearest_neighbors == DIM) break;
      }

      if (d_num_nearest_neighbors < DIM) {
         // if they are greater than distance == 1, it doesn't 
         // matter which one we pick.
         for (n = 0; n < num_outside_neighbors; n++) {
            dist = 0.;
            hier::Index<DIM> diff_index(0);
            for (i = 0; i < DIM; i++) {
               diff_index(i) = d_index(i) - 
                  outside_neighbors[n](i);
               dist += (double)diff_index(i) * (double)diff_index(i);
            }
            if (dist > 1.0) {
               setNearestNeighborNode(outside_neighbors[n]);
            }
            if (d_num_nearest_neighbors == DIM) break;
         }
      }
               
      if (d_num_nearest_neighbors < DIM) {
         TBOX_ERROR("BoundaryNode::setNearestNeighborNodes()"
                    << "\nDid not find DIM nearest neighbors!"
                    << endl);
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
template<int DIM> void
BoundaryNode<DIM>::setClosestBoundaryPoint(const double* location)
{
   for (int i = 0; i < DIM; i++) {
      d_closest_boundary_point[i] = location[i];
   }
}

template<int DIM> void
BoundaryNode<DIM>::setClosestBoundaryPoint(const double location,
                                           const int i) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(i < DIM);
#endif
   d_closest_boundary_point[i] = location;
}

/*
*************************************************************************
*                                                                       *
* Set distance and normal to boundary.
*                                                                       *
*************************************************************************
*/
template<int DIM> void
BoundaryNode<DIM>::setDistanceToBoundary(
   tbox::Pointer<hier::Patch<DIM> > &patch)
{
   if (d_distance_to_boundary < 0.) {
      const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
         patch->getPatchGeometry();
      const double* dx  = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const hier::Index<DIM> ifirst =patch->getBox().lower();
   
      double node_loc[DIM];
      double offset;
      double dist = 0.;
      double distsq = 0.;
   
      for (int i = 0; i < DIM; i++) {
         offset = (double)(d_index(i)-ifirst(i));
         node_loc[i] = xlo[i] + offset*dx[i];
         dist = node_loc[i] - d_closest_boundary_point[i];
         distsq += dist*dist;
      }

      d_distance_to_boundary = sqrt(distsq);

      if (d_distance_to_boundary < s_on_boundary_threshold) {
         d_on_boundary = true;
      }
   }
}

template<int DIM> void
BoundaryNode<DIM>::setDistanceToBoundary(const double dist)
{
   d_distance_to_boundary = dist;
}

template<int DIM> void
BoundaryNode<DIM>::setNormalToBoundary(
   tbox::Pointer<hier::Patch<DIM> > &patch)
{
   if (d_distance_to_boundary < 0.) {
      setDistanceToBoundary(patch);
   }

   const tbox::Pointer<geom::CartesianPatchGeometry<DIM> > pgeom = 
      patch->getPatchGeometry();
   const double* dx  = pgeom->getDx();
   const double* xlo = pgeom->getXLower();
   const hier::Index<DIM> ifirst =patch->getBox().lower();
   
   double node_loc[DIM];
   double offset;
   double xdist;

   for (int i = 0; i < DIM; i++) {
      offset = (double)(d_index(i)-ifirst(i));
      node_loc[i] = xlo[i] + offset*dx[i];
      xdist = d_closest_boundary_point[i] - node_loc[i];
      d_normal_to_boundary[i] = xdist/d_distance_to_boundary;
   }
}

template<int DIM> void
BoundaryNode<DIM>::setNormalToBoundary(const double* normal)
{
   for (int i = 0; i < DIM; i++) {
      d_normal_to_boundary[i] = normal[i];
   }
}

template<int DIM> void
BoundaryNode<DIM>::setNormalToBoundary(const double normal,
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
template<int DIM> void
BoundaryNode<DIM>::initializeBoundaryNodeData()
{
   int i;

   d_num_nearest_neighbors = -1;
   d_num_outside_neighbors = -1;
   d_nearest_neighbors.resizeArray(DIM);
   for (i = 0; i < DIM; i++) {
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

template<int DIM> void
BoundaryNode<DIM>::copyBoundaryNodeData(
   const appu::BoundaryNode<DIM>& bdry_node)
{
   d_on_boundary = bdry_node.d_on_boundary;
   d_num_nearest_neighbors = bdry_node.d_num_nearest_neighbors;
   d_num_outside_neighbors = bdry_node.d_num_outside_neighbors;
   d_nearest_neighbors.resizeArray(d_num_nearest_neighbors);

   int i;
   for (i = 0; i < d_num_nearest_neighbors; i++) {
      d_nearest_neighbors[i] = bdry_node.d_nearest_neighbors[i];
   }
   for (i = 0; i < DIM; i++) {
      d_closest_boundary_point[i] = bdry_node.d_closest_boundary_point[i];
      d_normal_to_boundary[i] = bdry_node.d_normal_to_boundary[i];
   }
   d_distance_to_boundary = bdry_node.d_distance_to_boundary;

}

}
}
#endif
