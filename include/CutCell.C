//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/apputils/embedded_boundary/CutCell.C $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1846 $
// Modified:    $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description: Cut cell struct for embedded boundary implementations.
//

#ifndef included_appu_CutCell_C
#define included_appu_CutCell_C

#include "CutCell.h"

#include "CellIndex.h"
#include "tbox/IOStream.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#define CUTCELL_VERSION 1

namespace SAMRAI {
    namespace appu {

/*
*************************************************************************
*                                                                       *
* Initialization for static data members.                               *
*                                                                       *
*************************************************************************
*/
template<int DIM> 
bool CutCell<DIM>::s_enable_boundary_node_storage = false;

/*
*************************************************************************
*                                                                       *
* Static function to set whether or not to compute and store boundary   *
* node info.                                                            *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
void CutCell<DIM>::enableBoundaryNodeStorage()
{
   s_enable_boundary_node_storage = true;
}

/*
*************************************************************************
*                                                                       *
* Static function that returns whether boundary node data is enabled.   *
*                                                                       *
*************************************************************************
*/
template<int DIM> 
bool CutCell<DIM>::boundaryNodesEnabled()
{
   return(s_enable_boundary_node_storage);
}

/*
*************************************************************************
*                                                                       *
* Default Constructor                                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> CutCell<DIM>::CutCell()
{
   hier::Index<DIM> dummy(-999);
   d_index = pdat::CellIndex<DIM>(dummy);
   initializeCutCellData();
}


/*
*************************************************************************
*                                                                       *
* Construct a boundary cell, given the cell index.
*                                                                       *
*************************************************************************
*/

template<int DIM> CutCell<DIM>::CutCell(
   const pdat::CellIndex<DIM>& cut_cell) :
   d_index(cut_cell)
{
   initializeCutCellData();
}


/*
*************************************************************************
*                                                                       *
* Copy Constructor                                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> CutCell<DIM>::CutCell(
   const appu::CutCell<DIM>& cut_cell) :
   d_index(cut_cell.d_index)
{
   copyCutCellData(cut_cell);
}


/*
*************************************************************************
*                                                                       *
* Assignment operator                                                   *
*                                                                       *
*************************************************************************
*/
#if 0
template<int DIM> CutCell<DIM>& CutCell<DIM>::operator=(
   const appu::CutCell<DIM>& cut_cell)
{
   d_index = cut_cell.d_index;
   copyCutCellData(cut_cell);
   return(*this);
}
#endif

template<int DIM> CutCell<DIM>& CutCell<DIM>::copy(
   const appu::CutCell<DIM>& cut_cell)
{
   d_index = cut_cell.d_index;
   copyCutCellData(cut_cell);
   return(*this);
}



/*
*************************************************************************
*                                                                       *
* Destructor
*                                                                       *
*************************************************************************
*/

template<int DIM> CutCell<DIM>::~CutCell()
{
}

/*
*************************************************************************
*  
*  Return cell index
*                                                                       *
*************************************************************************
*/
template<int DIM> pdat::CellIndex<DIM> 
CutCell<DIM>::getIndex() const 
{
   return(d_index);
}

/*
*************************************************************************
*  
*  Return volume fraction
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getVolume() const 
{
   return(d_vol_fraction);
}

/*
*************************************************************************
*  
*  Return pointer to cell area fraction vector (dimension 2*DIM).
*                                                                       *
*************************************************************************
*/
template<int DIM> const double* 
CutCell<DIM>::getArea() const 
{
   return(d_area_fraction);
}

/*
*************************************************************************
*  
*  Return cell area fraction for face i.
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getArea(const int i) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < 2*DIM);
#endif
   return(d_area_fraction[i]);
}

/*
*************************************************************************
*  
*  Return pointer to normal vector (dimension DIM)
*                                                                       *
*************************************************************************
*/
template<int DIM> const double* 
CutCell<DIM>::getNormal() const
{
   return(d_normal);
}

/*
*************************************************************************
*  
*  Return normal component for direction i.
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getNormal(const int i) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < DIM);
#endif
   return(d_normal[i]);
}

/*
*************************************************************************
*  
*  Return area of the exposed cut surface.
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getFrontArea() const 
{
   return(d_front_area);
}

/*
*************************************************************************
*  
*  Return pointer to front centroid vector (dimension DIM).
*                                                                       *
*************************************************************************
*/
template<int DIM> const double* 
CutCell<DIM>::getFrontCentroid() const 
{
   return(d_front_centroid);
}

/*
*************************************************************************
*  
*  Return front centroid component for direction i.
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getFrontCentroid(const int i) const 
{
   return(d_front_centroid[i]);
}




/*
*************************************************************************
*  
*  Return volume of cells surrounding the cut-cell
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getSurrVolume() const 
{
   return(d_surr_vol);
}


/*
*************************************************************************
*  
*  Return the new base with components (i,j)
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getNewBase(const int i, const int j) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < DIM);
   TBOX_ASSERT(j < DIM);
#endif
   return(d_newbase[i][j]);
}

/*
*************************************************************************
*  
* Get the number of boundary nodes.
*                                                                       *
*************************************************************************
*/
template<int DIM> int
CutCell<DIM>::getNumberBoundaryNodes() const 
{
   return(d_num_boundary_nodes);
}

/*
*************************************************************************
*  
* Get the boundary node at location i
*                                                                       *
*************************************************************************
*/
template<int DIM> BoundaryNode<DIM>
CutCell<DIM>::getBoundaryNode(const int i) const 
{
   if (!s_enable_boundary_node_storage) {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
                 << "\nBoundary node storage is not enabled.  Use the"
                 << "\nappu::CutCell<DIM>::enableBoundaryNodeStorage()"
                 << "\nmethod to enable boundary node storage." << std::endl);
   } 
   return(d_boundary_nodes[i]);
}


/*
*************************************************************************
*  
*  To be removed...
*                                                                       *
*************************************************************************
*/
template<int DIM> double 
CutCell<DIM>::getFluxFront(const int m) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(m < DIM+2);
#endif
   return(d_flux_front[m]);
}

/*
*************************************************************************
*  
*  Set volume
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setVolume(const double volume) 
{
   d_vol_fraction = volume;
}

/*
*************************************************************************
*  
*  Set cell area fraction for face i (dimension 2*DIM).
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setArea(const double area,
                            const int i) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < 2*DIM);
#endif
   d_area_fraction[i] = area;
}

/*
*************************************************************************
*  
*  Set normal component for direction i (dimension DIM).
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setNormal(const double normal,
                              const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < DIM);
#endif
   d_normal[i] = normal;
}

/*
*************************************************************************
*  
*  Set pointer to normal vector
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setNormal(const double* normal) 
{
   int i;
   for (i =0; i < DIM; i++) {
      d_normal[i] = normal[i];
   }
}

/*
*************************************************************************
*  
*  Set volume of cells surrounding the cut-cell
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setSurrVolume(const double surrvol)
{
   d_surr_vol = surrvol;
}

/*
*************************************************************************
*  
*  Set the frontal area of the cut-cell
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setFrontArea(const double area)
{
   d_front_area = area;
}

/*
*************************************************************************
*  
*  Set the frontal centroid
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setFrontCentroid(const double loc,
                               const int i)
{
   d_front_centroid[i] = loc;
}

/*
*************************************************************************
*  
*  Set whether cut cell is split
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setSplit()
{
   d_is_split = true;
}

/*
*************************************************************************
*  
*  Explicitly set the new base vector with components (i,j)
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setNewBase(const double base,
                         const int i, 
                         const int j)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < DIM);
   TBOX_ASSERT(j < DIM);
#endif
   d_newbase[i][j] = base;
}

/*
*************************************************************************
*  
*  Set new base using supplied vector
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setNewBase(const double* base)
{
   double norm = 0.0;
   int i;
   
   for  (i = 0; i < DIM; i++) {
      norm = norm + base[i]*base[i];
   }
   norm = sqrt(norm);

   /*
    * For the case where the base is all zero, render the normal and new
    * base zero as well.
    */
   if (tbox::MathUtilities<double>::equalEps(norm,0.0)) {
      norm = tbox::MathUtilities<double>::getMax();
   }
      
   for (i = 0; i < DIM; i++) {
      d_normal[i] = base[i]/norm;
      d_newbase[0][i]= d_normal[i];
   }

   if (DIM==2) {
      d_newbase[1][0] = -d_normal[1];
      d_newbase[1][1] = d_normal[0];
   } else if (DIM==3) {
      int min_dir = 0;
      if (tbox::MathUtilities<double>::Abs(d_normal[1]) < 
          tbox::MathUtilities<double>::Abs(d_normal[0])) {
         min_dir = 1;
      }
      if (tbox::MathUtilities<double>::Abs(d_normal[DIM-1]) < 
          tbox::MathUtilities<double>::Abs(d_normal[min_dir])) {
         min_dir = 2;
      }

      norm = 0.0;
      for  (i = 0; i < DIM; i++) {
         d_newbase[1][i] = d_normal[i]*d_normal[min_dir];
      }
      d_newbase[1][min_dir] = d_normal[min_dir]*d_normal[min_dir]-1;

      for  (i = 0; i < DIM; i++) {
         norm = norm + d_newbase[1][i]*d_newbase[1][i];
      }
      norm = sqrt(norm);

      for  (i = 0; i < DIM; i++) {
         d_newbase[1][i] = d_newbase[1][i]/norm;
      }

      d_newbase[DIM-1][0] = d_normal[1]*d_newbase[1][DIM-1] - 
         d_normal[DIM-1]*d_newbase[1][1];
      d_newbase[DIM-1][1] = d_normal[DIM-1]*d_newbase[1][0] - 
         d_normal[0]*d_newbase[1][DIM-1];
      d_newbase[DIM-1][DIM-1] = d_normal[0]*d_newbase[1][1] - 
         d_normal[1]*d_newbase[1][0];
   }   

}

/*
*************************************************************************
*  
*  Set new base using cells normal vector.
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setNewBase()
{
   setNewBase(d_normal);
}

/*
*************************************************************************
*  
* Add the boundary node to the array of boundary nodes maintained
* by the cut cell.
*                                                                       *
*************************************************************************
*/
template<int DIM> void
CutCell<DIM>::setBoundaryNode(const BoundaryNode<DIM>& node) 
{
   if (s_enable_boundary_node_storage) {

      if (d_num_boundary_nodes >= DIM*DIM) {
         TBOX_ERROR("CutCell::setBoundaryNode()"
                    << "\nNumber of boundary nodes set exceeds max allowed"
                    << "(DIM*DIM)." << std::endl);
      } else {
         d_boundary_nodes[d_num_boundary_nodes] = node;
         d_num_boundary_nodes++;
      }
   } else {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
                 << "\nBoundary node storage is not enabled.  Use the"
                 << "\nappu::CutCell<DIM>::enableBoundaryNodeStorage()"
                 << "\nmethod to enable boundary node storage." << std::endl);
   } 
   
}

/*
*************************************************************************
*  
* Add the boundary node to the array of boundary nodes maintained
* by the cut cell.
*                                                                       *
*************************************************************************
*/
template<int DIM> void
CutCell<DIM>::setBoundaryNode(const BoundaryNode<DIM>& node, 
                              const int i) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < DIM*DIM);
#endif
   if (s_enable_boundary_node_storage) {
      d_boundary_nodes[i] = node;
   } else {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
                 << "\nBoundary node storage is not enabled.  Use the"
                 << "\nappu::CutCell<DIM>::enableBoundaryNodeStorage()"
                 << "\nmethod to enable boundary node storage." << std::endl);
   } 
}

/*
*************************************************************************
*  
* Add the boundary node to the particular index i of the array of 
* boundary nodes maintained by the cut cell.
*                                                                       *
*************************************************************************
*/
template<int DIM> void
CutCell<DIM>::setBoundaryNode(const pdat::NodeIndex<DIM>& node) 
{
   if (s_enable_boundary_node_storage) {
      appu::BoundaryNode<DIM> bn(node);
      setBoundaryNode(bn);
   } else {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
                 << "\nBoundary node storage is not enabled.  Use the"
                 << "\nappu::CutCell<DIM>::enableBoundaryNodeStorage()"
                 << "\nmethod to enable boundary node storage." << std::endl);
   } 
}


/*
*************************************************************************
*  
*  To be removed (eventually)
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::setFluxFront(const int m, const double front)
{
   d_flux_front[m] = front;
}

/*
*************************************************************************
*
* Print volume and area data.
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::printVolumeAndAreas(std::ostream &os) const
{
   tbox::pout << std::flush;
   int i;
   
   os << "index = "<< d_index << std::endl;
   os << "volume fraction = "<< d_vol_fraction << std::endl;
   os << "area fraction = ";
   for (i = 0; i < 2*DIM; i++) {
      os << d_area_fraction[i] << " ";
   }
   os << "surr vol = "<< d_surr_vol << std::endl;
   os << std::endl;
   os << std::endl;
}

/*
*************************************************************************
*
* Print normal components.
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::printNormal(std::ostream &os) const
{
   tbox::pout << std::flush;
   os << "index = "<< d_index << std::endl;
   os << "normal = ";
   for (int i = 0; i < DIM; i++) {
      os << d_normal[i] << " ";
   }
   os << "front area = "<< d_front_area << std::endl;
   os << std::endl;
}

/*
*************************************************************************
*
* Print boundary node information
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::printBoundaryNodes(std::ostream &os) const
{
   if (s_enable_boundary_node_storage) {
      tbox::pout << std::flush;
      int i,j;
      
      os << "cut cell index = " << d_index << std::endl;
      os << "   front centroid:  ";
      for (i = 0; i < DIM; i++) {
         os << d_front_centroid[i] << "\t";
      }
      os << std::endl;
      
      os << "   number boundary nodes: " << d_num_boundary_nodes << std::endl;
      for (i = 0; i < d_num_boundary_nodes; i++) {
         BoundaryNode<DIM> bn = getBoundaryNode(i);
         pdat::NodeIndex<DIM> bnode = bn.getIndex();
         os << "   boundary node: " << i << "\t" << bnode << std::endl;
         os << "      closest boundary point loc: " << "\t";
         for (i = 0; i < DIM; i++) {
            os << bn.getClosestBoundaryPoint(i) << "\t";
         }
         os << std::endl;
         os << "      normal to boundary: " << "\t";
         for (i = 0; i < DIM; i++) {
            os << bn.getNormalToBoundary(i) << "\t";
         }
         os << std::endl;
         double dist_to_boundary = bn.getDistanceToBoundary();
         os << "      distance to boundary: " << "\t" 
            << dist_to_boundary << std::endl;
         bool on_boundary = bn.getNodeOnBoundary();
         os << "      on boundary?: " << "\t" << on_boundary << std::endl;
         
         int nn = bn.getNumberNearestNeighborNodes();
         os << "      number nearest neighbors: " << nn << std::endl;
         for (j = 0; j < nn; j++) {
            pdat::NodeIndex<DIM> bnode_nbr = bn.getNearestNeighborNode(j);
            os << "      nearest neighbor: " << j 
               << "\t" << bnode_nbr << std::endl;
         }
      }
   } else {
      os << "Boundary Node data NOT COMPUTED" << std::endl;
   }
   os << std::endl;
}


/*
*************************************************************************
*
* Print ALL data.
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::printAll(std::ostream &os) const
{
   os << std::flush;
   int i,j;

   os << "ptr_boundary_cell = "<< (CutCell<DIM>*) this << std::endl;
   printVolumeAndAreas(os);
   printNormal(os);
   
   for (i = 0; i < DIM; i++) {
     os << "newbase[" << i << "] = ";
     for (j=0; j < DIM; j++) {
	os << d_newbase[i][j] << " ";
     }
     os << std::endl;
   }

   printBoundaryNodes(os);
   os << std::endl;
}

/*
*************************************************************************
*                                                                       *
* The copySourceItem() method is used to copy CutCell data in the  *
* SAMRAI communication infrastructure. This method is required in order * 
* for CutCell to be a templated data type for IndexData<DIM>      * 
* i.e. IndexData<DIM,CutCell<DIM>>.                                   * 
*                                                                       *
*************************************************************************
*/
template<int DIM> void 
CutCell<DIM>::copySourceItem(hier::Index<DIM>& index,
                        const hier::IntVector<DIM>& src_offset,
                        appu::CutCell<DIM>& src_item)
{
   NULL_USE(src_offset);

   /*
    * Copy src_item data into *this.  Note that we don't do
    * anything with the src_offset.  This is because we have
    * access to the index already.
    */
   d_index           = (pdat::CellIndex<DIM>) index;
   copyCutCellData(src_item);
   
}

/*
*************************************************************************
*                                                                       *
* The getDataStreamSize(), packStream(), and unpackStream() methods     *
* are required to template CutCell as IndexData<DIM> type - i.e.        *
* IndexData<DIM,CutCell<DIM>>.  They are used to communicate          *
* CutCell data.                                              *
*                                                                       * 
* The getDataStreamSize() method specifies how many bytes of data       *
* will be packed in the packStream() method.                            *
*                                                                       *
*************************************************************************
*/

template<int DIM> size_t 
CutCell<DIM>::getDataStreamSize()
{
   /*
    * #bytes = 
    *   d_index           (int[DIM]) +
    *   d_boundary_nodes  (int[1 + DIM*DIM*(DIM+1+DIM*DIM)] + 
    *                      double[DIM*DIM*(1+2*DIM)])
    *   d_vol_fraction    (double) +
    *   d_surr_vol        (double) +
    *   d_front_area      (double) +
    *   d_normal          (double[DIM]) +
    *   d_front_centroid  (double[DIM]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM]) 
    */

   size_t bytes =
      (DIM + 1 + DIM*DIM*(DIM+1+DIM*DIM)) * 
      tbox::AbstractStream::sizeofInt() +
      (3 + 4*DIM + DIM+2 + DIM*DIM + DIM*DIM) * 
      tbox::AbstractStream::sizeofDouble();

   return(bytes);
}

/*
*************************************************************************
*                                                                       *
*  Pack message stream.                                                 *
*                                                                       *
*************************************************************************
*/  
template<int DIM> void 
CutCell<DIM>::packStream(tbox::AbstractStream& stream)
{
   int i,j,k;

   /*
    * #ints = 
    *   d_index                (int[DIM]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)]) 
    */
   int int_buff_size = DIM + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += DIM*DIM*(DIM+1+DIM*DIM);
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;
   for (i = 0; i < DIM; i++) {
      ibuffer[counter] = d_index(i);
      counter++;
   }

   ibuffer[counter] = d_num_boundary_nodes;
   counter++;

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (j = 0; j < DIM*DIM; j++) {
         pdat::NodeIndex<DIM> index = d_boundary_nodes[j].getIndex();
         for (i = 0; i < DIM; i++) {
            ibuffer[counter] = index(i);
            counter++;
         }
         bool on_boundary = d_boundary_nodes[j].getNodeOnBoundary();
         int iosb = 0;
         if (on_boundary) iosb = 1;
         ibuffer[counter] = iosb;
         counter++;
         
         // the maximum number of nearest neighbor nodes is DIM
         for (k = 0; k < DIM; k++) {
            index = d_boundary_nodes[j].getNearestNeighborNode(k);
            for (i = 0; i < DIM; i++) {
               ibuffer[counter] = index(i);
               counter++;
            }
         }
      }
   }
      
   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(int_buff_size == counter);
#endif
   stream.pack(ibuffer, counter);

   /*
    * #doubles =  
    *   d_vol_fraction    (double) +
    *   d_surr_vol        (double) +
    *   d_front_area      (double) +
    *   d_normal          (double[DIM]) +
    *   d_front_centroid  (double[DIM]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */ 

   int dbl_buff_size = 3 + 4*DIM + DIM+2 + DIM*DIM;
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += DIM*DIM*(1+2*DIM);
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;

   dbuffer[0] = d_vol_fraction;
   dbuffer[1] = d_surr_vol;
   dbuffer[2] = d_front_area;
   counter += 3;

   for (i = 0; i < DIM; i++) {
      dbuffer[counter] = d_normal[i];
      dbuffer[counter+1] = d_front_centroid[i];
      counter += 2;
   }

   for (i = 0; i < 2*DIM; i++) {
      dbuffer[counter] = d_area_fraction[i];
      counter++;
   }

   for (i = 0; i < DIM+2; i++) {
      dbuffer[counter] = d_flux_front[i];
      counter++;
   }

   for (i = 0; i < DIM; i++) {
      for (j=0; j < DIM; j++) {
         dbuffer[counter] = d_newbase[i][j];
         counter++;
      }
   }

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (i = 0; i < DIM*DIM; i++) {
         dbuffer[counter] = d_boundary_nodes[i].getDistanceToBoundary();
         counter++;
         for (j = 0; j < DIM; j++) {
            dbuffer[counter] = d_boundary_nodes[i].getNormalToBoundary(j);
            dbuffer[counter+1] = d_boundary_nodes[i].getClosestBoundaryPoint(j);
            counter += 2;
         }
      }
   }
   
   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dbl_buff_size == counter);
#endif
   stream.pack(dbuffer, counter);
}

/*
*************************************************************************
*                                                                       *
*  Unpack message stream.                                               *
*                                                                       *
*************************************************************************
*/  
template<int DIM> void 
CutCell<DIM>::unpackStream(tbox::AbstractStream& stream,
                      const hier::IntVector<DIM>& offset)
{
   int i,j,k;

   /*
    * #ints = 
    *   d_index                (int[DIM]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)]) 
    */
   int int_buff_size = DIM + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += DIM*DIM*(DIM+1+DIM*DIM);
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;
   stream.unpack(ibuffer, int_buff_size);
   pdat::CellIndex<DIM> cell_index;
   for (i = 0; i < DIM; i++) {
      cell_index(i) = ibuffer[counter];
      counter++;
   }
   d_index = cell_index + offset;

   int d_num_boundary_nodes = ibuffer[counter];
   counter++;

   if (s_enable_boundary_node_storage) {

      for (j = 0; j < d_num_boundary_nodes; j++) {
         pdat::NodeIndex<DIM> node_index;
         for (i = 0; i < DIM; i++) {
            node_index(i) = ibuffer[counter];
            counter++;
         }
         
         BoundaryNode<DIM> bn(node_index);
         
         int iosb = ibuffer[counter];
         counter++;
         if (iosb == 1) bn.setNodeOnBoundary();
         
         for (k = 0; k < DIM; k++) {
            for (i = 0; i < DIM; i++) {
               node_index(i) = ibuffer[counter];
               counter++;
            }
            bn.setNearestNeighborNode(node_index);
         }
         d_boundary_nodes[j] = bn;
      }
   }

   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(int_buff_size == counter);
#endif

   /*
    * #doubles =  
    *   d_vol_fraction    (double) +
    *   d_surr_vol        (double) +
    *   d_front_area      (double) +
    *   d_normal          (double[DIM]) +
    *   d_front_centroid  (double[DIM]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */ 

   int dbl_buff_size = 3 + 4*DIM + DIM+2 + DIM*DIM;
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += DIM*DIM*(1+2*DIM);
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;
   stream.unpack(dbuffer, dbl_buff_size);
   d_vol_fraction = dbuffer[0];
   d_surr_vol     = dbuffer[1];
   d_front_area   = dbuffer[2];
   counter += 3;

   for (i = 0; i < DIM; i++) {
      d_normal[i] = dbuffer[counter];
      d_front_centroid[i] = dbuffer[counter+1];
      counter += 2;
   }

   for (i = 0; i < 2*DIM; i++) {
      d_area_fraction[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < DIM+2; i++) {
      d_flux_front[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < DIM; i++) {
      for (j=0; j < DIM; j++) {
         d_newbase[i][j] = dbuffer[counter];
         counter++;
      }
   } 

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (i = 0; i < DIM*DIM; i++) {
         d_boundary_nodes[i].setDistanceToBoundary(dbuffer[counter]);
         counter++;
         for (j = 0; j < DIM; j++) {
            d_boundary_nodes[i].setNormalToBoundary(dbuffer[counter],j);
            d_boundary_nodes[i].setClosestBoundaryPoint(dbuffer[counter+1],j);
            counter += 2;
         }
      } 
   }
   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dbl_buff_size == counter);
#endif

}

/*
*************************************************************************
*                                                                       *
* The putToDatabase() and getFromDatabase() methods are required to     *
* template CutCell as IndexData<DIM> type                         *
* i.e. IndexData<DIM,CutCell<DIM>>.  They are used to write/read      *
* CutCell data to/from the restart database.                       * 
*                                                                       *
* The following writes data to the restart database.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
CutCell<DIM>::putToDatabase(
   tbox::Pointer<tbox::Database>& database)
{
   int i,j,k;

   /*
    * #ints = 
    *   d_index                (int[DIM]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)]) 
    */
   int int_buff_size = DIM + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += DIM*DIM*(DIM+1+DIM*DIM);
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;
   for (i = 0; i < DIM; i++) {
      ibuffer[counter] = d_index(i);
      counter++;
   }

   ibuffer[counter] = d_num_boundary_nodes;
   counter++;

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (j = 0; j < DIM*DIM; j++) {
         pdat::NodeIndex<DIM> index = d_boundary_nodes[j].getIndex();
         for (i = 0; i < DIM; i++) {
            ibuffer[counter] = index(i);
            counter++;
         }
         bool on_boundary = d_boundary_nodes[j].getNodeOnBoundary();
         int iosb = 0;
         if (on_boundary) iosb = 1;
         ibuffer[counter] = iosb;
         counter++;

         // the maximum number of nearest neighbor nodes is DIM
         for (k = 0; k < DIM; k++) {
            index = d_boundary_nodes[j].getNearestNeighborNode(k);
            for (i = 0; i < DIM; i++) {
               ibuffer[counter] = index(i);
               counter++;
            }
         }
      }
   }

   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(int_buff_size == counter);
#endif
   database->putIntegerArray("ibuffer", ibuffer, int_buff_size);


   /*
    * #doubles =  
    *   d_vol_fraction    (double) +
    *   d_surr_vol        (double) +
    *   d_front_area      (double) +
    *   d_normal          (double[DIM]) +
    *   d_front_centroid  (double[DIM]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */ 

   int dbl_buff_size = 3 + 4*DIM + DIM+2 + DIM*DIM;
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += DIM*DIM*(1+2*DIM);
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;

   dbuffer[0] = d_vol_fraction;
   dbuffer[1] = d_surr_vol;
   dbuffer[2] = d_front_area;
   counter += 3;

   for (i = 0; i < DIM; i++) {
      dbuffer[counter] = d_normal[i];
      dbuffer[counter+1] = d_front_centroid[i];
      counter += 2;
   }

   for (i = 0; i < 2*DIM; i++) {
      dbuffer[counter] = d_area_fraction[i];
      counter++;
   }

   for (i = 0; i < DIM+2; i++) {
      dbuffer[counter] = d_flux_front[i];
      counter++;
   }

   for (i = 0; i < DIM; i++) {
      for (j=0; j < DIM; j++) {
         dbuffer[counter] = d_newbase[i][j];
         counter++;
      }
   }

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (i = 0; i < DIM*DIM; i++) {
         dbuffer[counter] = d_boundary_nodes[i].getDistanceToBoundary();
         counter++;
      }
   }

   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dbl_buff_size == counter);
#endif
   database->putDoubleArray("dbuffer", dbuffer, dbl_buff_size);

}



/*
*************************************************************************
*                                                                       *
*  Read data from restart.                                              *
*                                                                       *
*************************************************************************
*/  

template<int DIM> void 
CutCell<DIM>::getFromDatabase(
   tbox::Pointer<tbox::Database>& database)
{

   int i,j,k;
   /*
    * #ints = 
    *   d_index                (int[DIM]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)]) 
    */
   int int_buff_size = DIM + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += DIM*DIM*(DIM+1+DIM*DIM);
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;

   database->getIntegerArray("ibuffer", ibuffer, int_buff_size);
   pdat::CellIndex<DIM> cell_index;
   for (i = 0; i < DIM; i++) {
      cell_index(i) = ibuffer[counter];
      counter++;
   }
   d_index = cell_index;

   int d_num_boundary_nodes = ibuffer[counter];
   counter++;

   if (s_enable_boundary_node_storage) {
      for (j = 0; j < d_num_boundary_nodes; j++) {
         pdat::NodeIndex<DIM> node_index;
         for (i = 0; i < DIM; i++) {
            node_index(i) = ibuffer[counter];
            counter++;
         }

         BoundaryNode<DIM> bn(node_index);
         
         int iosb = ibuffer[counter];
         counter++;
         if (iosb == 1) bn.setNodeOnBoundary();
         
         for (k = 0; k < DIM; k++) {
            for (i = 0; i < DIM; i++) {
               node_index(i) = ibuffer[counter];
               counter++;
            }
            bn.setNearestNeighborNode(node_index);
         }
         d_boundary_nodes[j] = bn;
      }
   }

   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(int_buff_size == counter);
#endif

   /*
    * #doubles =  
    *   d_vol_fraction    (double) +
    *   d_surr_vol        (double) +
    *   d_front_area      (double) +
    *   d_normal          (double[DIM]) +
    *   d_front_centroid  (double[DIM]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */ 

   int dbl_buff_size = 3 + 4*DIM + DIM+2 + DIM*DIM;
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += DIM*DIM*(1+2*DIM);
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;

   database->getDoubleArray("dbuffer", dbuffer, dbl_buff_size);
   d_vol_fraction = dbuffer[0];
   d_surr_vol     = dbuffer[1];
   d_front_area   = dbuffer[2];
   counter += 3;

   for (i = 0; i < DIM; i++) {
      d_normal[i] = dbuffer[counter];
      d_front_centroid[i] = dbuffer[counter+1];
      counter += 2;
   }

   for (i = 0; i < 2*DIM; i++) {
      d_area_fraction[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < DIM+2; i++) {
      d_flux_front[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < DIM; i++) {
      for (j=0; j < DIM; j++) {
         d_newbase[i][j] = dbuffer[counter];
         counter++;
      }
   } 

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (i = 0; i < DIM*DIM; i++) {
         d_boundary_nodes[i].setDistanceToBoundary(dbuffer[counter]);
         counter++;
      } 
   }
      
   /*
    * The counter should equal the buffer size.  Otherwise, there is an
    * error.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(dbl_buff_size == counter);
#endif

}

/*
*************************************************************************
*                                                                       *
* Initialize data in a new cut cell.
*                                                                       *
*************************************************************************
*/
template<int DIM> void
CutCell<DIM>::initializeCutCellData()
{
   d_vol_fraction = tbox::MathUtilities<double>::getSignalingNaN();
   d_surr_vol = tbox::MathUtilities<double>::getSignalingNaN();
   d_front_area = tbox::MathUtilities<double>::getSignalingNaN();

   int i;
   for (i = 0; i < DIM; i++) {
      d_normal[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_front_centroid[i] = tbox::MathUtilities<double>::getSignalingNaN();
   }

   int j;
   for (i = 0; i < DIM; i++) {
      for (j=0; j < DIM; j++) {
         d_newbase[i][j] = tbox::MathUtilities<double>::getSignalingNaN();
      }
   }

   for (i = 0; i < 2*DIM; i++) {
      d_area_fraction[i] = tbox::MathUtilities<double>::getSignalingNaN();
   }

   for (i = 0; i < DIM+2; i++) {
      d_flux_front[i] = 0.;
   }

   d_num_boundary_nodes = 0;
   if (s_enable_boundary_node_storage) {
      d_boundary_nodes.resizeArray(DIM*DIM);
      for (i = 0; i < DIM*DIM; i++) {
         BoundaryNode<DIM> bn;
         d_boundary_nodes[i] = bn;
      }
      
      d_is_split = false;
   }
   
}

/*
*************************************************************************
*                                                                       *
* Copy data from supplied cut cell                                      *
*                                                                       *
*************************************************************************
*/

template<int DIM> void
CutCell<DIM>::copyCutCellData(
   const appu::CutCell<DIM>& cut_cell)
{
   d_vol_fraction    = cut_cell.d_vol_fraction;
   d_surr_vol        = cut_cell.d_surr_vol;
   d_front_area      = cut_cell.d_front_area;

   int i;
   for (i = 0; i < DIM; i++) {
      d_normal[i]           = cut_cell.d_normal[i];
      d_front_centroid[i]   = cut_cell.d_front_centroid[i];
      for (int j = 0; j < DIM; j++) {
         d_newbase[i][j]    = cut_cell.d_newbase[i][j];
      }
   }
   for (i = 0; i < 2*DIM; i++) {
      d_area_fraction[i]    = cut_cell.d_area_fraction[i];
   }

   for (i = 0; i < DIM+2; i++) {
      d_flux_front[i] = cut_cell.d_flux_front[i];
   }

   d_num_boundary_nodes = cut_cell.d_num_boundary_nodes;
   if (s_enable_boundary_node_storage) {
      d_boundary_nodes.resizeArray(DIM*DIM);
      for (i = 0; i < d_num_boundary_nodes; i++) {
         d_boundary_nodes[i] = cut_cell.d_boundary_nodes[i];
      }
      for (i = d_num_boundary_nodes+1; i < DIM*DIM; i++) {   
         BoundaryNode<DIM> bn;
         d_boundary_nodes[i] = bn;
      }
   }
}

}
}
#endif
