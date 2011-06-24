/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Cut cell struct for embedded boundary implementations. 
 *
 ************************************************************************/

#ifndef included_appu_CutCell_C
#define included_appu_CutCell_C

#include "SAMRAI/appu/CutCell.h"

#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace appu {

const int CutCell::CUTCELL_VERSION = 1;

/*
 *************************************************************************
 *                                                                       *
 * Initialization for static data members.                               *
 *                                                                       *
 *************************************************************************
 */
bool CutCell::s_enable_boundary_node_storage = false;

/*
 *************************************************************************
 *                                                                       *
 * Static function to set whether or not to compute and store boundary   *
 * node info.                                                            *
 *                                                                       *
 *************************************************************************
 */

void CutCell::enableBoundaryNodeStorage()
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
bool CutCell::boundaryNodesEnabled()
{
   return s_enable_boundary_node_storage;
}

/*
 *************************************************************************
 *                                                                       *
 * Constructor                                                           *
 *                                                                       *
 *************************************************************************
 */

CutCell::CutCell():
   d_index(pdat::CellIndex(hier::Index(
      tbox::Dimension::getInvalidDimension())))
{
   initializeCutCellData();
}

/*
 *************************************************************************
 *                                                                       *
 * Constructor                                                           *
 *                                                                       *
 *************************************************************************
 */

CutCell::CutCell(
   const tbox::Dimension& dim):
   d_index(pdat::CellIndex(hier::Index(dim, -999)))
{
   initializeCutCellData();
}

/*
 *************************************************************************
 *                                                                       *
 * Construct a boundary cell, given the cell index.
 *                                                                       *
 *************************************************************************
 */

CutCell::CutCell(
   const pdat::CellIndex& cut_cell):
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

CutCell::CutCell(
   const appu::CutCell& cut_cell):
   d_index(cut_cell.d_index)
{
   copyCutCellData(cut_cell);
}

/*
 *************************************************************************
 *
 * copy
 *
 *************************************************************************
 */

CutCell& CutCell::copy(
   const appu::CutCell& cut_cell)
{
   d_index = cut_cell.d_index;
   copyCutCellData(cut_cell);
   return *this;
}

/*
 *************************************************************************
 *                                                                       *
 * Destructor
 *                                                                       *
 *************************************************************************
 */

CutCell::~CutCell()
{
}

/*
 *************************************************************************
 *
 *  Return cell index
 *                                                                       *
 *************************************************************************
 */
pdat::CellIndex
CutCell::getIndex() const
{
   return d_index;
}

/*
 *************************************************************************
 *
 *  Return volume fraction
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getVolume() const
{
   return d_vol_fraction;
}

/*
 *************************************************************************
 *
 *  Return pointer to cell area fraction vector (dimension 2*DIM).
 *                                                                       *
 *************************************************************************
 */
const double *
CutCell::getArea() const
{
   return d_area_fraction;
}

/*
 *************************************************************************
 *
 *  Return cell area fraction for face i.
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getArea(
   const int i) const
{
   TBOX_ASSERT(i < 2 * d_index.getDim().getValue());

   return d_area_fraction[i];
}

/*
 *************************************************************************
 *
 *  Return pointer to normal vector (dimension DIM)
 *                                                                       *
 *************************************************************************
 */
const double *
CutCell::getNormal() const
{
   return d_normal;
}

/*
 *************************************************************************
 *
 *  Return normal component for direction i.
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getNormal(
   const int i) const
{
   TBOX_ASSERT(i < d_index.getDim().getValue());

   return d_normal[i];
}

/*
 *************************************************************************
 *
 *  Return area of the exposed cut surface.
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getFrontArea() const
{
   return d_front_area;
}

/*
 *************************************************************************
 *
 *  Return pointer to front centroid vector (dimension DIM).
 *                                                                       *
 *************************************************************************
 */
const double *
CutCell::getFrontCentroid() const
{
   return d_front_centroid;
}

/*
 *************************************************************************
 *
 *  Return front centroid component for direction i.
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getFrontCentroid(
   const int i) const
{
   return d_front_centroid[i];
}

/*
 *************************************************************************
 *
 *  Return volume of cells surrounding the cut-cell
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getSurrVolume() const
{
   return d_surr_vol;
}

/*
 *************************************************************************
 *
 *  Return the new base with components (i,j)
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getNewBase(
   const int i,
   const int j) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(i < d_index.getDim().getValue());
   TBOX_ASSERT(j < d_index.getDim().getValue());
#endif
   return d_newbase[i][j];
}

/*
 *************************************************************************
 *
 * Get the number of boundary nodes.
 *                                                                       *
 *************************************************************************
 */
int
CutCell::getNumberOfBoundaryNodes() const
{
   return d_num_boundary_nodes;
}

/*
 *************************************************************************
 *
 * Get the boundary node at location i
 *                                                                       *
 *************************************************************************
 */
BoundaryNode
CutCell::getBoundaryNode(
   const int i) const
{
   if (!s_enable_boundary_node_storage) {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
         << "\nBoundary node storage is not enabled.  Use the"
         << "\nappu::CutCell::enableBoundaryNodeStorage()"
         << "\nmethod to enable boundary node storage." << std::endl);
   }
   return d_boundary_nodes[i];
}

/*
 *************************************************************************
 *
 *  To be removed...
 *                                                                       *
 *************************************************************************
 */
double
CutCell::getFluxFront(
   const int m) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(m < d_index.getDim().getValue() + 2);
#endif
   return d_flux_front[m];
}

/*
 *************************************************************************
 *
 *  Set volume
 *                                                                       *
 *************************************************************************
 */
void
CutCell::setVolume(
   const double volume)
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
void
CutCell::setArea(
   const double area,
   const int i)
{

   TBOX_ASSERT(i < 2 * d_index.getDim().getValue());

   d_area_fraction[i] = area;
}

/*
 *************************************************************************
 *
 *  Set normal component for direction i (dimension DIM).
 *                                                                       *
 *************************************************************************
 */
void
CutCell::setNormal(
   const double normal,
   const int i)
{
   TBOX_ASSERT(i < d_index.getDim().getValue());

   d_normal[i] = normal;
}

/*
 *************************************************************************
 *
 *  Set pointer to normal vector
 *                                                                       *
 *************************************************************************
 */
void
CutCell::setNormal(
   const double* normal)
{
   int i;
   for (i = 0; i < d_index.getDim().getValue(); i++) {
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
void
CutCell::setSurrVolume(
   const double surrvol)
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
void
CutCell::setFrontArea(
   const double area)
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
void
CutCell::setFrontCentroid(
   const double loc,
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
void
CutCell::setSplit()
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
void
CutCell::setNewBase(
   const double base,
   const int i,
   const int j)
{
   TBOX_ASSERT(i < d_index.getDim().getValue());
   TBOX_ASSERT(j < d_index.getDim().getValue());

   d_newbase[i][j] = base;
}

/*
 *************************************************************************
 *
 *  Set new base using supplied vector
 *                                                                       *
 *************************************************************************
 */
void
CutCell::setNewBase(
   const double* base)
{
   const tbox::Dimension& dim(d_index.getDim());

   double norm = 0.0;
   int i;

   for (i = 0; i < dim.getValue(); i++) {
      norm = norm + base[i] * base[i];
   }
   norm = sqrt(norm);

   /*
    * For the case where the base is all zero, render the normal and new
    * base zero as well.
    */
   if (tbox::MathUtilities<double>::equalEps(norm, 0.0)) {
      norm = tbox::MathUtilities<double>::getMax();
   }

   for (i = 0; i < dim.getValue(); i++) {
      d_normal[i] = base[i] / norm;
      d_newbase[0][i] = d_normal[i];
   }

   if ((dim == tbox::Dimension(2))) {
      d_newbase[1][0] = -d_normal[1];
      d_newbase[1][1] = d_normal[0];
   } else if ((dim == tbox::Dimension(3))) {
      int min_dir = 0;
      if (tbox::MathUtilities<double>::Abs(d_normal[1]) <
          tbox::MathUtilities<double>::Abs(d_normal[0])) {
         min_dir = 1;
      }
      if (tbox::MathUtilities<double>::Abs(d_normal[dim.getValue() - 1]) <
          tbox::MathUtilities<double>::Abs(d_normal[min_dir])) {
         min_dir = 2;
      }

      norm = 0.0;
      for (i = 0; i < dim.getValue(); i++) {
         d_newbase[1][i] = d_normal[i] * d_normal[min_dir];
      }
      d_newbase[1][min_dir] = d_normal[min_dir] * d_normal[min_dir] - 1;

      for (i = 0; i < dim.getValue(); i++) {
         norm = norm + d_newbase[1][i] * d_newbase[1][i];
      }
      norm = sqrt(norm);

      for (i = 0; i < dim.getValue(); i++) {
         d_newbase[1][i] = d_newbase[1][i] / norm;
      }

      d_newbase[dim.getValue() - 1][0] = d_normal[1] * d_newbase[1][dim.getValue() - 1]
         - d_normal[dim.getValue() - 1] * d_newbase[1][1];
      d_newbase[dim.getValue() - 1][1] = d_normal[dim.getValue() - 1] * d_newbase[1][0]
         - d_normal[0] * d_newbase[1][dim.getValue() - 1];
      d_newbase[dim.getValue() - 1][dim.getValue() - 1] = d_normal[0] * d_newbase[1][1]
         - d_normal[1] * d_newbase[1][0];
   }

}

/*
 *************************************************************************
 *
 *  Set new base using cells normal vector.
 *                                                                       *
 *************************************************************************
 */
void
CutCell::setNewBase()
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
void
CutCell::setBoundaryNode(
   const BoundaryNode& node)
{
   const tbox::Dimension& dim(d_index.getDim());

   if (s_enable_boundary_node_storage) {

      if (d_num_boundary_nodes >= dim.getValue() * dim.getValue()) {
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
         << "\nappu::CutCell::enableBoundaryNodeStorage()"
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
void
CutCell::setBoundaryNode(
   const BoundaryNode& node,
   const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim(d_index.getDim());
   TBOX_ASSERT(i < dim.getValue() * dim.getValue());
#endif

   if (s_enable_boundary_node_storage) {
      d_boundary_nodes[i] = node;
   } else {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
         << "\nBoundary node storage is not enabled.  Use the"
         << "\nappu::CutCell::enableBoundaryNodeStorage()"
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
void
CutCell::setBoundaryNode(
   const pdat::NodeIndex& node)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_index, node);

   if (s_enable_boundary_node_storage) {
      appu::BoundaryNode bn(node);
      setBoundaryNode(bn);
   } else {
      TBOX_ERROR("appu::CutCell::setBoundaryNode()"
         << "\nBoundary node storage is not enabled.  Use the"
         << "\nappu::CutCell::enableBoundaryNodeStorage()"
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
void
CutCell::setFluxFront(
   const int m,
   const double front)
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
void
CutCell::printVolumeAndAreas(
   std::ostream& os) const
{
   tbox::pout << std::flush;
   int i;

   os << "index = " << d_index << std::endl;
   os << "volume fraction = " << d_vol_fraction << std::endl;
   os << "area fraction = ";
   for (i = 0; i < 2 * d_index.getDim().getValue(); i++) {
      os << d_area_fraction[i] << " ";
   }
   os << "surr vol = " << d_surr_vol << std::endl;
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
void
CutCell::printNormal(
   std::ostream& os) const
{
   tbox::pout << std::flush;
   os << "index = " << d_index << std::endl;
   os << "normal = ";
   for (int i = 0; i < d_index.getDim().getValue(); i++) {
      os << d_normal[i] << " ";
   }
   os << "front area = " << d_front_area << std::endl;
   os << std::endl;
}

/*
 *************************************************************************
 *
 * Print boundary node information
 *                                                                       *
 *************************************************************************
 */
void
CutCell::printBoundaryNodes(
   std::ostream& os) const
{
   if (s_enable_boundary_node_storage) {
      const tbox::Dimension& dim(d_index.getDim());

      tbox::pout << std::flush;
      int i, j;

      os << "cut cell index = " << d_index << std::endl;
      os << "   front centroid:  ";
      for (i = 0; i < dim.getValue(); i++) {
         os << d_front_centroid[i] << "\t";
      }
      os << std::endl;

      os << "   number boundary nodes: " << d_num_boundary_nodes << std::endl;
      for (i = 0; i < d_num_boundary_nodes; i++) {
         BoundaryNode bn = getBoundaryNode(i);
         pdat::NodeIndex bnode = bn.getIndex();
         os << "   boundary node: " << i << "\t" << bnode << std::endl;
         os << "      closest boundary point loc: " << "\t";
         for (i = 0; i < dim.getValue(); i++) {
            os << bn.getClosestBoundaryPoint(i) << "\t";
         }
         os << std::endl;
         os << "      normal to boundary: " << "\t";
         for (i = 0; i < dim.getValue(); i++) {
            os << bn.getNormalToBoundary(i) << "\t";
         }
         os << std::endl;
         double dist_to_boundary = bn.getDistanceToBoundary();
         os << "      distance to boundary: " << "\t"
            << dist_to_boundary << std::endl;
         bool on_boundary = bn.getNodeOnBoundary();
         os << "      on boundary?: " << "\t" << on_boundary << std::endl;

         int nn = bn.getNumberOfNearestNeighborNodes();
         os << "      number nearest neighbors: " << nn << std::endl;
         for (j = 0; j < nn; j++) {
            pdat::NodeIndex bnode_nbr = bn.getNearestNeighborNode(j);
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
void
CutCell::printAll(
   std::ostream& os) const
{
   const tbox::Dimension& dim(d_index.getDim());

   os << std::flush;
   int i, j;

   os << "ptr_boundary_cell = " << (CutCell *)this << std::endl;
   printVolumeAndAreas(os);
   printNormal(os);

   for (i = 0; i < dim.getValue(); i++) {
      os << "newbase[" << i << "] = ";
      for (j = 0; j < dim.getValue(); j++) {
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
 * i.e. IndexData<CutCell>.                                   *
 *                                                                       *
 *************************************************************************
 */
void
CutCell::copySourceItem(
   hier::Index& index,
   const hier::IntVector& src_offset,
   appu::CutCell& src_item)
{
   NULL_USE(src_offset);

   /*
    * Copy src_item data into *this.  Note that we don't do
    * anything with the src_offset.  This is because we have
    * access to the index already.
    */
   d_index = (pdat::CellIndex)index;
   copyCutCellData(src_item);

}

/*
 *************************************************************************
 *                                                                       *
 * The getDataStreamSize(), packStream(), and unpackStream() methods     *
 * are required to template CutCell as IndexData<DIM> type - i.e.        *
 * IndexData<CutCell>.  They are used to communicate          *
 * CutCell data.                                              *
 *                                                                       *
 * The getDataStreamSize() method specifies how many bytes of data       *
 * will be packed in the packStream() method.                            *
 *                                                                       *
 *************************************************************************
 */

size_t
CutCell::getDataStreamSize()
{
   /*
    * #bytes =
    *   d_index           (int[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_boundary_nodes  (int[1 + DIM*DIM*(DIM+1+DIM*DIM)] +
    *                      double[DIM*DIM*(1+2*DIM)])
    *   d_vol_fraction    (double) +
    *   d_surr_vol        (double) +
    *   d_front_area      (double) +
    *   d_normal          (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_front_centroid  (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    */

   int dim_value = d_index.getDim().getValue();

   size_t bytes =
      (dim_value + 1 + dim_value * dim_value * (dim_value + 1 + dim_value * dim_value))
      * tbox::MessageStream::getSizeof<int>()
      + (3 + 4 * dim_value + dim_value + 2 + dim_value * dim_value + dim_value * dim_value)
      * tbox::MessageStream::getSizeof<double>();

   return bytes;
}

/*
 *************************************************************************
 *                                                                       *
 *  Pack message stream.                                                 *
 *                                                                       *
 *************************************************************************
 */
void
CutCell::packStream(
   tbox::MessageStream& stream)
{
   const tbox::Dimension& dim(d_index.getDim());

   int i, j, k;

   /*
    * #ints =
    *   d_index                (int[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)])
    */
   int int_buff_size = dim.getValue() + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += dim.getValue() * dim.getValue()
         * (dim.getValue() + 1 + dim.getValue() * dim.getValue());
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;
   for (i = 0; i < dim.getValue(); i++) {
      ibuffer[counter] = d_index(i);
      counter++;
   }

   ibuffer[counter] = d_num_boundary_nodes;
   counter++;

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (j = 0; j < dim.getValue() * dim.getValue(); j++) {
         pdat::NodeIndex index = d_boundary_nodes[j].getIndex();
         for (i = 0; i < dim.getValue(); i++) {
            ibuffer[counter] = index(i);
            counter++;
         }
         bool on_boundary = d_boundary_nodes[j].getNodeOnBoundary();
         int iosb = 0;
         if (on_boundary) iosb = 1;
         ibuffer[counter] = iosb;
         counter++;

         // the maximum number of nearest neighbor nodes is DIM
         for (k = 0; k < dim.getValue(); k++) {
            index = d_boundary_nodes[j].getNearestNeighborNode(k);
            for (i = 0; i < dim.getValue(); i++) {
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
    *   d_normal          (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_front_centroid  (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */

   int dbl_buff_size = 3 + 4 * dim.getValue() + dim.getValue() + 2 + dim.getValue() * dim.getValue();
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += dim.getValue() * dim.getValue() * (1 + 2 * dim.getValue());
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;

   dbuffer[0] = d_vol_fraction;
   dbuffer[1] = d_surr_vol;
   dbuffer[2] = d_front_area;
   counter += 3;

   for (i = 0; i < dim.getValue(); i++) {
      dbuffer[counter] = d_normal[i];
      dbuffer[counter + 1] = d_front_centroid[i];
      counter += 2;
   }

   for (i = 0; i < 2 * dim.getValue(); i++) {
      dbuffer[counter] = d_area_fraction[i];
      counter++;
   }

   for (i = 0; i < dim.getValue() + 2; i++) {
      dbuffer[counter] = d_flux_front[i];
      counter++;
   }

   for (i = 0; i < dim.getValue(); i++) {
      for (j = 0; j < dim.getValue(); j++) {
         dbuffer[counter] = d_newbase[i][j];
         counter++;
      }
   }

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is DIM*DIM
      for (i = 0; i < dim.getValue() * dim.getValue(); i++) {
         dbuffer[counter] = d_boundary_nodes[i].getDistanceToBoundary();
         counter++;
         for (j = 0; j < dim.getValue(); j++) {
            dbuffer[counter] = d_boundary_nodes[i].getNormalToBoundary(j);
            dbuffer[counter + 1] = d_boundary_nodes[i].getClosestBoundaryPoint(
                  j);
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
void
CutCell::unpackStream(
   tbox::MessageStream& stream,
   const hier::IntVector& offset)
{
   const tbox::Dimension& dim(d_index.getDim());

   int i, j, k;

   /*
    * #ints =
    *   d_index                (int[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)])
    */
   int int_buff_size = dim.getValue() + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += dim.getValue() * dim.getValue()
         * (dim.getValue() + 1 + dim.getValue() * dim.getValue());
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;
   stream.unpack(ibuffer, int_buff_size);
   pdat::CellIndex cell_index(dim);
   for (i = 0; i < dim.getValue(); i++) {
      cell_index(i) = ibuffer[counter];
      counter++;
   }
   d_index = cell_index + offset;

   int d_num_boundary_nodes = ibuffer[counter];
   counter++;

   if (s_enable_boundary_node_storage) {

      for (j = 0; j < d_num_boundary_nodes; j++) {
         pdat::NodeIndex node_index(dim);
         for (i = 0; i < dim.getValue(); i++) {
            node_index(i) = ibuffer[counter];
            counter++;
         }

         BoundaryNode bn(node_index);

         int iosb = ibuffer[counter];
         counter++;
         if (iosb == 1) bn.setNodeOnBoundary();

         for (k = 0; k < dim.getValue(); k++) {
            for (i = 0; i < dim.getValue(); i++) {
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
    *   d_normal          (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_front_centroid  (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */

   int dbl_buff_size = 3 + 4 * dim.getValue() + dim.getValue() + 2 + dim.getValue() * dim.getValue();
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += dim.getValue() * dim.getValue() * (1 + 2 * dim.getValue());
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;
   stream.unpack(dbuffer, dbl_buff_size);
   d_vol_fraction = dbuffer[0];
   d_surr_vol = dbuffer[1];
   d_front_area = dbuffer[2];
   counter += 3;

   for (i = 0; i < dim.getValue(); i++) {
      d_normal[i] = dbuffer[counter];
      d_front_centroid[i] = dbuffer[counter + 1];
      counter += 2;
   }

   for (i = 0; i < 2 * dim.getValue(); i++) {
      d_area_fraction[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < dim.getValue() + 2; i++) {
      d_flux_front[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < dim.getValue(); i++) {
      for (j = 0; j < dim.getValue(); j++) {
         d_newbase[i][j] = dbuffer[counter];
         counter++;
      }
   }

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is dim*dim
      for (i = 0; i < dim.getValue() * dim.getValue(); i++) {
         d_boundary_nodes[i].setDistanceToBoundary(dbuffer[counter]);
         counter++;
         for (j = 0; j < dim.getValue(); j++) {
            d_boundary_nodes[i].setNormalToBoundary(dbuffer[counter], j);
            d_boundary_nodes[i].setClosestBoundaryPoint(dbuffer[counter + 1], j);
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
 * i.e. IndexData<CutCell>.  They are used to write/read      *
 * CutCell data to/from the restart database.                       *
 *                                                                       *
 * The following writes data to the restart database.                    *
 *                                                                       *
 *************************************************************************
 */

void
CutCell::putToDatabase(
   tbox::Pointer<tbox::Database>& database)
{
   const tbox::Dimension& dim(d_index.getDim());

   int i, j, k;

   /*
    * #ints =
    *   d_index                (int[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)])
    */
   int int_buff_size = dim.getValue() + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += dim.getValue() * dim.getValue()
         * (dim.getValue() + 1 + dim.getValue() * dim.getValue());
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;
   for (i = 0; i < dim.getValue(); i++) {
      ibuffer[counter] = d_index(i);
      counter++;
   }

   ibuffer[counter] = d_num_boundary_nodes;
   counter++;

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is dim*dim
      for (j = 0; j < dim.getValue() * dim.getValue(); j++) {
         pdat::NodeIndex index = d_boundary_nodes[j].getIndex();
         for (i = 0; i < dim.getValue(); i++) {
            ibuffer[counter] = index(i);
            counter++;
         }
         bool on_boundary = d_boundary_nodes[j].getNodeOnBoundary();
         int iosb = 0;
         if (on_boundary) iosb = 1;
         ibuffer[counter] = iosb;
         counter++;

         // the maximum number of nearest neighbor nodes is dim
         for (k = 0; k < dim.getValue(); k++) {
            index = d_boundary_nodes[j].getNearestNeighborNode(k);
            for (i = 0; i < dim.getValue(); i++) {
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
    *   d_normal          (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_front_centroid  (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */

   int dbl_buff_size = 3 + 4 * dim.getValue() + dim.getValue() + 2 + dim.getValue() * dim.getValue();
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += dim.getValue() * dim.getValue() * (1 + 2 * dim.getValue());
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;

   dbuffer[0] = d_vol_fraction;
   dbuffer[1] = d_surr_vol;
   dbuffer[2] = d_front_area;
   counter += 3;

   for (i = 0; i < dim.getValue(); i++) {
      dbuffer[counter] = d_normal[i];
      dbuffer[counter + 1] = d_front_centroid[i];
      counter += 2;
   }

   for (i = 0; i < 2 * dim.getValue(); i++) {
      dbuffer[counter] = d_area_fraction[i];
      counter++;
   }

   for (i = 0; i < dim.getValue() + 2; i++) {
      dbuffer[counter] = d_flux_front[i];
      counter++;
   }

   for (i = 0; i < dim.getValue(); i++) {
      for (j = 0; j < dim.getValue(); j++) {
         dbuffer[counter] = d_newbase[i][j];
         counter++;
      }
   }

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is dim*dim
      for (i = 0; i < dim.getValue() * dim.getValue(); i++) {
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

void
CutCell::getFromDatabase(
   tbox::Pointer<tbox::Database>& database)
{
   const tbox::Dimension& dim(d_index.getDim());

   int i, j, k;
   /*
    * #ints =
    *   d_index                (int[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_num_boundary_nodes   (int) +
    *   d_boundary_nodes       (int[DIM*DIM*(DIM+1+DIM*DIM)])
    */
   int int_buff_size = dim.getValue() + 1;
   if (s_enable_boundary_node_storage) {
      int_buff_size += dim.getValue() * dim.getValue()
         * (dim.getValue() + 1 + dim.getValue() * dim.getValue());
   }
   int* ibuffer = new int[int_buff_size];
   int counter = 0;

   database->getIntegerArray("ibuffer", ibuffer, int_buff_size);
   pdat::CellIndex cell_index(dim);
   for (i = 0; i < dim.getValue(); i++) {
      cell_index(i) = ibuffer[counter];
      counter++;
   }
   d_index = cell_index;

   int d_num_boundary_nodes = ibuffer[counter];
   counter++;

   if (s_enable_boundary_node_storage) {
      for (j = 0; j < d_num_boundary_nodes; j++) {
         pdat::NodeIndex node_index(dim);
         for (i = 0; i < dim.getValue(); i++) {
            node_index(i) = ibuffer[counter];
            counter++;
         }

         BoundaryNode bn(node_index);

         int iosb = ibuffer[counter];
         counter++;
         if (iosb == 1) bn.setNodeOnBoundary();

         for (k = 0; k < dim.getValue(); k++) {
            for (i = 0; i < dim.getValue(); i++) {
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
    *   d_normal          (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_front_centroid  (double[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]) +
    *   d_area_fraction   (double[2*DIM]) +
    *   d_flux_front      (double[DIM+2]) +
    *   d_newbase         (double[DIM*DIM])
    *   d_boundary_nodes  (double[DIM*DIM*(1+2*DIM)])
    */

   int dbl_buff_size = 3 + 4 * dim.getValue() + dim.getValue() + 2 + dim.getValue() * dim.getValue();
   if (s_enable_boundary_node_storage) {
      dbl_buff_size += dim.getValue() * dim.getValue() * (1 + 2 * dim.getValue());
   }
   double* dbuffer = new double[dbl_buff_size];
   counter = 0;

   database->getDoubleArray("dbuffer", dbuffer, dbl_buff_size);
   d_vol_fraction = dbuffer[0];
   d_surr_vol = dbuffer[1];
   d_front_area = dbuffer[2];
   counter += 3;

   for (i = 0; i < dim.getValue(); i++) {
      d_normal[i] = dbuffer[counter];
      d_front_centroid[i] = dbuffer[counter + 1];
      counter += 2;
   }

   for (i = 0; i < 2 * dim.getValue(); i++) {
      d_area_fraction[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < dim.getValue() + 2; i++) {
      d_flux_front[i] = dbuffer[counter];
      counter++;
   }

   for (i = 0; i < dim.getValue(); i++) {
      for (j = 0; j < dim.getValue(); j++) {
         d_newbase[i][j] = dbuffer[counter];
         counter++;
      }
   }

   if (s_enable_boundary_node_storage) {
      // the maximum number of boundary nodes is dim*dim
      for (i = 0; i < dim.getValue() * dim.getValue(); i++) {
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
void
CutCell::initializeCutCellData()
{
   const tbox::Dimension& dim(d_index.getDim());

   if ((dim == tbox::Dimension(1)) || (dim > tbox::Dimension(3))) {
      TBOX_ERROR("CutCell : (dim == tbox::Dimension(1)) or > 3 not implemented");
   }

   d_vol_fraction = tbox::MathUtilities<double>::getSignalingNaN();
   d_surr_vol = tbox::MathUtilities<double>::getSignalingNaN();
   d_front_area = tbox::MathUtilities<double>::getSignalingNaN();

   int i;
   for (i = 0; i < dim.getValue(); i++) {
      d_normal[i] = tbox::MathUtilities<double>::getSignalingNaN();
      d_front_centroid[i] = tbox::MathUtilities<double>::getSignalingNaN();
   }

   int j;
   for (i = 0; i < dim.getValue(); i++) {
      for (j = 0; j < dim.getValue(); j++) {
         d_newbase[i][j] = tbox::MathUtilities<double>::getSignalingNaN();
      }
   }

   for (i = 0; i < 2 * dim.getValue(); i++) {
      d_area_fraction[i] = tbox::MathUtilities<double>::getSignalingNaN();
   }

   for (i = 0; i < dim.getValue() + 2; i++) {
      d_flux_front[i] = 0.;
   }

   d_num_boundary_nodes = 0;
   if (s_enable_boundary_node_storage) {
      d_boundary_nodes.clear();
      d_boundary_nodes.resizeArray(dim.getValue() * dim.getValue(), BoundaryNode(dim));
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

void
CutCell::copyCutCellData(
   const appu::CutCell& cut_cell)
{
   const tbox::Dimension& dim(d_index.getDim());

   d_vol_fraction = cut_cell.d_vol_fraction;
   d_surr_vol = cut_cell.d_surr_vol;
   d_front_area = cut_cell.d_front_area;

   int i;
   for (i = 0; i < dim.getValue(); i++) {
      d_normal[i] = cut_cell.d_normal[i];
      d_front_centroid[i] = cut_cell.d_front_centroid[i];
      for (int j = 0; j < dim.getValue(); j++) {
         d_newbase[i][j] = cut_cell.d_newbase[i][j];
      }
   }
   for (i = 0; i < 2 * dim.getValue(); i++) {
      d_area_fraction[i] = cut_cell.d_area_fraction[i];
   }

   for (i = 0; i < dim.getValue() + 2; i++) {
      d_flux_front[i] = cut_cell.d_flux_front[i];
   }

   d_num_boundary_nodes = cut_cell.d_num_boundary_nodes;
   if (s_enable_boundary_node_storage) {
      d_boundary_nodes.resizeArray(dim.getValue() * dim.getValue(), BoundaryNode(dim));
      for (i = 0; i < d_num_boundary_nodes; i++) {
         d_boundary_nodes[i] = cut_cell.d_boundary_nodes[i];
      }
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
