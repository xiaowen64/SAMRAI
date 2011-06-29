/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Cut cell class for embedded boundary implementations 
 *
 ************************************************************************/

#ifndef included_appu_CutCell
#define included_appu_CutCell

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryNode.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Database.h"

namespace SAMRAI {
namespace appu {

/*!
 * @brief The CutCell struct holds data and methods to define a cut-cell
 * on an irregular boundary.
 *
 * Information maintained by the struct includes the following:
 *
 *
 *     - INDEX                 - index (i,j,k) of the cell
 *     - VOL FRACTION          - volume fraction of the cell
 *     - AREA FRACTION[2*DIM]  - areas of the faces of the cell
 *     - NORMAL[tbox::Dimension::MAXIMUM_DIMENSION_VALUE]           - normal to the cut plane through the cell
 *     - FRONT AREA            - area of exposed plane cutting the cell
 *     - SURROUNDING VOL       - volume of cells surrounding the cut cell
 *     - NEW BASE              - a basis aligned with the cell cut face
 *     - BOUNDARY NODE[DIM**2] - nodes on the corners of the cut cell, which
 *                               contain information for immersed boundary
 *                               calculations
 */

class CutCell
{
public:
   /*!
    * Static function to set whether storage will be allocated for boundary
    * node information.
    *
    */
   static void
   enableBoundaryNodeStorage();

   /*!
    * Static function that returns whether storage is allocated for boundary
    * node information.
    */
   static bool
   boundaryNodesEnabled();

   /*!
    * Create a new ``empty'' CutCell.
    */
   CutCell();

   /*!
    * Create a new ``empty'' CutCell.
    */
   explicit CutCell(
      const tbox::Dimension& dim);

   /*!
    * Create a new cut cell with specified cell index.
    */
   explicit CutCell(
      const pdat::CellIndex& cut_cell);

   /*!
    * The copy constructor copies the data of the argument cell.
    */
   CutCell(
      const appu::CutCell& cut_cell);

   /*!
    * The assignment operator copies the data of the argument cell.
    */
#if 0
   CutCell&
   operator = (
      const appu::CutCell& cut_cell);
#endif
   CutCell&
   copy(
      const appu::CutCell& cut_cell);

   /*!
    * The destructor for CutCell.
    */
   ~CutCell();

   /*!
    * Returns the index (i,j,k) of the cell.
    */
   pdat::CellIndex
   getIndex() const;

   /*!
    * Returns the volume fraction of cell.
    */
   double
   getVolume() const;

   /*!
    * Returns the area fraction of the cell (dimension 2*DIM)
    *   0,1 - Xlower,Xupper
    *   2,3 - Ylower,Yupper
    *   4,5 - Zlower,Zupper
    */
   const double *
   getArea() const;

   /*!
    * Returns the area fraction of the cell for face i
    *   i = 0,1 - Xlower,Xupper
    *   i = 2,3 - Ylower,Yupper
    *   i = 4,5 - Zlower,Zupper
    */
   double
   getArea(
      const int i) const;

   /*!
    * Returns the normal vector (dimension DIM).
    */
   const double *
   getNormal() const;

   /*!
    * Returns the ith component of the normal vector.
    */
   double
   getNormal(
      const int i) const;

   /*!
    * Returns the frontal area of exposed cut surface.
    */
   double
   getFrontArea() const;

   /*!
    * Returns the front centroid vector (dimension DIM).
    */
   const double *
   getFrontCentroid() const;

   /*!
    * Returns the front centroid location for direction i.
    */
   double
   getFrontCentroid(
      const int i) const;

   /*!
    * Returns the volume of the cells surrounding the cut cell.
    */
   double
   getSurrVolume() const;

   /*!
    * Returns the new base with components (i,j).
    */
   double
   getNewBase(
      const int i,
      const int j) const;

   /*!
    * Return number of boundary nodes.
    */
   int
   getNumberOfBoundaryNodes() const;

   /*!
    * Return boundary node class for location i.
    */
   BoundaryNode
   getBoundaryNode(
      const int i) const;

   /*!
    * TO BE REMOVED (eventually)
    */
   double
   getFluxFront(
      const int m) const;

   /*!
    * Sets the volume fraction of cell.
    */
   void
   setVolume(
      const double volume);

   /*!
    * Sets the area fraction of the cell for face i
    *   i = 0,1 - Xlower,Xupper
    *   i = 2,3 - Ylower,Yupper
    *   i = 4,5 - Zlower,Zupper
    */
   void
   setArea(
      const double area,
      const int i);

   /*!
    * Sets the normal vector for dimension i.
    */
   void
   setNormal(
      const double normal,
      const int i);

   /*!
    * Sets the normal vector (dimension DIM).
    */
   void
   setNormal(
      const double* normal);

   /*!
    * Sets the volume of the cells surrounding the cut cell.
    */
   void
   setSurrVolume(
      const double surrvol);

   /*!
    * Sets the front area.
    */
   void
   setFrontArea(
      const double area);

   /*!
    * Sets the front centroid location for direction i.
    */
   void
   setFrontCentroid(
      const double centroid,
      const int i);

   /*!
    * Set split cell.
    */
   void
   setSplit();

   /*!
    * Explicitly set the new base vector with components (i,j).
    */
   void
   setNewBase(
      const double base,
      const int i,
      const int j);

   /*!
    * Set new base vector using supplied "base" vector.
    */
   void
   setNewBase(
      const double* base);

   /*!
    * Set new base vector using cell's normal vector.
    */
   void
   setNewBase();

   /*!
    * TO BE REMOVED (eventually)
    */
   void
   setFluxFront(
      const int m,
      const double front);

   /*!
    * Add the supplied boundary node to the array of boundary nodes
    * maintained by this cut cell.
    */
   void
   setBoundaryNode(
      const BoundaryNode& node);

   /*!
    * Overwrite index i of the boundary node array maintained by
    * the cut cell with the supplied boundary node.
    */
   void
   setBoundaryNode(
      const BoundaryNode& node,
      const int i);

   /*!
    * Add the supplied index as an uninitialized boundary node.
    */
   void
   setBoundaryNode(
      const pdat::NodeIndex& node);

   /*!
    * Print volume and area data for the cell.
    */
   void
   printVolumeAndAreas(
      std::ostream& os) const;

   /*!
    * Print normal data for the cell.
    */
   void
   printNormal(
      std::ostream& os) const;

   /*!
    * Print boundary nodes for the cell.
    */
   void
   printBoundaryNodes(
      std::ostream& os) const;

   /*!
    * Print all data in the struct.
    */
   void
   printAll(
      std::ostream& os) const;

   /*!
    * The copySourceItem() method allows CutCell to be a templated
    * data type for IndexData - i.e. IndexData<CutCell>.  In
    * addition to this method, the other methods that must be defined are
    * getDataStreamSize(), packStream(), unpackStream() for communication,
    * putToDatabase(), getFromDatabase for restart.  These are described
    * below.
    */
   void
   copySourceItem(
      hier::Index& index,
      const hier::IntVector& src_offset,
      appu::CutCell& src_item);

   /*!
    * The following functions enable parallel communication with
    * CutCells. They are used in SAMRAI communication infrastructure to
    * specify the number of bytes of data stored in each CutCell object,
    * and to pack and unpack the data to the specified stream.
    */
   size_t
   getDataStreamSize();
   void
   packStream(
      tbox::MessageStream& stream);
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::IntVector& offset);

   /*!
    * These functions are used to read/write CutCell data to/from
    * restart.
    */
   void
   getFromDatabase(
      tbox::Pointer<tbox::Database>& database);
   void
   putToDatabase(
      tbox::Pointer<tbox::Database>& database);

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int CUTCELL_VERSION;

   /*
    * Specifies whether or not storage will be allocated for boundary node
    * data.
    */
   static bool s_enable_boundary_node_storage;

   /*
    * Initialize data in a new cut cell.
    */
   void
   initializeCutCellData();

   /*
    * Copy data from supplied cut cell.
    */
   void
   copyCutCellData(
      const appu::CutCell& cut_cell);

   /*
    * Index of CutCell
    */
   pdat::CellIndex d_index;

   /*
    * Volume and area fractions.
    */
   double d_vol_fraction;
   double d_area_fraction[2 * tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Surface normal and exposed (wetted) area of cut region.
    */
   double d_normal[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double d_front_area;
   double d_front_centroid[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   bool d_is_split;

   /*
    * Volume of surrounding cells.
    */
   double d_surr_vol;

   /*
    * New base.  Not sure exactly if we need to store this here, but
    * keep it for now.
    */
   double d_newbase[tbox::Dimension::MAXIMUM_DIMENSION_VALUE][tbox::Dimension::
                                                        MAXIMUM_DIMENSION_VALUE];

   /*
    * Array of BoundaryNodes on this cut cell.
    */
   int d_num_boundary_nodes;
   tbox::Array<appu::BoundaryNode> d_boundary_nodes;

   /*
    * Flux front - EVENTUALLY REMOVE!!
    */
   double d_flux_front[2 + tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

};

}
}
#endif
