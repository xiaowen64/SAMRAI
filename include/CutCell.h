//
// File:        CutCell.h
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 726 $
// Modified:    $Date: 2005-11-10 17:19:40 -0800 (Thu, 10 Nov 2005) $
// Description: Cut cell class for embedded boundary implementations
//

#ifndef included_appu_CutCell
#define included_appu_CutCell

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_AbstractStream
#include "tbox/AbstractStream.h"
#endif
#ifndef included_appu_BoundaryNode
#include "BoundaryNode.h"
#endif
#ifndef included_pdat_CellIndex
#include "CellIndex.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_hier_Index
#include "Index.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_IOStream
#include "tbox/IOStream.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif


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
 *     - NORMAL[DIM]           - normal to the cut plane through the cell
 *     - FRONT AREA            - area of exposed plane cutting the cell
 *     - SURROUNDING VOL       - volume of cells surrounding the cut cell
 *     - NEW BASE              - a basis aligned with the cell cut face
 *     - BOUNDARY NODE[DIM**2] - nodes on the corners of the cut cell, which
 *                               contain information for immersed boundary
 *                               calculations
 */


template<int DIM> class CutCell
{
public:
   /*!
    * Static function to set whether storage will be allocated for boundary 
    * node information.
    *
    * @param  enable_bn bool indicating whether boundary node data storage
    *                   should be enabled.
    */
   static void enableBoundaryNodeStorage();

   /*!
    * Static function that returns whether storage is allocated for boundary 
    * node information.
    */
   static bool boundaryNodesEnabled();

   /*!
    * Create a new ``empty'' CutCell.
    */
   CutCell();

   /*!
    * Create a new cut cell with specified cell index.
    */
   CutCell(const pdat::CellIndex<DIM>& cut_cell);

   /*!
    * The copy constructor copies the data of the argument cell.
    */
   CutCell(const appu::CutCell<DIM>& cut_cell);
 
   /*!
    * The assignment operator copies the data of the argument cell.
    */
#if 0
   CutCell& operator=(const appu::CutCell<DIM>& cut_cell);
#endif
   CutCell& copy(const appu::CutCell<DIM>& cut_cell);

   /*!
    * The destructor for CutCell.
    */
   ~CutCell();

   /*!
    * Returns the index (i,j,k) of the cell.
    */
   pdat::CellIndex<DIM> getIndex() const;

   /*!
    * Returns the volume fraction of cell.
    */
   double getVolume() const;


   /*!
    * Returns the area fraction of the cell (dimension 2*DIM)
    *   0,1 - Xlower,Xupper
    *   2,3 - Ylower,Yupper
    *   4,5 - Zlower,Zupper
    */
   const double* getArea() const;

   /*!
    * Returns the area fraction of the cell for face i 
    *   i = 0,1 - Xlower,Xupper
    *   i = 2,3 - Ylower,Yupper
    *   i = 4,5 - Zlower,Zupper
    */
   double getArea(const int i) const;

   /*!
    * Returns the normal vector (dimension DIM).
    */
   const double* getNormal() const;
   
   /*!
    * Returns the ith component of the normal vector.
    */
   double getNormal(const int i) const;

   /*!
    * Returns the frontal area of exposed cut surface.
    */
   double getFrontArea() const;

   /*!
    * Returns the front centroid vector (dimension DIM).
    */
   const double* getFrontCentroid() const;

   /*!
    * Returns the front centroid location for direction i.
    */
   double getFrontCentroid(const int i) const;
   
   /*!
    * Returns the volume of the cells surrounding the cut cell.
    */
   double getSurrVolume() const;


   /*!
    * Returns the new base with components (i,j).
    */
   double getNewBase(const int i, const int j) const;

   /*!
    * Return number of boundary nodes.
    */
   int getNumberBoundaryNodes() const;

   /*!
    * Return boundary node class for location i.
    */
   BoundaryNode<DIM> getBoundaryNode(const int i) const;

   /*!
    * TO BE REMOVED (eventually)
    */
   double getFluxFront(const int m) const;

   /*!
    * Sets the volume fraction of cell.
    */
   void setVolume(const double volume);

   /*!
    * Sets the area fraction of the cell for face i 
    *   i = 0,1 - Xlower,Xupper
    *   i = 2,3 - Ylower,Yupper
    *   i = 4,5 - Zlower,Zupper
    */
   void setArea(const double area,
                const int i);
   
   /*!
    * Sets the normal vector for dimension i.
    */
   void setNormal(const double normal,
                  const int i);

   /*!
    * Sets the normal vector (dimension DIM).
    */
   void setNormal(const double* normal);

   /*!
    * Sets the volume of the cells surrounding the cut cell.
    */
   void setSurrVolume(const double surrvol);

   /*!
    * Sets the front area.
    */
   void setFrontArea(const double area);

   /*!
    * Sets the front centroid location for direction i.
    */
   void setFrontCentroid(const double centroid,
                         const int i);

   /*!
    * Set split cell.
    */
   void setSplit();

   /*!
    * Explicitly set the new base vector with components (i,j).
    */
   void setNewBase(const double base,
                   const int i, 
                   const int j);

   /*!
    * Set new base vector using supplied "base" vector.
    */
   void setNewBase(const double* base);

   /*!
    * Set new base vector using cell's normal vector.
    */
   void setNewBase();

   /*!
    * TO BE REMOVED (eventually)
    */
   void setFluxFront(const int m, const double front);
   
   /*!
    * Add the supplied boundary node to the array of boundary nodes
    * maintained by this cut cell.
    */
   void setBoundaryNode(const BoundaryNode<DIM>& node);

   /*!
    * Overwrite index i of the boundary node array maintained by 
    * the cut cell with the supplied boundary node.
    */
   void setBoundaryNode(const BoundaryNode<DIM>& node,
                        const int i);

   /*!
    * Add the supplied index as an uninitialized boundary node. 
    */
   void setBoundaryNode(const pdat::NodeIndex<DIM>& node);

   /*!
    * Print volume and area data for the cell.
    */
   void printVolumeAndAreas(ostream& os) const;

   /*!
    * Print normal data for the cell.
    */
   void printNormal(ostream& os) const;

   /*!
    * Print boundary nodes for the cell.
    */
   void printBoundaryNodes(ostream& os) const;

   /*!
    * Print all data in the struct.
    */
   void printAll(ostream& os) const;

   /*!
    * The copySourceItem() method allows CutCell to be a templated
    * data type for IndexData - i.e. IndexData<DIM,CutCell<DIM>>.  In
    * addition to this method, the other methods that must be defined are
    * getDataStreamSize(), packStream(), unpackStream() for communication,
    * putToDatabase(), getFromDatabase for restart.  These are described
    * below. 
    */
   void copySourceItem(hier::Index<DIM>& index, 
                       const hier::IntVector<DIM>& src_offset, 
                       appu::CutCell<DIM>& src_item);

   /*!
    * The following functions enable parallel communication with 
    * CutCells. They are used in SAMRAI communication infrastructure to 
    * specify the number of bytes of data stored in each CutCell object,
    * and to pack and unpack the data to the specified stream.
    */
   size_t getDataStreamSize();
   void packStream(tbox::AbstractStream& stream);
   void unpackStream(tbox::AbstractStream& stream, 
                     const hier::IntVector<DIM>& offset);
   

   /*!
    * These functions are used to read/write CutCell data to/from 
    * restart.
    */
   void getFromDatabase(tbox::Pointer<tbox::Database>& database);
   void putToDatabase(tbox::Pointer<tbox::Database>& database);
   
 
private:

   /*
    * Specifies whether or not storage will be allocated for boundary node
    * data.
    */
   static bool s_enable_boundary_node_storage;

   /*
    * Initialize data in a new cut cell.
    */
   void initializeCutCellData();

   /*
    * Copy data from supplied cut cell.
    */
   void copyCutCellData(const appu::CutCell<DIM>& cut_cell);

   /*
    * Index of CutCell
    */
   pdat::CellIndex<DIM> d_index;

   /*
    * Volume and area fractions.
    */
   double d_vol_fraction;
   double d_area_fraction[2*DIM];

   /*
    * Surface normal and exposed (wetted) area of cut region.
    */
   double d_normal[DIM];
   double d_front_area;
   double d_front_centroid[DIM];
   
   bool d_is_split;
   

   /*
    * Volume of surrounding cells.
    */
   double d_surr_vol;

   /*
    * New base.  Not sure exactly if we need to store this here, but
    * keep it for now.
    */
   double d_newbase[DIM][DIM];

   /*
    * Array of BoundaryNodes on this cut cell.
    */
   int d_num_boundary_nodes;
   tbox::Array< appu::BoundaryNode<DIM> > d_boundary_nodes;

   /*
    * Flux front - EVENTUALLY REMOVE!!
    */
   double d_flux_front[2+DIM];

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CutCell.C"
#endif

