//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/BoundaryBox.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Box representing a portion of the AMR index space
//

#ifndef included_hier_BoundaryBox
#define included_hier_BoundaryBox

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif


namespace SAMRAI {
   namespace hier {


/**
 * Class BoundaryBox<DIM> is a class that is used to facilitate the filling
 * of ghost layers around a patch.  Objects of this type are created to be
 * part of a PatchGeometry<DIM> object.  The BoundaryBox<DIM> consists
 * of a Box<DIM>, a boundary type, and a location index.  The Box<DIM> is of
 * width 1 in at least
 * one dirction and is located just outside of a patch boundary, intersecting
 * the patch boundary's face, edge, or corner (node).  The boundary type
 * location index indicates the location of the boundary box in relation to the
 * location of the patch. 
 *
 * @see hier::Box
 * @see hier::PatchGeometry
 * @see hier::BoundaryLookupTable
 */

template<int DIM> class BoundaryBox
{
public:
   /**
    * The default constructor creates an undefined boundary box, with
    * invalid values.
    */
   BoundaryBox();

   /**
    * Create a boundary box from a Box<DIM> and integers that indicate
    * the boundary type and the location index
    */
   BoundaryBox(const Box<DIM>& box, 
                     const int bdry_type,
                     const int location_index);

   /**
    * The copy constructor copies the data of the argument box.
    */
   BoundaryBox(const BoundaryBox<DIM>& boundary_box);

   /**
    * The destructor for BoundaryBox.
    */
   ~BoundaryBox<DIM>();

   /**
    * Return a reference to the Box<DIM> member of the boundary box
    */
   const Box<DIM>& getBox() const;

   /**
    * Return the boundary type (codimension) of the boundary box.
    *
    * \verbatim
    * Convention:
    * ===========
    *
    * 1d
    * --
    * 1 = node
    *
    * 2d
    * --
    * 1 = edge
    * 2 = node
    *
    * 3d
    * --
    * 1 = face
    * 2 = edge
    * 3 = node
    * \endverbatim
    */ 
   int getBoundaryType() const;

  /**
   * Return the location index for the boundary box.  The location index
   * is an integer which indicates the location of the boundary box in relation
   * to the location of the associated patch.  That is, the location index
   * tells whether the boundary box is ``above'' or ``below'' the patch in
   * each coordinate direction.  The conventions for the location index depend
   * on the dimension of the problem and the boundary type of the BoundaryBox.
   *
   * \verbatim
   * Conventions: 
   * ============
   *
   * 1d
   * --
   * node (codimension 1):
   * x_lo : 0
   * x_hi : 1
   *
   * 2d
   * --
   * edge (codimension 1):
   * x_lo: 0
   * x_hi: 1
   * y_lo: 2
   * y_hi: 3
   *
   * node (codimension 2):
   * x_lo, y_lo: 0
   * x_hi, y_lo: 1
   * x_lo, y_hi: 2
   * x_hi, y_hi: 3
   *
   * 3d
   * --
   *
   * face (codimension 1):
   * x_lo: 0
   * x_hi: 1
   * y_lo: 2
   * y_hi: 3
   * z_lo: 4
   * z_hi: 5
   *
   * edge (codimension 2):
   * y_lo, z_lo: 0
   * y_hi, z_lo: 1
   * y_lo, z_hi: 2
   * y_hi, z_hi: 3
   * x_lo, z_lo: 4
   * x_lo, z_hi: 5
   * x_hi, z_lo: 6
   * x_hi, z_hi: 7
   * x_lo, y_lo: 8
   * x_hi, y_lo: 9
   * x_lo, y_hi: 10
   * x_hi, y_hi: 11
   *
   * node (codimension 3):
   * x_lo, y_lo, z_lo: 0
   * x_hi, y_lo, z_lo: 1
   * x_lo, y_hi, z_lo: 2
   * x_hi, y_hi, z_lo: 3
   * x_lo, y_lo, z_hi: 4
   * x_hi, y_lo, z_hi: 5
   * x_lo, y_hi, z_hi: 6
   * x_hi, y_hi, z_hi: 7
   *
   * \endverbatim 
   */
   int getLocationIndex() const;

   /**
    * Set the multiblock singularity flag to the argument value.
    */
   void setIsMultiblockSingularity(bool is_mblk_singularity);

   /**
    * Get the value of the multiblock singularity flag.  This should
    * always return false when not running a problem on a multiblock domain.
    */
   bool getIsMultiblockSingularity() const;

   /**
    * The assignment operator copies all data components.
    */
   BoundaryBox<DIM>& operator=(const BoundaryBox<DIM>& boundary_box);

private:

   Box<DIM> d_box;
   int d_bdry_type;
   int d_location_index;
   bool d_is_mblk_singularity;
};


}
}

#ifndef DEBUG_NO_INLINE
#include "BoundaryBox.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoundaryBox.C"
#endif
