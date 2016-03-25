//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/patches/BoundaryLookupTable.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	 Lookup table to aid in BoundaryBox construction
//

#ifndef included_hier_BoundaryLookupTable
#define included_hier_BoundaryLookupTable

#include "SAMRAI_config.h"
#include "IntVector.h"
#include "tbox/Array.h"

namespace SAMRAI {
   namespace hier {


/*!
 * Class BoundaryLookupTable<DIM> is a singleton class that maintains 
 * a table that organizes all of the possible boundary region cases 
 * for a patch.  It is used primarily by the hier::GridGeometry<DIM> 
 * during the construction of physical boundary boxes for patches and
 * by the hier::PatchGeometry<DIM> class to determine box regions 
 * to be filled during a physical boundary fill.
 *
 * This class is useful for any situation where enumerating the
 * cases for boundary regions around a box is needed. The main advantage
 * of using this class is that such calculations can be programmed in
 * a dimension-independent way.
 *
 * @see hier::BoundaryBox
 * @see hier::GridGeometry
 * @see hier::PatchGeometry
 */

template<int DIM> class BoundaryLookupTable 
{
public:
   /*!
    * @brief Return pointer to singleton instance of the boundary 
    * lookup table. 
    *
    * Note that when the database is accessed for the first time, the
    * Singleton instance is registered with the ShutdownRegistry
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    * 
    * @return  tbox::Pointer to lookup table instance.
    */
   static BoundaryLookupTable<DIM>* getLookupTable();

   /*!
    * @brief Deallocate the BoundaryLookupTable<DIM> instance.
    *
    * It is not necessary to call this function at program termination,
    * since it is automatically called by the ShutdownRegistry class.
    */
   static void freeLookupTable();

   /*!
    * @brief Set the lookup table to use boundary box location
    * numbering scheme used in SAMRAI prior to version 2.0, when
    * all spatial geometry dependent classes became templated on
    * the spatial dimension.
    *
    * To use the older numbering scheme for backward compatibility, 
    * call this function with the argument set to true.  Otherwise,
    * the new numbering scheme is applied.
    *
    * In the newer (version 2.0 and beyond) the location numbering 
    * scheme for boundary boxes of codimension 2 in 3 spatial 
    * dimensions changed so that the numbering scheme generalizes
    * consistently to all spatial dimensions. 
    * 
    * @param use_original bool argument true if original location
    *                     index numbering scheme is desired
    */ 
   static void setUsingOriginalLocations(bool use_original);
   
   /*!
    * @brief Get array of active directions for specific boundary 
    * region case.  Such active directions refer to those coordinate 
    * directions in which the boundary region would have to be shifted 
    * to be contained in the corresponding box region (whose boundary
    * we are interested in).
    *
    * @return  const reference to integer array of length codim 
    *          containing the active directions for boundary case.
    *
    * @param loc integer location index of boundary region
    * @param codim integer codimension of boundary region
    */
   const tbox::Array<int>& getDirections(int loc, int codim) const;

   /*!
    * @brief Get array of maximum number of locations for each 
    * codimension boundary case. 
    *
    * @return integer array of length DIM, each entry of which indicates
    *         the maximum number of boundary locations for each 
    *         codimension
    */
   const tbox::Array<int>& getMaxLocationIndices() const;

   /*!
    * @brief Determines if given boundary information indicates a
    * a lower boundary region (i.e., the associated box region 
    * contains higher values along the axis in the coordinate 
    * direction than the boundary region).
    *
    * @return bool true if the boundary type of codimension codim indexed 
    * by loc is a lower boundary in the specified dimension; 
    * return false if the boundary is an upper boundary.
    *
    * @param loc integer location index of boundary region
    * @param codim integer codimension of boundary region
    * @param index integer spatial dimension identifier
    */
   bool isLower(int loc, int codim, int index) const;

   /*!
    * @brief Determines if given boundary information indicates a
    * an upper boundary region (i.e., the associated box region 
    * contains lower values along the axis in the coordinate 
    * direction than the boundary region).
    *
    * @return bool true if the boundary type of codimension codim indexed 
    * by loc is an upper boundary in the specified dimension; 
    * return false if the boundary is a lower boundary.
    *
    * @param loc integer location index of boundary region
    * @param codim integer codimension of boundary region
    * @param index integer spatial dimension identifier
    */
   bool isUpper(int loc, int codim, int index) const;

   /*!
    * @brief Map boundary box location index between the original 
    * numbering and new scheme.  In particular, for codimension 2 
    * in 3 spatial dimensions, the value of the argument is mapped 
    * from the original numbering scheme to the new scheme, or vice versa.
    *
    * @return integer location index for other numbering scheme.
    * 
    * @param loc integer location index of boundary region
    */
   int mapLocationIndex(int loc) const;

   /*!
    * @brief Get array of boundary direction IntVectors.
    *
    * For any codimension, there is a particular number of valid boundary
    * locations.  This function returns an array of IntVectors that provide
    * information about where each boundary location lies in relation to
    * a patch.  The array's length is the number of valid locations for
    * the given codimension, and the array is indexed by the location id's
    * that are set up by this BoundaryLookupTable class.
    *
    * For a particular location, each element of the IntVector tells whether
    * the location is on the lower or upper side, or neither, of the patch in
    * a specific coordinate direction.  -1 indicates lower, 0 indicates
    * neither, and 1 indicates upper.
    *
    * @return       Array of IntVectors, one element for each valid location 
    * @param codim  codimension
    */ 
   const tbox::Array< IntVector<DIM> >& getBoundaryDirections(int codim) const;

protected:
   /**
    * The constructor for BoundaryLookupTable<DIM> is protected. 
    * Consistent with the definition of a Singleton class, only the 
    * database object has access to the constructor for the class. 
    *
    * The constructor initializes the state of lookup table contents.
    */
   BoundaryLookupTable();

   /**
    * The destructor for BoundaryLookupTable<DIM> is protected. See the
    * comments for the constructor.
    *
    * The destructor deallocates lookup table contents.
    */
   ~BoundaryLookupTable<DIM>();

private:

   /*
    * Private member function that recursively computes the entries in the
    * lookup table for a given codimension.
    */
   void buildTable(int *table, int codim, int ibeg, int (&work)[DIM], int &lvl, int *&ptr);

   /*
    * Private member function that builds table of intvectors that stores
    * whether a boundary location is upper, lower, or neither in each
    * coordinate direction
    */
   void buildBoundaryDirectionVectors();

   /*
    * Static data members used to control access to and destruction of
    * singleton variable database instance.
    */
   static BoundaryLookupTable<DIM>* s_lookup_table_instance;
   static bool s_registered_callback;

   /*
    * Static member that tells whether original location index scheme is
    * being used
    */
   static bool s_using_original_locations;

   /*
    * Data member array used to store the number of combinations for each
    * codimension.
    */
   tbox::Array<int> d_ncomb;   

   /*
    * Data member array use to store the number of possible location indices
    * for each codimension.
    */
   tbox::Array<int> d_max_li;   

   /*
    * Data member used to store the lookup table.
    */
   tbox::Array< tbox::Array<int> > d_table[DIM];

   tbox::Array< tbox::Array< hier::IntVector<DIM> > > d_bdry_dirs; 
};

}
}

#ifndef DEBUG_NO_INLINE
#include "BoundaryLookupTable.I"
#endif

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoundaryLookupTable.C"
#endif
