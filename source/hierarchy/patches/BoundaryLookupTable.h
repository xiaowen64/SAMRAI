//
// File:	BoundaryLookupTable.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	 Lookup table to aid in BoundaryBox construction
//

#ifndef included_hier_BoundaryLookupTable
#define included_hier_BoundaryLookupTable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif

namespace SAMRAI {
   namespace hier {


/*!
 * Class BoundaryLookupTable<DIM> is a singleton class that contains a table
 * that organizes all of the possible directions where a physical boundary
 * can exist around a domain.  It is used by GridGeometry<DIM> during
 * the construction of boundary boxes and by PatchGeometry<DIM> to determine
 * the box regions that need to be filled during a physical boundary fill.
 *
 * @see hier::BoundaryBox
 * @see hier::GridGeometry
 * @see hier::PatchGeometry
 */

template<int DIM> class BoundaryLookupTable 
{
public:
   /*!
    * @brief Return pointer to singleton instance of the boundary lookup table. 
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
    * @brief Set the lookup table to return original numbering scheme
    *
    * For codimension 2 in 3-dimensional problems, the numbering scheme
    * for the location indices of BoundaryBox has been changed.  To use
    * the original numbering scheme for backward compatibility, call this
    * static function with the argument set to true;
    *
    * @param use_original bool argument set to true if using original scheme
    */ 
   static void setUsingOriginalLocations(const bool use_original);
   
   /*!
    * @brief Return array of active directions for specific case.
    *
    * Returns integer array of length codim of the active directions for 
    * the boundary of codimension codim indexed by loc.
    *
    * @param loc Location index being used
    * @param codim Codimension being used
    */
   const tbox::Array<int>& getDirections(const int loc, const int codim) const;

   /*!
    * @brief Returns array containd maximum locations for each codimension.
    *
    * Returns integer array of length DIM of the maximum limits, for each 
    * codimension, of the location indicies.
    */
   const tbox::Array<int>& getMaxLocationIndices() const;

   /*!
    * @brief Determines if boundary is lower boundary
    *
    * Returns true if the boundary type of codimension codim indexed by loc
    * is a lower boundary in the specified dimension.
    *
    * @param loc Location index of boundary being tested
    * @param codim Codimension of boundary being tested
    * @param dim dimension identifier
    */
   bool isLower(const int loc, const int codim, const int dim) const;

   /*!
    * @brief Determines if boundary is upper boundary
    *
    * Returns true if the boundary type of codimension codim indexed by loc
    * is a upper boundary in the specified dimension.
    *
    * @param loc Location index of boundary being tested
    * @param codim Codimension of boundary being tested
    * @param dim dimension identifier
    */
   bool isUpper(const int loc, const int codim, const int index) const;

   /*!
    * @brief execute the mapping between original numbering and new scheme
    *
    * For codimension 2 in 3 dimensions, maps the value of the argument
    * from the original numbering scheme to the new scheme, or vice versa,
    * and returns the mapped value.
    *
    * @param loc location index to be mapped
    */
   int mapLocationIndex(const int loc) const;

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
   void buildTable(int *table, int codim, int ibeg);

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
