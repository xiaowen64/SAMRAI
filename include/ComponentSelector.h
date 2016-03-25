//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/ComponentSelector.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Simple bit vector of a fixed length (128 bits)
//

#ifndef included_hier_ComponentSelector
#define included_hier_ComponentSelector

#include "SAMRAI_config.h"
#include "tbox/DescribedClass.h"
#include "tbox/PIO.h"

namespace SAMRAI {
   namespace hier {

/*!
 * @brief Class ComponentSelector implements a simple bit vector of a fixed
 * length and is typically used to apply operations on subsets of entries 
 * in the patch data array owned by a patch (e.g., allocate/deallocate).  
 * All ComponentSelector objects have the same bit vector length that is
 * established by the SAMRAIManager utility. See the documentation
 * of the SAMRAIManager utility for information about changing this 
 * maximum value.
 * 
 * @see tbox::SAMRAIManager
 */

class ComponentSelector : public tbox::DescribedClass
{
public:
   /*!
    * Create a component selector and initialize all bits to the specified 
    * default boolean flag value.  If no default value is provided, then all 
    * bits are set to false.
    */
   ComponentSelector(const bool flag = false);

   /*!
    * Copy construct a create a component selector identical to the argument.
    */
   ComponentSelector(const ComponentSelector& flags);

   /*!
    * The destructor for a component selector does nothing interesting.
    */
   ~ComponentSelector();

   /*!
    * Return total number of flags (i.e., bits) in this component selector.
    */
   int getSize() const;

   /*!
    * Set all bit settings in this component selector to those in the
    * argument component selector.
    */
   ComponentSelector& operator=(const ComponentSelector& flags);

   /*!
    * Compare two component selectors for equality in all bit positions.
    * If all bits are logically equal, then the return value is true;
    * otherwise, the return value is false.
    */
   bool operator==(const ComponentSelector& flags) const;

   /*!
    * Compare two component selectors for inequality in any bit position.
    * If any two bits are logically unequal, then true is returned; 
    * otherwise, false is returned.
    */
   bool operator!=(const ComponentSelector& flags) const;

   /*!
    * Generate and return a component selector set to the bitwise logical OR 
    * of this component selector and the argument component selector.
    */
   ComponentSelector operator|(const ComponentSelector& flags) const;

   /*!
    * Generate and return a component selector set to the bitwise logical AND
    * of this component selector and the argument component selector.
    */
   ComponentSelector operator&(const ComponentSelector& flags) const;

   /*!
    * Generate and return a component selector set to the bitwise logical
    * negation of this component selector.
    */
   ComponentSelector operator!() const;

   /*!
    * Generate and return a component selector set to the bitwise logical
    * AND of this component selector and the bitwise logical negation of 
    * the argument component selector.
    */
   ComponentSelector andNot(const ComponentSelector& flags) const;

   /*!
    * Set all bit settings in this component selector to the bitwise logical OR
    * of this component selector and the argument component selector.
    */
   ComponentSelector& operator|=(const ComponentSelector& flags);

   /*!
    * Set all bit settings in this component selector to the bitwise logical AND
    * of this component selector and the argument component selector.
    */
   ComponentSelector& operator&=(const ComponentSelector& flags);

   /*!
    * Check whether the specified bit vector position is true.  If so, 
    * return true; otherwise, return false.
    *
    * When assertion checking is active, an assertion will result if the
    * given position is out-of-bounds.
    */
   bool isSet(const int i) const;

   /*!
    * Set the specified bit vector position to true.  
    *
    * When assertion checking is active, an assertion will result if the
    * given position is out-of-bounds.
    */
   void setFlag(const int i);

   /*!
    * Set the specified bit vector position to false.
    *
    * When assertion checking is active, an assertion will result if the
    * given position is out-of-bounds.
    */
   void clrFlag(const int i);

   /*!
    * Set all bit vector positions to true.
    */
   void setAllFlags();

   /*!
    * Set all bit vector positions to false.
    */
   void clrAllFlags();

   /*!
    * @brief Print the bitvector data to the specified output stream.
    */
   virtual void printClassData( std::ostream &os = tbox::plog ) const;

private:
   static int s_bits_per_long;

   int d_num_bitvector_longs;
   unsigned long* d_vector;
};

}
}

#ifndef DEBUG_NO_INLINE
#include "ComponentSelector.I"
#endif
#endif
