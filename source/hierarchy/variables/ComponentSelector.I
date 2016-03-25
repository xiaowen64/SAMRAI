//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/ComponentSelector.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Simple bit vector of a fixed length (128 bits)
//

#include "tbox/SAMRAIManager.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {

inline
ComponentSelector::ComponentSelector(const bool flag)
{
   d_num_bitvector_longs =
      tbox::SAMRAIManager::getMaxNumberPatchDataEntries() / s_bits_per_long;
   int num_mod = 
      tbox::SAMRAIManager::getMaxNumberPatchDataEntries() % s_bits_per_long;
   if (num_mod != 0) {
      d_num_bitvector_longs += 1; 
   }
   d_vector = new unsigned long[d_num_bitvector_longs];

   unsigned long v = (flag ? ~((unsigned long) 0) : ((unsigned long) 0));
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] = v;
   }
}

inline
ComponentSelector::ComponentSelector(
   const ComponentSelector& flags)
{
   d_num_bitvector_longs = flags.d_num_bitvector_longs;
   d_vector = new unsigned long[d_num_bitvector_longs];

   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] = flags.d_vector[i];
   }
}

inline
ComponentSelector::~ComponentSelector()
{
   if (d_vector) {
      d_num_bitvector_longs = 0;
      delete [] d_vector;
      d_vector = (unsigned long*)0;
   }
}

inline
int ComponentSelector::getSize() const
{
   return(d_num_bitvector_longs*s_bits_per_long);
}

inline
ComponentSelector&
ComponentSelector::operator=(const ComponentSelector& flags)
{
   d_num_bitvector_longs = flags.d_num_bitvector_longs;
   d_vector = new unsigned long[d_num_bitvector_longs];

   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] = flags.d_vector[i];
   }
   return(*this);
}

inline
bool
ComponentSelector::operator==(const ComponentSelector& flags) const
{
   bool ret_val = true;
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      ret_val = ret_val && (d_vector[i] == flags.d_vector[i]);
   }
   return(ret_val);
}

inline
bool
ComponentSelector::operator!=(const ComponentSelector& flags) const
{
   bool ret_val = false;
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      ret_val = ret_val || (d_vector[i] != flags.d_vector[i]);
   }
   return(ret_val);
}

inline
ComponentSelector
ComponentSelector::operator|(const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      tmp.d_vector[i] = d_vector[i] | flags.d_vector[i];
   }
   return(tmp);
}

inline
ComponentSelector
ComponentSelector::operator&(const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      tmp.d_vector[i] = d_vector[i] & flags.d_vector[i];
   }
   return(tmp);
}

inline
ComponentSelector
ComponentSelector::operator!() const
{
   ComponentSelector tmp;
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      tmp.d_vector[i] = ~(d_vector[i]);
   }
   return(tmp);
}

inline
ComponentSelector
ComponentSelector::andNot(const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      tmp.d_vector[i] = d_vector[i] & ~(flags.d_vector[i]);
   }
   return(tmp);
}

inline
ComponentSelector&
ComponentSelector::operator|=(const ComponentSelector& flags)
{
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] |= flags.d_vector[i];
   }
   return(*this);
}

inline
ComponentSelector&
ComponentSelector::operator&=(const ComponentSelector& flags)
{
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] &= flags.d_vector[i];
   }
   return(*this);
}

inline
bool ComponentSelector::isSet(const int i) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < getSize()));
#endif
   const int index = ((unsigned) i)/s_bits_per_long;
   const int bitps = ((unsigned) i) & (s_bits_per_long - 1);
   return((d_vector[index] & (((unsigned long) 1) << bitps)) != 0);
}

inline
void ComponentSelector::setFlag(const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < getSize()));
#endif
   const int index  = ((unsigned) i)/s_bits_per_long;
   const int bitps  = ((unsigned) i) & (s_bits_per_long - 1);
   d_vector[index] |= (((unsigned long) 1) << bitps);
}

inline
void ComponentSelector::clrFlag(const int i)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((i >= 0) && (i < getSize()));
#endif
   const int index  = ((unsigned) i)/s_bits_per_long;
   const int bitps  = ((unsigned) i) & (s_bits_per_long - 1);
   d_vector[index] &= ~(((unsigned long) 1) << bitps);
}

inline
void ComponentSelector::setAllFlags()
{
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] = ~((unsigned long) 0);
   }
}

inline
void ComponentSelector::clrAllFlags()
{
   for (int i = 0; i < d_num_bitvector_longs; i++) {
      d_vector[i] = (unsigned long) 0;
   }
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
