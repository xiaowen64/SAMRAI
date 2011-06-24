/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple bit vector of a fixed length (128 bits) 
 *
 ************************************************************************/

#ifndef included_hier_ComponentSelector_C
#define included_hier_ComponentSelector_C

#include "SAMRAI/hier/ComponentSelector.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/ComponentSelector.I"
#endif

namespace SAMRAI {
namespace hier {

bool
ComponentSelector::any() const {
   std::vector<std::bitset<C_BITSET_SIZE> >::const_iterator iter;
   bool set = false;
   for (iter = d_bit_vector.begin(); iter != d_bit_vector.end() && !set;
        ++iter) {
      set = iter->any();
   }
   return set;
}

bool
ComponentSelector::none() const {
   return !any();
}

int ComponentSelector::_findMaxIndex(
      const std::vector<std::bitset<C_BITSET_SIZE> >& bits) const
{
   bool bits_set = false;
   int max_index = -1;
   for (size_t i = 0; i < bits.size() && !bits_set; ++i)  {
      bits_set |= bits[i].any();
   }

   if (bits_set) {
      int j = C_BITSET_SIZE - 1;
      while (!bits[_index(j)].test(_element(j))) {
         --j;
      }
      max_index = j;
   }
   return max_index;
}


void ComponentSelector::printClassData(
   std::ostream& os) const
{
   int i;
   const int number_of_bits = getSize();
   for (i = 0; i < number_of_bits; ++i) {
      os << " | Bit " << i << " = " << isSet(i);
   }
   os << "|\n";
}

}
}

#endif
