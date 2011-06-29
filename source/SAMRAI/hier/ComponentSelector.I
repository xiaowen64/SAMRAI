/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple bit vector.  The size of the bit vector depends on
 *                The number of PatchData entries.  It is implemented as a 
 *                vector of bitsets to take advantage of the speed of 
 *                processing elements in a bitset.
 *
 ************************************************************************/
#include "SAMRAI/tbox/SAMRAIManager.h"
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
namespace hier {

SAMRAI_INLINE_KEYWORD
ComponentSelector::ComponentSelector(
   const bool flag) : d_max_bit_index(-1)
{
   int pd_entries = tbox::SAMRAIManager::getMaxNumberPatchDataEntries();
   int num_bitset_elements = pd_entries / C_BITSET_SIZE;
   int num_mod = pd_entries % C_BITSET_SIZE;

   if (num_mod != 0) {
      num_bitset_elements++;
   }  

   std::bitset<C_BITSET_SIZE> l_bits;
   d_bit_vector.resize(num_bitset_elements, l_bits);   

   if (flag) {
      // use the bitset "set" operation to set each vector element's
      // bitset values to "true".
      for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
         d_bit_vector[vi].set();
      }
      d_max_bit_index =
         (static_cast<int>(d_bit_vector.size()) * C_BITSET_SIZE) - 1;
   }
}

SAMRAI_INLINE_KEYWORD
ComponentSelector::ComponentSelector(
   const ComponentSelector& flags):
   tbox::DescribedClass()
{
   d_bit_vector = flags.d_bit_vector;
   d_max_bit_index = flags.d_max_bit_index;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector::~ComponentSelector()
{
}

SAMRAI_INLINE_KEYWORD
int ComponentSelector::getSize() const
{
   return static_cast<int>(d_bit_vector.size()) * C_BITSET_SIZE;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector&
ComponentSelector::operator = (
   const ComponentSelector& flags)
{
   d_bit_vector = flags.d_bit_vector;
   d_max_bit_index = flags.d_max_bit_index;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool
ComponentSelector::operator == (
   const ComponentSelector& flags) const
{
   return d_bit_vector == flags.d_bit_vector;
}

SAMRAI_INLINE_KEYWORD
bool
ComponentSelector::operator != (
   const ComponentSelector& flags) const
{
   return d_bit_vector != flags.d_bit_vector;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector
ComponentSelector::operator | (
   const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      tmp.d_bit_vector[vi] = d_bit_vector[vi] | flags.d_bit_vector[vi];
   }
   tmp.d_max_bit_index = 
      tbox::MathUtilities<int>::Max(
         d_max_bit_index, flags.d_max_bit_index);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector
ComponentSelector::operator& (
   const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      tmp.d_bit_vector[vi] = d_bit_vector[vi] & flags.d_bit_vector[vi];
   }
   tmp.d_max_bit_index = _findMaxIndex(tmp.d_bit_vector);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector
ComponentSelector::operator ! () const
{
   ComponentSelector tmp;
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      tmp.d_bit_vector[vi] = ~(d_bit_vector[vi]);
   }
   tmp.d_max_bit_index = _findMaxIndex(tmp.d_bit_vector);
   return tmp;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector&
ComponentSelector::operator |= (
   const ComponentSelector& flags)
{
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      d_bit_vector[vi] |= flags.d_bit_vector[vi];
   }
   d_max_bit_index = 
      tbox::MathUtilities<int>::Max(
         d_max_bit_index, flags.d_max_bit_index);
   return *this;
}

SAMRAI_INLINE_KEYWORD
ComponentSelector&
ComponentSelector::operator &= (
   const ComponentSelector& flags)
{
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      d_bit_vector[vi] &= flags.d_bit_vector[vi];
   }
   d_max_bit_index = _findMaxIndex(d_bit_vector);
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool ComponentSelector::isSet(
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < getSize()));
   return d_bit_vector[_index(i)].test(_element(i));
}

SAMRAI_INLINE_KEYWORD
void ComponentSelector::setFlag(
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < getSize()));
   d_bit_vector[_index(i)].set(_element(i));
}

SAMRAI_INLINE_KEYWORD
void ComponentSelector::clrFlag(
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < getSize()));
   d_bit_vector[_index(i)].reset(_element(i));
   d_max_bit_index = _findMaxIndex(d_bit_vector);
}

SAMRAI_INLINE_KEYWORD
void ComponentSelector::setAllFlags()
{
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      (d_bit_vector[vi]).set();
   }
   d_max_bit_index = C_BITSET_SIZE * static_cast<int>(d_bit_vector.size()) - 1;
}

SAMRAI_INLINE_KEYWORD
void ComponentSelector::clrAllFlags()
{
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      (d_bit_vector[vi]).reset();
   }
   d_max_bit_index = -1;
}

SAMRAI_INLINE_KEYWORD
int ComponentSelector::getMaxIndex() const
{
   return d_max_bit_index;
}

SAMRAI_INLINE_KEYWORD
int ComponentSelector::_index(const int i) const
{
   return i / C_BITSET_SIZE; 
}

SAMRAI_INLINE_KEYWORD
int ComponentSelector::_element(const int i) const
{
   return i % C_BITSET_SIZE; 
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
