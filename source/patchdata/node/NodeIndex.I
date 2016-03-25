//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/node/NodeIndex.I $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Index for node centered patch data types
//

#ifdef DEBUG_CHECK_ASSERTIONS
#include "tbox/Utilities.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace pdat {

template<int DIM> inline
NodeIndex<DIM>::NodeIndex()
{
   setOffsets();
}

template<int DIM> inline
NodeIndex<DIM>::NodeIndex(const hier::Index<DIM>& rhs, const Corner corner)
{
   setOffsets();
   hier::IntVector<DIM>::operator=(rhs+s_offsets[(int) corner]);
}

template<int DIM> inline
NodeIndex<DIM>::NodeIndex(const hier::Index<DIM>& rhs,
                                 const hier::IntVector<DIM>& corner)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < DIM; i++) {
      TBOX_ASSERT(corner(i) == 0 || corner(i) == 1);
   }
#endif
   setOffsets();
   hier::IntVector<DIM>::operator=(rhs+corner);
}

template<int DIM> inline
NodeIndex<DIM>::NodeIndex(const NodeIndex<DIM>& rhs) : hier::Index<DIM>(rhs)
{
   setOffsets();
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator=(const NodeIndex<DIM>& rhs)
{
   hier::Index<DIM>::operator=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM>::~NodeIndex()
{
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator+=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM> NodeIndex<DIM>::operator+(const hier::IntVector<DIM>& rhs) const
{
   NodeIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator+=(const int rhs)
{
   hier::Index<DIM>::operator+=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM> NodeIndex<DIM>::operator+(const int rhs) const
{
   NodeIndex<DIM> tmp = *this;
   tmp += rhs;
   return(tmp);
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator-=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM> NodeIndex<DIM>::operator-(const hier::IntVector<DIM>& rhs) const
{
   NodeIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator-=(const int rhs)
{
   hier::Index<DIM>::operator-=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM> NodeIndex<DIM>::operator-(const int rhs) const
{
   NodeIndex<DIM> tmp = *this;
   tmp -= rhs;
   return(tmp);
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator*=(const hier::IntVector<DIM>& rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM> NodeIndex<DIM>::operator*(const hier::IntVector<DIM>& rhs) const
{
   NodeIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
NodeIndex<DIM>& NodeIndex<DIM>::operator*=(const int rhs)
{
   hier::Index<DIM>::operator*=(rhs);
   return(*this);
}

template<int DIM> inline
NodeIndex<DIM> NodeIndex<DIM>::operator*(const int rhs) const
{
   NodeIndex<DIM> tmp = *this;
   tmp *= rhs;
   return(tmp);
}

template<int DIM> inline
bool NodeIndex<DIM>::operator==(const NodeIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator==(rhs));
}

template<int DIM> inline
bool NodeIndex<DIM>::operator!=(const NodeIndex<DIM>& rhs) const
{
   return(((hier::Index<DIM> *) this)->operator!=(rhs));
}

template<int DIM> inline
void NodeIndex<DIM>::setOffsets()
{

   if (!s_offsets_are_set) {
      for (int i = 0; i < (1 << DIM); i++) {
         hier::IntVector<DIM> offset(0);

         offset(0) = i%2;
         for (int j = 1; j < DIM; j++) {
            offset(j) = (i/(1 << j))%2;
         }
         s_offsets[i] = offset;
      }
      s_offsets_are_set = true;
   }
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
