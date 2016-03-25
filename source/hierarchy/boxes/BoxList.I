//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxList.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	A list of boxes with basic domain calculus operations
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {


template<int DIM> inline
BoxList<DIM>::BoxList() : tbox::List< Box<DIM> >()
{
}

template<int DIM> inline
int BoxList<DIM>::getNumberOfBoxes() const
{
   return(this->getNumberItems());
}

template<int DIM> inline
void BoxList<DIM>::unionBoxes(const Box<DIM>& box)
{
   if (!box.empty()) {
      appendItem(box);
   }
}

template<int DIM> inline
void BoxList<DIM>::unionBoxes(const BoxList<DIM>& boxes)
{
   copyItems(boxes);
}

template<int DIM> inline
BoxList<DIM>::~BoxList()
{
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

