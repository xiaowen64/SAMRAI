//
// File:	BoxList.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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

