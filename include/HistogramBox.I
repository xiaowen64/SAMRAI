//
// File:        HistogramBox.I
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Histogram class for computing tagged cell signatures.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace mesh {

template<int DIM> inline 
const hier::Box<DIM>& HistogramBox<DIM>::getBox() const
{
   return(d_box);
}

template<int DIM> inline
tbox::Array<int>& HistogramBox<DIM>::histogram(const int id)
{
   return(d_histogram[id]);
}

template<int DIM> inline 
const tbox::Array<int>& HistogramBox<DIM>::histogram(const int id) const
{
   return(d_histogram[id]);
}

template<int DIM> inline 
int& HistogramBox<DIM>::histogramElt(const int id, const int ic)
{
   return( d_histogram[id][ic-(d_box.lower()(id))] );
}

template<int DIM> inline
const int& HistogramBox<DIM>::histogramElt(const int id, const int ic) const
{
   return( d_histogram[id][ic-(d_box.lower()(id))] );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
