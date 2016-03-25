//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mesh/clustering/HistogramBox.I $
// Package:     SAMRAI mesh generation
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
