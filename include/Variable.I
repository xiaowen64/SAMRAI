//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/Variable.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Base class for application-level variables
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace hier {

template<int DIM> inline
int Variable<DIM>::getCurrentMaximumInstanceNumber()
{
   return(s_instance_counter);
}

template<int DIM> inline
int Variable<DIM>::getInstanceIdentifier() const
{
   return(d_instance);
}

template<int DIM> inline
const std::string& Variable<DIM>::getName() const
{
   return(d_name);
}

template<int DIM> inline
void Variable<DIM>::setPatchDataFactory(
   tbox::Pointer< PatchDataFactory<DIM> > factory)
{
   d_factory = factory;
}

template<int DIM> inline
tbox::Pointer< PatchDataFactory<DIM> >
Variable<DIM>::getPatchDataFactory() const
{
   return(d_factory);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
