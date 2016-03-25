//
// File:	Variable.I
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
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
const string& Variable<DIM>::getName() const
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
