//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/variables/VariableContext.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Simple integer id and namestring variable context
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace hier {

inline
int VariableContext::getCurrentMaximumInstanceNumber()
{
   return(s_instance_counter);
}

inline
int VariableContext::getIndex() const
{
   return(d_index);    
}

inline
const std::string& VariableContext::getName() const
{
   return(d_name);
}

inline
bool VariableContext::operator==(const VariableContext& other) const
{
   return( (d_index == other.d_index) );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
