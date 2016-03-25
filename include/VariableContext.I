//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/VariableContext.I $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
