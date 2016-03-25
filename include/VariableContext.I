//
// File:        VariableContext.I	
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
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
const string& VariableContext::getName() const
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
