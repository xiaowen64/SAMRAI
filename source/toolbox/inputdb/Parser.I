//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/inputdb/Parser.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Parser that reads the input database grammar
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {

inline
Parser *Parser::getParser()
{
   return(s_default_parser);
}

inline
Pointer<Database>& Parser::getScope()
{
   return(d_scope_stack.getFirstItem());
}

inline
void Parser::enterScope(const std::string& name)
{
   d_scope_stack.addItem(d_scope_stack.getFirstItem()->putDatabase(name));
}

inline
void Parser::leaveScope()
{
   d_scope_stack.removeFirstItem();
}

inline
int Parser::getNumberErrors() const
{
   return(d_errors);
}

inline
int Parser::getNumberWarnings() const
{
   return(d_warnings);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

