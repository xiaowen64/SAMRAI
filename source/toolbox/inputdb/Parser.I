//
// File:	Parser.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
void Parser::enterScope(const string& name)
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

