/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Parser that reads the input database grammar
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
Parser*
Parser::getParser()
{
   return s_default_parser;
}

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<Database>&
Parser::getScope()
{
   return d_scope_stack.getFirstItem();
}

SAMRAI_INLINE_KEYWORD
void
Parser::enterScope(
   const std::string& name)
{
   d_scope_stack.addItem(d_scope_stack.getFirstItem()->putDatabase(name));
}

SAMRAI_INLINE_KEYWORD
void
Parser::leaveScope()
{
   d_scope_stack.removeFirstItem();
}

SAMRAI_INLINE_KEYWORD
int
Parser::getNumberErrors() const
{
   return d_errors;
}

SAMRAI_INLINE_KEYWORD
int
Parser::getNumberWarnings() const
{
   return d_warnings;
}

SAMRAI_INLINE_KEYWORD
void
Parser::setLine(
   const std::string& line)
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();
   pd.d_linebuffer = line;
}

}
}
